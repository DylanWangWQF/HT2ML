#include <string.h>
#include <NTL/matrix.h>
#include <vector>
#include <chrono>

#include "dispatcher.h"
#include "trace.h"

//HTparams* params
ecall_dispatcher::ecall_dispatcher()
{
    // e_HTparams = params;
}

int ecall_dispatcher::enclave_init(uint8_t* hecontext, size_t context_len)
{
    int ret = 0;
    // if (context_len == 0)
    // {
    //     TRACE_ENCLAVE("passed hecontext is null string, context_len =0");
    //     ret = 1;
    //     goto exit;
    // }
    auto start= chrono::steady_clock::now();

    TRACE_ENCLAVE("Enclave: receive the HE params");
    stringstream ess;
    string e_hecontext = string(hecontext, hecontext + context_len);
    // cout << "Enclave: check hecontext: " << hecontext << endl;
    // cout << "Enclave: check context_len, containing context, sk and pk: " << context_len << endl;
    ess << e_hecontext;
    TRACE_ENCLAVE("Enclave: readfrom the string and reconstruct the context: ");
    e_context = Context::readPtrFrom(ess);
    TRACE_ENCLAVE("Enclave: Setup Seckey by reading from ss!");
    activeSecKey = make_unique<SecKey>(SecKey::readFrom(ess, *e_context));
    TRACE_ENCLAVE("Enclave: Setup PubKey by reading from ss!");
    activePubKey = make_unique<PubKey>(PubKey::readFrom(ess, *e_context));

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Enclave: Total Setup Time for HE = " << timeElapsed << " s" << endl;
    cout << "------------------------------------------------------------------------" << endl;

exit:
    TRACE_ENCLAVE("Enclave: free memory in enclave_init()!");
    e_hecontext.shrink_to_fit();
    return ret;
}

int ecall_dispatcher::nonLinearLayer(uint8_t* ectxt, size_t ectxt_len, size_t num_ectxt, uint8_t** octxt, size_t* octxt_len)
{
    stringstream css;
    // vector<Ctxt> ReceivedCtxts(num_ectxt, Ctxt(*activePubKey)); // num_ectxt = num_channel

    // TRACE_ENCLAVE("Enclave: Receiving Ctxt start!");
    string etemp = string(ectxt, ectxt + ectxt_len);
    css << etemp;
    vector<vector<double>> pooling_output(64); // 64 images, 6 * 4 * 4 for each image
    for (int i = 0; i < num_ectxt; ++i)
    {
        Ctxt ctemp(*activePubKey);
        ctemp.Ctxt::read(css);
        ProcessNonLinear(pooling_output, ctemp);
    }

    // output 2 of 64 * 64 (96 * 10)
    css.str(std::string());
    css.clear();
    // TRACE_ENCLAVE("Enclave: Transforming Ctxt start!");
    for (int k = 0; k < 2; k++)
    {
        vector<complex<double>> cmsg(64 * 64);
        for (int i = 0; i < 64; i++)
        {
            for(int j = (64 * k); j < (64 * (k + 1)); j++)
            {
                if (j < 96)
                {
                    cmsg[64 * (j - 64 * k) + i].real(pooling_output[i][j]);
                }
                else
                {
                    cmsg[64 * (j - 64 * k) + i].real(0.0);
                }
            }
        }
        Ctxt wtemp(*activePubKey);
        PtxtArray pa(*e_context, cmsg);
        pa.encrypt(wtemp);
        wtemp.writeTo(css);
    }

    string otemp = css.str();
    uint8_t* host_buf = (uint8_t*) oe_host_malloc(size_t(otemp.size() + 1));
    memcpy(host_buf, (uint8_t*)otemp.c_str(), otemp.size() + 1);
    *octxt = host_buf;
    *octxt_len = otemp.size();

    return 0;
}

int ecall_dispatcher::multipleCtxtsTransform(uint8_t* ectxt, size_t ectxt_len, size_t num_ectxt, uint8_t** octxt, size_t* octxt_len)
{
    stringstream css;
    vector<Ctxt> ReceivedCtxts(num_ectxt, Ctxt(*activePubKey));

    // TRACE_ENCLAVE("Enclave: Receiving Ctxt start!");
    string etemp = string(ectxt, ectxt + ectxt_len);
    css << etemp;
    for (size_t i = 0; i < num_ectxt; ++i)
    {
        ReceivedCtxts[i].Ctxt::read(css);
        RefreshRmat(ReceivedCtxts[i]);
    }

    css.str(std::string());
    css.clear();
    // TRACE_ENCLAVE("Enclave: Transforming Ctxt start!");
    for (size_t i = 0; i < num_ectxt; ++i)
    {
        ReceivedCtxts[i].writeTo(css);
    }
    string otemp = css.str();
    uint8_t* host_buf = (uint8_t*) oe_host_malloc(size_t(otemp.size() + 1));
    memcpy(host_buf, (uint8_t*)otemp.c_str(), otemp.size() + 1);
    *octxt = host_buf;
    *octxt_len = otemp.size();

    return 0;
}

int ecall_dispatcher::singleCtxtTransform(uint8_t* ectxt, size_t ectxt_len, uint8_t** octxt, size_t* octxt_len)
{
    stringstream css;
    Ctxt ReceivedCtxts(*activePubKey);

    // TRACE_ENCLAVE("Enclave: Receiving Ctxt start!");
    string etemp = string(ectxt, ectxt + ectxt_len);
    css << etemp;
    ReceivedCtxts.Ctxt::read(css);
    RefreshRmat(ReceivedCtxts);

    css.str(std::string());
    css.clear();
    // TRACE_ENCLAVE("Enclave: Transforming Ctxt start!");
    ReceivedCtxts.writeTo(css);

    string otemp = css.str();
    uint8_t* host_buf = (uint8_t*) oe_host_malloc(size_t(otemp.size() + 1));
    memcpy(host_buf, (uint8_t*)otemp.c_str(), otemp.size() + 1);
    *octxt = host_buf;
    *octxt_len = otemp.size();

    return 0;
}

void ecall_dispatcher::RefreshRmat(Ctxt& ctxt)
{
    vector<complex<double>> cmsg;
    PtxtArray pa1(*e_context);
    pa1.decryptComplex(ctxt, *activeSecKey);
    pa1.store(cmsg);

    PtxtArray pa2(*e_context, cmsg);
    pa2.encrypt(ctxt);
    return;
}

double ecall_dispatcher::sigmoid(double x)
{
    return (1 / (1 + exp(-x)));
}

void ecall_dispatcher::ProcessNonLinear(vector<vector<double>>& pooling_output, Ctxt& ctxt)
{
    vector<complex<double>> cmsg;
    PtxtArray pa1(*e_context);
    pa1.decryptComplex(ctxt, *activeSecKey);
    pa1.store(cmsg);

    vector<vector<double>> mat(64,vector<double>(64));; // 64 * 64, vector<double> => image 
    for(long i = 0; i < 64; i++){
        for(long j = 0; j < 64; j++){
            mat[i][j] = cmsg[64 * j + i].real();
        }
        // i-th image
        for (int m = 0; m < 4; m++)
        {
            for (int n = 0; n < 4; n++)
            {
                double temp = sigmoid(mat[i][(m * 16 + n * 2)]) + sigmoid(mat[i][(m * 16 + n * 2 + 1)]) + sigmoid(mat[i][(m * 16 + n * 2 + 8)]) + sigmoid(mat[i][(m * 16 + n * 2 + 9)]);
                pooling_output[i].push_back(temp / 4);
            }
        }
        
    }
}

void ecall_dispatcher::close()
{
    TRACE_ENCLAVE("Enclave: release context and HTparams!");
    // release context and keys
    delete e_context;
    TRACE_ENCLAVE("ecall_dispatcher::close");
}