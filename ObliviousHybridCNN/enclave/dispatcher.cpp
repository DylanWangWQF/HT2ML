#include <string.h>
#include <NTL/matrix.h>
#include <vector>
#include <chrono>

#include "include/obl/obl_primitives.h"
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

    // TRACE_ENCLAVE("Enclave: receive the HE params");
    stringstream ess;
    string e_hecontext = string(hecontext, hecontext + context_len);
    // cout << "Enclave: check hecontext: " << hecontext << endl;
    // cout << "Enclave: check context_len, containing context, sk and pk: " << context_len << endl;
    ess << e_hecontext;
    // TRACE_ENCLAVE("Enclave: readfrom the string and reconstruct the context: ");
    e_context = Context::readPtrFrom(ess);
    // TRACE_ENCLAVE("Enclave: Setup Seckey by reading from ss!");
    activeSecKey = make_unique<SecKey>(SecKey::readFrom(ess, *e_context));
    // TRACE_ENCLAVE("Enclave: Setup PubKey by reading from ss!");
    activePubKey = make_unique<PubKey>(PubKey::readFrom(ess, *e_context));

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Enclave: Total Setup Time for HE inside enclave = " << timeElapsed << " s" << endl;
    cout << "------------------------------------------------------------------------" << endl;

exit:
    TRACE_ENCLAVE("Enclave: free memory in enclave_init()!");
    e_hecontext.shrink_to_fit();
    return ret;
}

int ecall_dispatcher::nonLinearLayer(uint8_t* ectxt, size_t ectxt_len, uint8_t** octxt, size_t* octxt_len)
{
    // test obl primitives
    // cout << endl << "Obl greater from (5.9, 0.0) = " << ObliviousChoose((ObliviousGreater(5.9, 0.0)), 5.9, 0.0) << endl;
    // vector<double> dylan;
    // dylan.push_back(2.9);
    // dylan.push_back(0.7);
    // dylan.push_back(4.9);
    // dylan.push_back(1.2);
    // ObliviousSort(dylan.begin(), dylan.end());
    // for (int i = 0; i < dylan.size(); i++)
    // {
    //     cout << dylan[i] << " ";
    // }
    // cout << endl << endl;

    auto start= chrono::steady_clock::now();
    stringstream css;
    // vector<Ctxt> ReceivedCtxts(num_ectxt, Ctxt(*activePubKey)); // num_ectxt = num_channel

    // TRACE_ENCLAVE("Enclave: Receiving Ctxt start!");
    string etemp = string(ectxt, ectxt + ectxt_len);
    css << etemp;
    vector<vector<double>> pooling_output(64); // 64 images, 6 * 12 * 12 for each image
    for (int m = 0; m < 6; m++)
    {
        // vector<vector<double>> mat(64,vector<double>(576));
        vector<vector<double>> mat(64);
        for (size_t n = 0; n < 9; n++)
        {
            Ctxt ctemp(*activePubKey);
            ctemp.Ctxt::read(css);
            vector<complex<double>> cmsg;
            PtxtArray pa1(*e_context);
            pa1.decryptComplex(ctemp, *activeSecKey);
            pa1.store(cmsg);

            for(long i = 0; i < 64; i++){
                for(long j = 0; j < 64; j++){
                    mat[i].push_back(cmsg[64 * j + i].real());
                }
            }
        }
        ProcessNonLinear(pooling_output, mat);
    }

    // output 14 of 64 * 64
    css.str(std::string());
    css.clear();
    // TRACE_ENCLAVE("Enclave: Transforming Ctxt start!");
    for (int k = 0; k < 14; k++)
    {
        vector<complex<double>> cmsg(64 * 64);
        for (int i = 0; i < 64; i++)
        {
            for(int j = (64 * k); j < (64 * (k + 1)); j++)
            {
                if (j < 864) // 6 * 12 * 12
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

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Enclave: Non Linear Layer Runtime inside enclave = " << timeElapsed << " s" << endl;
    cout << "------------------------------------------------------------------------" << endl;

    return 0;
}

double ecall_dispatcher::relu(double x)
{
    return ObliviousChoose((ObliviousGreater(x, 0.0)), x, 0.0);
}

void ecall_dispatcher::ProcessNonLinear(vector<vector<double>>& pooling_output, vector<vector<double>>& mat)
{
    for(long i = 0; i < 64; i++){
        // i-th image
        for (int m = 0; m < 12; m++)
        {
            for (int n = 0; n < 12; n++)
            {
                vector<double> temp_sort;
                temp_sort.push_back(relu(mat[i][(m * 48 + n * 2)]));
                temp_sort.push_back(relu(mat[i][(m * 48 + n * 2 + 1)]));
                temp_sort.push_back(relu(mat[i][(m * 48 + n * 2 + 24)]));
                temp_sort.push_back(relu(mat[i][(m * 48 + n * 2 + 25)]));
                ObliviousSort(temp_sort.begin(), temp_sort.end());
                pooling_output[i].push_back(temp_sort[0]);
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