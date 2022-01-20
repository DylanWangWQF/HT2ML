#include <string.h>
#include <NTL/RR.h>
#include <vector>
#include <chrono>

#include "dispatcher.h"
#include "trace.h"

ecall_dispatcher::ecall_dispatcher() {}

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
    TRACE_ENCLAVE("Enclave: receive the HE context and keys");
    stringstream ess;
    string e_hecontext = string(hecontext, hecontext + context_len);
    ess << e_hecontext;
    TRACE_ENCLAVE("Enclave: reconstruct the context: ");
    e_context = Context::readPtrFrom(ess);
    e_context->printout();

    TRACE_ENCLAVE("Enclave: reconstruct the SecKey!");
    activeSecKey = make_unique<SecKey>(SecKey::readFrom(ess, *e_context));

    TRACE_ENCLAVE("Enclave: reconstruct the PubKey!");
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

int ecall_dispatcher::multipleCtxtsTransform(uint8_t* ectxt, size_t ectxt_len, size_t num_ectxt, uint8_t** octxt, size_t* octxt_len)
{
    stringstream css;
    vector<Ctxt> ReceivedCtxts(num_ectxt, Ctxt(*activePubKey));

    TRACE_ENCLAVE("Enclave: Receiving Ctxt start!");
    string etemp = string(ectxt, ectxt + ectxt_len);
    css << etemp;
    for (size_t i = 0; i < num_ectxt; ++i)
    {
        ReceivedCtxts[i].Ctxt::read(css);
        RefreshRmat(ReceivedCtxts[i]);
    }

    css.str(std::string());
    css.clear();
    TRACE_ENCLAVE("Enclave: Transforming Ctxt start!");
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

    TRACE_ENCLAVE("Enclave: Receiving Ctxt start!");
    string etemp = string(ectxt, ectxt + ectxt_len);
    css << etemp;
    ReceivedCtxts.Ctxt::read(css);
    RefreshRmat(ReceivedCtxts);

    css.str(std::string());
    css.clear();
    TRACE_ENCLAVE("Enclave: Transforming Ctxt start!");
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

void ecall_dispatcher::close()
{
    TRACE_ENCLAVE("Enclave: release context and HTparams!");
    // release context and keys
    delete e_context;
    // activePubKey.reset();
    // activeSecKey.reset();
    TRACE_ENCLAVE("ecall_dispatcher::close");
}