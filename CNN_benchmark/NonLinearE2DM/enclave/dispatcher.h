#pragma once

#include <cstring>
#include <sstream>
#include <math.h>
#include <helib/helib.h>
#include <openenclave/enclave.h>

using namespace std;
using namespace helib;
using namespace NTL;

class ecall_dispatcher
{
  private:
    // HE context pointer
    Context* e_context;
    // HE secret key
    unique_ptr<SecKey> activeSecKey;
    // HE public key
    unique_ptr<PubKey> activePubKey;


  public:
    // ecall_dispatcher(HTparams* params);
    ecall_dispatcher();
    int enclave_init(uint8_t* hecontext, size_t context_len);
    int nonLinearLayer(uint8_t* ectxt, size_t ectxt_len, size_t num_ectxt, uint8_t** octxt, size_t* octxt_len);
    void close();
  private:
    double sigmoid(double x);
    void ProcessNonLinear(vector<vector<double>>& pooling_output, Ctxt& ctxt);
};
