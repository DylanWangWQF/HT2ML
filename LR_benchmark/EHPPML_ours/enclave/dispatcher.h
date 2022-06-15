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

    long MatrixDim = 16;


  public:
    // ecall_dispatcher(HTparams* params);
    ecall_dispatcher();
    int enclave_init(uint8_t* hecontext, size_t context_len);
    int MatrixOperation(uint8_t* ectxt, size_t ectxt_len, uint8_t** octxt, size_t* octxt_len);
    void close();
  private:
    int** matMul(int** A, int** B);
    void LUP_Descomposition(int**& A, int**& L, int**& U, int*& P);
};
