#pragma once

#include <cstring>
#include <sstream>
#include <math.h>
#include <helib/helib.h>
#include <openenclave/enclave.h>

#define MatrixDim 16

// #define numAttr 2
// #define numAttr 4
#define numAttr 6
// #define numAttr 8
// #define numAttr 16

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

    //int MatrixDim = 16;
    // This param should be passed to enclave. Here, for simplicity, we just hard-coded it.
    // int numAttr = 6;


  public:
    // ecall_dispatcher(HTparams* params);
    ecall_dispatcher();
    int enclave_init(uint8_t* hecontext, size_t context_len);
    int MatrixOperation(uint8_t* ectxt, size_t ectxt_len, uint8_t** octxt, size_t* octxt_len);
    void close();
  private:
    void getCofactor(int A[numAttr][numAttr], int temp[numAttr][numAttr], int p, int q, int n);
    int determinant(int A[numAttr][numAttr], int n);
    void adjoint(int A[numAttr][numAttr],int adj[numAttr][numAttr]);
    bool inverse(int A[numAttr][numAttr], int inverse[numAttr][numAttr]);
    // void display(int A[MatrixDim][MatrixDim]);
    // void mul(int A[MatrixDim][MatrixDim], int B[MatrixDim][MatrixDim], int C[MatrixDim][MatrixDim]);
};
