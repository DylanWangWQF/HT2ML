//
//  CKKSmatrix.h
//  ePPDSC
//
//  Created by Qifan Wang on 30/12/21.
//

#ifndef CKKSmatrix_h
#define CKKSmatrix_h

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <sys/time.h>

#include <NTL/BasicThreadPool.h>
#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include "NTL/RR.h"
#include <NTL/ZZX.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>

#include <helib/helib.h>

using namespace std;
using namespace NTL;
using namespace helib;

typedef struct ckksMatpar{
    //! real dimension of data
    long nrows;
    long ncols;
    
    //! power-of-two integer larger than dim
    long dim;
    long dim1;
    //!it can be used for denoting the number of rows (or columns) in rectangular matrix multiplication
    long subdim;
    
    long logdim;
    long neededslots; //! the number of slots (Power of Two)
    long sqrdim;     //! sqrt(d)
    long nbatching;  //! number of matrices in a single ciphertext (default = 1)
    
    long m, bits, prec, c;
    
}ckksMatpar;

// assign parameters
struct ckksParams
{
  const long m, bits, prec, c;
    ckksParams(long _m,
               long _bits,
               long _prec,
               long _c) :
      m(_m), bits(_bits), prec(_prec), c(_c) {}
    ckksParams(const ckksParams& other) :
    ckksParams(other.m,
               other.bits,
               other.prec,
               other.c) {}
  bool operator!=(ckksParams& other) const { return !(*this == other); }
  bool operator==(ckksParams& other) const
  {
    return m == other.m && bits == other.bits && prec == other.prec && c == other.c;
  }
};

struct ckksContextAndKeys
{
  const ckksParams params;
  Context context;
  SecKey secretKey;
  const PubKey& publicKey;
  const EncryptedArray& ea;

    ckksContextAndKeys(ckksParams& _params) :
        params(_params),
        context(ContextBuilder<CKKS>()
                    .m(params.m)
                    .bits(params.bits)
                    .precision(params.prec)
                    .c(params.c)
                    .build()),
        secretKey(context),
        publicKey((secretKey.GenSecKey(),
                   addSome1DMatrices(secretKey),
                   secretKey)),
        ea(context.getEA()){}
};

struct ckksMeta
{
    std::unique_ptr<ckksContextAndKeys> data;
    ckksMeta& operator()(ckksParams& params)
    {
        // Only change if nullptr or different.
        if (data == nullptr || data->params != params)
            data = std::make_unique<ckksContextAndKeys>(params);
        return *this;
    }
};

void readckksMatpar(ckksMatpar& ckksmatpar, long nrows, long ncols, const long subdim = 0, const long nbatching = 1);

void printRvector(vec_RR& vec, const long print_size = 0);

void printRmatrix(Mat<RR>& mat, const long print_size = 0);

RR getError(mat_RR Amat, mat_RR Bmat, const long nrows = 0, const long ncols = 0);

class CKKSmatrix{
public:
    ckksMatpar& ckksmatpar;
    ckksMeta& ckksmeta;
    CKKSmatrix(ckksMatpar& ckksmatpar, ckksMeta& ckksmeta) : ckksmatpar(ckksmatpar), ckksmeta(ckksmeta) {}
    
    // Encrypt and Decrypt matrix
    void encryptRmat(Ctxt& ctxt, mat_RR& mat);
    void decryptRmat(mat_RR& mat, Ctxt& ctxt);
    void encryptParallelRmat(Ctxt& ctxt, mat_RR*& mat, long nbatching);
    void decryptParallelRmat(mat_RR*& mat, Ctxt& ctxt);
    
    // Rotate vectors
    void msgleftRotate(vector<complex<double>>& res, vector<complex<double>> vals, long dim, long nrot);
    void msgrightRotate(vector<complex<double>>& res, vector<complex<double>> vals, long dim, long nrot);
    void msgleftRotateAndEqual(vector<complex<double>>& vals, long dim, long nrot);
    void msgrightRotateAndEqual(vector<complex<double>>& vals, long dim, long nrot);
    
    // Transposition
    void genTransPoly(vector<EncodedPtxt>& transpoly);
    void transpose(Ctxt& res, Ctxt& ctxt, vector<EncodedPtxt>& transpoly);
    void genTransPoly_Parallel(vector<EncodedPtxt>& transpoly);
    void transpose_Parallel(Ctxt& res, Ctxt& ctxt, vector<EncodedPtxt>& transpoly);
    
    // Shift
    void genShiftPoly(vector<EncodedPtxt>& shiftpoly, const long num = 0);
    void shiftBycols(Ctxt& res, Ctxt& ctxt, long k, vector<EncodedPtxt>& shiftpoly);
    void genShiftPoly_Parallel(vector<EncodedPtxt>& shiftpoly);
    void shiftBycols_Parallel(Ctxt& res, Ctxt& ctxt, long k, vector<EncodedPtxt>& shiftpoly);
    
    // multiplication preparation
    void genMultPoly(vector<vector<EncodedPtxt>>& Initpoly);
    void genMultPoly_Parallel(vector<vector<EncodedPtxt>>& Initpoly);
    void genMultBPoly(vector<EncodedPtxt>& Initpoly);
    
    void genInitCtxt(Ctxt& resA, Ctxt& resB, Ctxt& Actxt, Ctxt& Bctxt, vector<vector<EncodedPtxt>>& Initpoly);
    void genInitCtxt_Parallel(Ctxt& resA, Ctxt& resB, Ctxt& Actxt, Ctxt& Bctxt, vector<vector<EncodedPtxt>>& Initpoly);
    void genInitActxt(vector<Ctxt>& Actxts, Mat<RR>& mat);
    void genInitBctxt(Ctxt& resB, Ctxt& Bctxt, vector<EncodedPtxt>& Initpoly);
    void genInitRecActxt(vector<Ctxt>& Actxts, Mat<RR>& mat);
    
    // multiplication
    void HEmatmul_Hadamard(Ctxt& res, vector<Ctxt> Actxts, vector<Ctxt> Bctxts, long num);
    void HEmatmul(Ctxt& res, Ctxt& Actxt, Ctxt& Bctxt, vector<vector<EncodedPtxt>>& Initpoly, vector<EncodedPtxt>& shiftpoly);
    void HEmatmul_Parallel(Ctxt& res, Ctxt& Actxt, Ctxt& Bctxt, vector<vector<EncodedPtxt>>& Initpoly, vector<EncodedPtxt>& shiftpoly);
    void HErmatmul(Ctxt& res, Ctxt& Actxt, Ctxt& Bctxt, vector<vector<EncodedPtxt>>& Initpoly, vector<EncodedPtxt>& shiftpoly);
    void HEmatmul_preprocessing(Ctxt& res, vector<Ctxt>& Actxts, Ctxt& Bctxt, vector<EncodedPtxt>& Initpoly);
    void HErmatmul_preprocessing(Ctxt& res, vector<Ctxt>& Actxts, Ctxt& Bctxt, vector<EncodedPtxt>& Initpoly);
};

#endif /* CKKSmatrix_h */
