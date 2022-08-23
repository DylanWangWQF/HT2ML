//
//  hematrix.h
//  cloneHElib
//
//  Created by Qifan Wang on 19/04/21.
//

#ifndef hematrix_h
#define hematrix_h

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <sys/time.h>

#include <NTL/BasicThreadPool.h>
#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>

#include "../src/helib/helib.h"

using namespace std;
using namespace NTL;
using namespace helib;

//! structure for parameters
typedef struct HEMatpar{
    long nrows;      //! real dimension of data
    long ncols;
    
    long dim;        //! power-of-two integer larger than dim
    long dim1;
    long subdim;     //! it can be used for denoting the number of rows (or columns) in rectangular matrix multiplication
    
    long logdim;
    long neededslots;     //! the number of slots (Power of Two)
    long sqrdim;     //! sqrt(d)
    long nbatching;  //! number of matrices in a single ciphertext (default = 1)
    
    long m, p, r, bits, c;
    
}HEMatpar;

// assign parameters
struct Params
{
  const long m, p, r, bits, c;
  Params(long _m,
         long _p,
         long _r,
         long _bits,
         long _c) :
      m(_m), p(_p), r(_r), bits(_bits), c(_c) {}
  Params(const Params& other) :
      Params(other.m,
             other.p,
             other.r,
             other.bits,
             other.c) {}
  bool operator!=(Params& other) const { return !(*this == other); }
  bool operator==(Params& other) const
  {
    return m == other.m && p == other.p && r == other.r && bits == other.bits && c == other.c;
  }
};

struct ContextAndKeys
{
  const Params params;
  Context context;
  SecKey secretKey;
  const PubKey& publicKey;
  const EncryptedArray& ea;

  ContextAndKeys(Params& _params) :
        params(_params),
        context(ContextBuilder<BGV>()
                    .m(params.m)
                    .p(params.p)
                    .r(params.r)
                    .bits(params.bits)
                    .c(params.c)
                    .build()),
        secretKey(context),
        publicKey((secretKey.GenSecKey(),
                   addSome1DMatrices(secretKey),
                   secretKey)),
        ea(context.getEA()){}
};

// create object for params struct
//typedef struct Meta{
//    HEMatpar HEmatpar;
//    ContextAndKeys* cak = new ContextAndKeys(HEmatpar);
//}Meta;

struct Meta
{
  std::unique_ptr<ContextAndKeys> data;
  Meta& operator()(Params& params)
  {
    // Only change if nullptr or different.
    if (data == nullptr || data->params != params)
      data = std::make_unique<ContextAndKeys>(params);
    return *this;
  }
};

void readHEMatpar(HEMatpar& HEmatpar, long nrows, long ncols, const long subdim = 0, const long nbatching = 1);

class HEmatrix{
public:
    HEMatpar& HEmatpar;
    Meta& meta;
    
    HEmatrix(HEMatpar& HEmatpar, Meta& meta) : HEmatpar(HEmatpar), meta(meta) {}
    
    /**
     @param[in] mat The input matrix
     @param[out] ctxt generate the ciphertext
    */
    void encryptZmat(Ctxt& ctxt, Mat<long>& mat);
        
    /**
     @param[in] ctxt The input ciphertext encrypting "(1 << logp) * mat"
     @param[out] mat The matrix
    */
    void decryptZmat(Mat<long>& mat, Ctxt& ctxt);
        
    /**
     @param[in] mat The multiple input matrices
     @param[in] nbatching The number of multiple matrices in a single ciphertext
     @param[out] ctxt The ciphertext encrypting multiple matrices with scale (1 << logp)
    */
    void encryptParallelZmat(Ctxt& ctxt, Mat<long>*& mat, long nbatching);
        
    /**
     @param[in] ctxt The input ciphertext encrypting "(1 << logp) * mat"
     @param[out] mat Multiple matrices
    */
    void decryptParallelZmat(Mat<long>*& mat, Ctxt& ctxt);

    /**
     @param[in] vals The input vector, vals = (vals[0],vals[1],...,vals[d-1])
     @param[in] dim The dim of of the input vector
     @param[in] nrot The number of steps to rotate by left (positive)
     @param[out] res The output vector s.t. res = (vals[nrot], vals[nrot+1],..., vals[d-1]), vals[0], vals[1],...., vals[nrot-1])
    */
    void msgleftRotate(vector<long>& res, vector<long> vals, long dim, long nrot);
    
    void msgleftRotate_Mul(vector<long>& res, vector<long> vals, long dim, long nrot);
        
    /**
     @param[in] vals The input vector, vals = (vals[0],vals[1],...,vals[d-1])
     @param[in] dim The dim of of the input vector
     @param[in] nrot The number of steps to rotate by right (positive)
     @param[out] res The output vector s.t res = (vals[d-nrot], vals[d-nrot+1],..., vals[d-1]), vals[0], vals[1],...., vals[d-nrot-1])
    */
    void msgrightRotate(vector<long>& res, vector<long> vals, long dim, long nrot);
    
    void msgrightRotate_Mul(vector<long>& res, vector<long> vals, long dim, long nrot);
        
    /**
     @param[in] vals The input vector to rotate
     @param[in] dim The dim of of the input vector
     @param[in] nrot The number of steps to rotate by left (positive)
    */
    void msgleftRotateAndEqual(vector<long>& vals, long dim, long nrot);
    
    /**
    @param[in] vals The input vector to rotate
    @param[in] dim The dim of of the input vector
    @param[in] nrot The number of steps to rotate by right (positive)
    */
    void msgrightRotateAndEqual(vector<long>& vals, long dim, long nrot);
    
    /**
    @param[out] transpoly The polynomials needed for transposition
    */
    void genTransPoly(zzX*& transpoly);
        
    /**
    @param[in] ctxt The input ciphertext
    @param[in] transpoly The polynomials
    @param[out] res The output ciphertext that encrypts the transpose of the corresponding input matrix consumes a constant-multiplication level
    */
    void transpose(Ctxt& res, Ctxt& ctxt, zzX*& transpoly);
        
    /**
    @param[out] transpoly The polynomials needed for transpositions of multiple matrices
    */
    void genTransPoly_Parallel(zzX*& transpoly);
        
    /**
    @param[in] ctxt The input ciphertext encrypting multiple matrices
    @param[in] transpoly The polynomials
    @param[out] res The output ciphertext that encrypts the transpose results of the multiple input matrices
    */
    void transpose_Parallel(Ctxt& res, Ctxt& ctxt, zzX*& transpoly);
        
    /**
    @param[in] num The number of steps to shift by columns
    @param[out] shiftpoly The polynomials needed for column shift by num position num = 0: shift by (d-1)
    */
    void genShiftPoly(zzX*& shiftpoly, const long num = 0);
        
    /**
    @param[in] ctxt The input ciphertext, Enc(m[0],...m[d-1] | m[d],...m[2d-1] | ... )
    @param[in] k The number of steps to shift by columns
    @param[in] shiftpoly The polynomials
    @param[out] res The output ciphertext, Enc(m[k] ... m[k-1] | m[d+k]...m[d+k-1] | ... ) consumes a constant-multiplication level
         */
    void shiftBycols(Ctxt& res, Ctxt& ctxt, long k, zzX*& shiftpoly);
       
    /**
    @param[out] shiftpoly The polynomials needed for column shift of multiple matrices
    */
    void genShiftPoly_Parallel(zzX*& shiftpoly);
        
    /**
         @param[in] ctxt The input ciphertext
         @param[in] k The number of steps to shift by columns
         @param[in] shiftpoly The polynomials for parallel shift-by-column operations
         @param[out] res The output ciphertext
    */
    void shiftBycols_Parallel(Ctxt& res, Ctxt& ctxt, long k, zzX*& shiftpoly);
        
        
    //-------------------------------------------
    // multiplication
    //-------------------------------------------
        
    /**
         @param[out] Initpoly The polynomials needed for the initial linear transformations to get A[0] and B[0]
    */
    void genMultPoly(zzX**& Initpoly);
        
    /**
         @param[out] Initpoly The polynomials needed for the initial linear transformations
    */
    void genMultPoly_Parallel(zzX**& Initpoly);
        
    /**
         @param[out] Initpoly The polynomials needed for the initial linear transformations to get B[0]
    */
    void genMultBPoly(zzX*& Initpoly);
    
    /**
     @param[out] Initpoly The polynomials needed for the initial linear transformations to get B[0]
     */
    //void genMultBPoly_Parallel(zzX*& Initpoly);
        
    /**
         @param[in] Actxt The input ciphertext encrypting a matrix A
         @param[in] Bctxt The input ciphertext encrypting a matrix B
         @param[in] Initpoly The polynomials needed for linear transformations of multiplication
         @param[out] resA The output ciphertext encrypting a permuated matrix A[0]
         @param[out] resB The output ciphertext encrypting a permuated matrix B[0]
    */
    void genInitCtxt(Ctxt& resA, Ctxt& resB, Ctxt& Actxt, Ctxt& Bctxt, zzX**& Initpoly);
        
    /**
         @param[in] Actxt The input ciphertext encrypting multiple matrices As
         @param[in] Bctxt The input ciphertext encrypting multiple matrices Bs
         @param[in] Initpoly The polynomials needed for linear transformations of parallel multiplication
         @param[out] resA The output ciphertext encrypting multiple permuated matrices A[0]'s
         @param[out] resB The output ciphertext encrypting multiple permuated matrices B[0]'s
    */
    void genInitCtxt_Parallel(Ctxt& resA, Ctxt& resB, Ctxt& Actxt, Ctxt& Bctxt, zzX**& Initpoly);
        
    /**
         @param[in] Actxts The input ciphertext encrypting a matrix A
         @param[out] mat The output ciphertexts encrypting the permuated matrices A[i] for 0 <= i < dim
    */
    void genInitActxt(vector<Ctxt>& Actxts, Mat<long>& mat);
    
    /**
         @param[in] Actxts The input ciphertext encrypting a matrix A
         @param[out] mat The output ciphertexts encrypting the permuated matrices A[i] for 0 <= i < dim
    */
    //void genInitActxt_Parallel(vector<Ctxt>& Actxts, Mat<ZZ>*& mat);
        
    /**
         @param[in] Bctxt The input ciphertext encrypting a matrix B
         @param[in] Initpoly The polynomials needed for linear transformations of multiplication
         @param[out] resB The output ciphertext encrypting a permuted matrix B[0]
    */
    void genInitBctxt(Ctxt& resB, Ctxt& Bctxt, zzX*& Initpoly);
        
    /**
         @param[in] Actxts The input ciphertext encrypting a wide rectangular matrix A
         @param[out] mat The output ciphertexts encrypting matrices Ai's
    */
    void genInitRecActxt(vector<Ctxt>& Actxts, Mat<long>& mat);
        
    /**
         @param[in] Actxts The input ciphertexts encrypting Ai's
         @param[in] Bctxts The input ciphertexts encrypting Bi's
         @param[in] num The number of ciphertexts to multiply
         @param[out] res The output ciphertext s.t. Actxts[0] * Bctxts[0] + ... + Actxts[d-1] * Bctxts[d-1] followed by a rescaling operation
    */
    void HEmatmul_Hadamard(Ctxt& res, vector<Ctxt>& Actxts, vector<Ctxt>& Bctxts, long num);
        
    /**
         @param[in] Actxt The input ciphertext encrypting a matrix A
         @param[in] Bctxt The input ciphertext encrypting a matrix B
         @param[in] Initpoly The polynomials needed for linear transformations of multiplication
         @param[in] shiftpoly The polynomials for shift-by-column operations
         @param[out] res The output ciphertext encrypting a matrix (A * B)
    */
    void HEmatmul(Ctxt& res, Ctxt& Actxt, Ctxt& Bctxt, zzX**& Initpoly, zzX*& shiftpoly);
        
    /**
         @param[in] Actxt The input ciphertext encrypting multiple matrices As
         @param[in] Bctxt The input ciphertext encrypting multiple matrices Bs
         @param[in] Initpoly The polynomials needed for linear transformations of parallel multiplication
         @param[in] shiftpoly The polynomials for parallel shift-by-column operations
         @param[out] res The output ciphertext encrypting matrices (As * Bs)
    */
    void HEmatmul_Parallel(Ctxt& res, Ctxt& Actxt, Ctxt& Bctxt, zzX**& Initpoly, zzX*& shiftpoly);
        
    /**
         @param[in] Actxt The input ciphertext encrypting a wide rectangular matrix A (an m*n matrix with m<n)
         @param[in] Bctxt The input ciphertext encrypting a square matrix B (an n*n matrix)
         @param[in] Initpoly The polynomials needed for linear transformations of multiplication
         @param[in] shiftpoly The polynomials for shift-by-column operations
         @param[out] res The output ciphertext encrypting a matrix (A * B)
    */
    void HErmatmul(Ctxt& res, Ctxt& Actxt, Ctxt& Bctxt, zzX**& Initpoly, zzX*& shiftpoly);
    
    void HErmatmul_Parallel(Ctxt& res, Ctxt& Actxt, Ctxt& Bctxt, zzX**& Initpoly, zzX*& shiftpoly);
        
    /**
         @param[in] Actxts The input ciphertext encrypting permuted matrices A[i] ("A" are given as fresh ciphertexts)
         @param[in] Bctxt The input ciphertext encrypting a matrix B
         @param[in] Initpoly The polynomials needed for linear transformations of multiplication
         @param[out] res The output ciphertext encrypting a matrix (A * B)
    */
    void HEmatmul_preprocessing(Ctxt& res, vector<Ctxt>& Actxts, Ctxt& Bctxt, zzX*& Initpoly);
        
    /**
     @param[in] Actxts The input ciphertext encrypting permuted matrices A[i] (a rectangular matrix "A" are given as fresh ciphertexts)
     @param[in] Bctxt The input ciphertext encrypting a matrix B
     @param[in] Initpoly The polynomials needed for linear transformations of multiplication
     @param[out] res The output ciphertext encrypting a matrix (A * B)
    */
    void HErmatmul_preprocessing(Ctxt& res, vector<Ctxt>& Actxts, Ctxt& Bctxt, zzX*& Initpoly);
};

#endif /* hematrix_h */
