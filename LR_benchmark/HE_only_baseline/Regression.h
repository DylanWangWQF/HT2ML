#ifndef Regression_h
#define Regression_h

#include "PtMatrix.h"
#include "PtMatrix.cpp"

#include <helib/helib.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/matrix.h>

#include <sys/time.h>
#include <chrono>
#include <vector>
#include <fstream>

using namespace std;
using namespace NTL;
using namespace helib;

struct LRparams
{
    const long m, p, r, bits, c;
    const vector<long> mvec;
    const vector<long> gens;
    const vector<long> ords;
    LRparams(long _m,
            long _p,
            long _r,
            long _bits,
            long _c,
            const vector<long>& _gens = {},
            const vector<long>& _ords = {},
            const vector<long>& _mvec = {}) :
      m(_m), p(_p), r(_r), bits(_bits), c(_c), gens(_gens), ords(_ords), mvec(_mvec) {}
    LRparams(const LRparams& other) :
    LRparams(other.m,
            other.p,
            other.r,
            other.bits,
            other.c,
            other.gens,
            other.ords,
            other.mvec) {}
    bool operator!=(LRparams& other) const { return !(*this == other); }
    bool operator==(LRparams& other) const
    {
        return m == other.m && p == other.p && r == other.r && bits == other.bits && c == other.c && gens == other.gens && ords == other.ords && mvec == other.mvec;
    }
};

struct LR_ContextAndKeys
{
  const LRparams params;
  Context context;
  SecKey secretKey;
  const PubKey& publicKey;
  const EncryptedArray& ea;

    LR_ContextAndKeys(LRparams& _params) :
        params(_params),
        context(ContextBuilder<BGV>()
                    .m(params.m)
                    .p(params.p)
                    .r(params.r)
                    .bits(params.bits)
                    .c(params.c)
                    .gens(params.gens)
                    .ords(params.ords)
                    .bootstrappable(!params.mvec.empty())
                    .mvec(params.mvec)
                    .thickboot()
                    .build()),
        secretKey(context),
        publicKey((secretKey.GenSecKey(),
                   addSome1DMatrices(secretKey),
                   addFrbMatrices(secretKey),
                   secretKey.genRecryptData(),
                   secretKey)),
        ea(context.getEA()){}
};

struct LRmeta
{
  std::unique_ptr<LR_ContextAndKeys> data;
    LRmeta& operator()(LRparams& params)
  {
    // Only change if nullptr or different.
    if (data == nullptr || data->params != params)
      data = std::make_unique<LR_ContextAndKeys>(params);
    return *this;
  }
};

void RegressPT(vector<long> &theta, long &det, Mat<long> &data, vector<long> &labels);

bool LoadData(Mat<long> &rawData, vector<long> &labels, long &dim, string &filename);

void BatchData(vector<vector<vector<long>>> &ptxtData, vector<vector<long>> &ptxtLabels, const Mat<long> &rawData, const vector<long> &labels, long p, long nslots);

class Regression{
public:
    LRmeta& meta;
    Regression(LRmeta& meta) : meta(meta) {}
    
    void EncryptData(const vector<vector<vector<long>>> &ptxtData, const vector<vector<long>> &ptxtLabels);
    void Regress(vector<Ctxt> &theta, Ctxt &det);
    void Clear();
    void CheckWithThinBoot(Ctxt& c);
    
    // utils (i.e., mul, inverse)
    void MultiplyMatrix(vector<vector<Ctxt>> &output, vector<vector<Ctxt>> &input1, vector<vector<Ctxt>> &input2, bool isTranspose1, bool isTranspose2);
    void MultiplyVector(vector<vector<Ctxt>> &output, vector<vector<Ctxt>> &input1, vector<Ctxt> &input2, bool isTranspose1, bool isTranspose2);
    void MultByTranspose(vector<vector<Ctxt>> &output, vector<vector<Ctxt>> &input, bool isTranspose);
    void Invert(Ctxt &det, vector<vector<Ctxt>> &adj, vector<vector<Ctxt>> &input, bool isTranspose);
    void Determinant(Ctxt& det, vector<bool> &usedRows, vector<bool> &usedCols, long dim, vector<vector<Ctxt>> &input, bool isTranspose);
    
    vector<Ctxt> labels;
    vector<vector<Ctxt>> data;
};

#endif /* Regression_h */
