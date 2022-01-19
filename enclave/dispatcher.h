#pragma once

#include <cstring>
#include <sstream>
#include <math.h>
#include <unordered_map>
#include <map>
#include <helib/helib.h>
#include <openenclave/enclave.h>

#include "mHT/HoeffdingTree.h"

using namespace std;
using namespace helib;
using namespace NTL;

typedef struct _HTparams
{
  size_t numFeatures;
  size_t numClasses;
  size_t numBits; // synathetic dataset, should be 16, 32, 64
  vector<size_t> numAllCategories;
  vector<unordered_map<size_t, string>> FeaCatMapping;
  double successProbability;
  size_t maxSamples;
  size_t checkInterval;
  size_t minSamples;
}HTparams;

class ecall_dispatcher
{
  private:
    // HE context pointer
    Context* e_context;
    // HE secret key
    unique_ptr<SecKey> activeSecKey;
    // HE public key
    unique_ptr<PubKey> activePubKey;
    // Receive Ctxts from the host 
    // vector<Ctxt> ReceivedCtxts;
    // Return Ctxts to the enclave
    // vector<Ctxt> ReturnedCtxts;
    // Dataset info used in train and evaluation
    HTparams* e_HTparams;
    // HTparams* e_HTparams;
    // store the current leaf nodes
    list<HoeffdingTree*> leafnodes;
    // Initialize is 2, rawTreeMat in the host
    size_t num_returnedMat = 0;
    // Mat dimension (hard coded)
    size_t dim = 16;
    // Store current tree matrix, set it as global for multiple classifications
    Mat<size_t> TreeMat;


  public:
    // ecall_dispatcher(HTparams* params);
    ecall_dispatcher();
    int enclave_init(uint8_t* hecontext, size_t context_len, uint8_t** rawCtxt, size_t* rawCtxt_len, size_t* num_rawCtxt);
    // int CtxtTransform( uint8_t* ectxt, uint8_t** octxt, size_t ectxt_len, size_t* octxt_len);
    int multipleCtxtsTransform(uint8_t* ectxt, size_t ectxt_len, size_t num_ectxt, size_t batchIdx, uint8_t** octxt, size_t* octxt_len, size_t* num_octxt);
    int HT_Classify(uint8_t* EvaCtxt, size_t EvaCtxt_len, size_t num_EvaCtxt);
    void close();
  private:
    int HTtrain(vector<Ctxt>& ReceivedCtxts, vector<Ctxt>& ReturnedCtxts, size_t batchIdx);
    // Decrypt the ReceivedCtxts vector, generate a complete Mat
    void DecryptMat(vector<Ctxt>& ReceivedCtxts, vector<map<size_t, Mat<size_t>>>& countedMat, vector<pair<size_t, size_t>>& indiceAndsize, size_t batchIdx);
    // TRansform the leaf strings to Mats, encrypt the Mats
    void EncryptMat(vector<Ctxt>& ReturnedCtxts, vector<string>& Allleafmatrix);
    // Helper function to initialize Mat with 0
    void MatInit(Mat<size_t>& matInit);
    // copy uint8_t array
    void copy_arr_to_enclave(uint8_t* dst[], size_t num, uint8_t* src[], size_t lengths[]);
    void free_array(uint8_t* arr[], size_t len);
    // void* oe_host_malloc(size_t size);
};
