//
//  main.cpp
//  foSDSC
//
//  Created by Qifan Wang on 12/06/21.
//

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <sys/time.h>
#include <chrono>

#include <NTL/BasicThreadPool.h>
#include <NTL/matrix.h>
#include <openenclave/host.h>
#include <helib/helib.h>

#include "hematrix.h"
#include "LoadMatrix.h"
#include "fosdsc_u.h"

oe_enclave_t* enclave = NULL;

// bool check_simulate_opt(int* argc, const char* argv[])
// {
//     for (int i = 0; i < *argc; i++)
//     {
//         if (strcmp(argv[i], "--simulate") == 0)
//         {
//             cout << "Running in simulation mode" << endl;
//             memmove(&argv[i], &argv[i + 1], (*argc - i) * sizeof(char*));
//             (*argc)--;
//             return true;
//         }
//     }
//     return false;
// }

int CtxtExchange(vector<Ctxt>& resCtxt, vector<Ctxt>& enclaveCtxt, Meta& meta, uint64_t& totalLength, size_t batchIdx){
    oe_result_t result;
    int ret = 0;
    stringstream css;

    // ctxt from the host
    string HostCtxtTemp = "";
    size_t num_ectxt = resCtxt.size();
    size_t ectxt_len = 0;
    uint8_t *ectxt = NULL;
    // ctxt from the enclave
    string EnclaveCtxtTemp = "";
    uint8_t *octxt = NULL;
    size_t octxt_len = 0;
    size_t num_octxt = 0;

    for (size_t i = 0; i < num_ectxt; i++)
    {
        resCtxt[i].writeTo(css);

        // Mat<long> mat;
        // vector<long> cmsg;
        // meta.data->ea.decrypt(resCtxt[i], meta.data->secretKey, cmsg);
        // mat.SetDims(16, 16);
    
        // long k = 0;
        // for(long i = 0; i < 16; ++i){
        //     for(long j = 0; j < 16; ++j){
        //         mat[i][j] = cmsg[k];
        //         k++;
        //     }
        // }
        // cout << "Decrypt resCtxt[" << i << "]: " << endl;
        // cout << mat << endl;
    }
    HostCtxtTemp = css.str();
    ectxt_len = HostCtxtTemp.length();
    cout << "Host: check ectxt_len: " << ectxt_len <<endl;
    ectxt = (uint8_t*) HostCtxtTemp.c_str();
    cout << "Host: check ectxt: " << ectxt <<endl;
    totalLength += ectxt_len;
    cout << "Host: check totalLength after adding ectxt_len: " << totalLength <<endl;

    cout << "Host: send/receive Ctxt to/from enclave:" << endl;
    result = multipleCtxtsTransform(enclave, &ret, ectxt, ectxt_len, num_ectxt, batchIdx, &octxt, &octxt_len, &num_octxt);
    if (result != OE_OK){
        cerr << "Host: transform HE ciphertext failed. result = " << result << endl;
        ret = 1;
        goto exit;
    }
    if (ret != 0)
    {
        cerr << "Host: transform HE ciphertext failed. ret = " << ret << endl;
        goto exit;
    }

    enclaveCtxt = vector<Ctxt>(num_octxt, Ctxt(meta.data->publicKey));
    css.str(std::string());
    css.clear();
    cout << "Host: check octxt_len: " << octxt_len <<endl;
    EnclaveCtxtTemp = std::string(octxt, octxt + octxt_len);
    css << EnclaveCtxtTemp;
    for (size_t i = 0; i < num_octxt; ++i)
    {
        cout << "Host: read " << i << "-th Ctxt from enclave!" << endl;
        enclaveCtxt[i].Ctxt::read(css);
    }
    totalLength += octxt_len;
    cout << "Host: check totalLength after adding octxt_len: " << totalLength <<endl;
    
exit:
    cout << "Host: free strings in CtxtExchange()!" << endl;
    HostCtxtTemp.shrink_to_fit();
    EnclaveCtxtTemp.shrink_to_fit();
    cout << "Host: free memory for ectxt and octxt!" << endl;
    free(octxt);
    return ret;
}

// 2 change places: nrows, Params
int main(int argc, const char * argv[]) {
    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    SetNumThreads(16);
    // Setup the encalve params
    oe_result_t result;
    int ret = 0;
    const uint32_t flags = OE_ENCLAVE_FLAG_DEBUG_AUTO;

    // Setup the Matrix params
    long nrows = 16;
    long ncols = nrows;
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols);

    // Setup the HE params
    cout << "Host: create the Context: " << endl;
    Params param(/*m=*/9472, /*p=*/127, /*r=*/1, /*bits=*/119, /*c=*/2);
    // Params param(/*m=*/16384, /*p=*/6143, /*r=*/1, /*bits=*/180, /*c=*/2);
    // Params param(/*m=*/16384, /*p=*/8191, /*r=*/1, /*bits=*/120, /*c=*/2);
    Meta meta;
    meta(param);

    cout << "Host: print Context contents: " << endl;
    meta.data->context.printout();
    HEmatrix HEmatrix(HEmatpar, meta);

    //setup stringstream
    stringstream ss;

    // Setup the passed context, pk and sk
    string ContextStringTemp = "";
    size_t context_len = 0;
    uint8_t* hecontext = NULL;

    // Receive the raw tree matrix
    string EnclaveTreeTemp = "";
    uint8_t *rawCtxt = NULL;
    size_t rawCtxt_len = 0;
    size_t num_rawCtxt = 0;

    //Send the test data Mat to enclave for classification
    string testTemp = "";
    size_t num_EvaCtxt = 0;
    size_t EvaCtxt_len = 0;
    uint8_t *EvaCtxt = NULL;

    Mat<long>* RawDmatrix; // store the raw data in Mat 
    long numMat = 0; // number of Mat in RawDmatrix
    size_t num_HostCtxt = 1; // number of Ctxt in HostCtxt vector
    size_t queueCtxt_Idx = 0;
    size_t resCtxt_Idx = 0;
    // store the ctxt
    queue<vector<Ctxt>> QCtxt;
    // queue<Ctxt> QCtxt;

    // polynomials used for Mat Multiplication, zzX* for preprocessing
    zzX* Initpoly;
    // zzX** Initpoly;
    // zzX* shiftpoly;

    // Setup the ctxt for following use
    Ctxt queueCtxt(meta.data->publicKey);
    vector<vector<Ctxt>> HostCtxt;
    vector<Ctxt> EnclaveCtxt;
    vector<Ctxt> resCtxt;
    // Actxts: used in Mul_preprocessing
    vector<Ctxt> Actxts;
    // tempActxts used in Mul_preprocessing
    vector<Ctxt> tempActxts;

    // record the communication cost and computational time
    uint64_t totalLength = 0;
    auto start= chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = 0.0;
    double ClientTimeTemp = 0.0;

    /*---------------------------------------*/
    //  Create the enclave
    /*---------------------------------------*/
    // Command Format: ./host/fosdsc_host ./enclave/enclave.signed

    // Check the simulation
    // if (check_simulate_opt(&argc, argv))
    // {
    //     flags |= OE_ENCLAVE_FLAG_SIMULATE;
    // }

    cout << "Host: create enclave for image:" << argv[1] << endl;
    result = oe_create_fosdsc_enclave(argv[1], OE_ENCLAVE_TYPE_SGX, flags, NULL, 0, &enclave);
    if (result != OE_OK)
    {
        cerr << "oe_create_fosdsc_enclave() failed with " << argv[0] << " " << result << endl;
        ret = 1;
        goto exit;
    }

    /*-------------------------------------------------------------------*/
    //  Call into the enclave to transform the HE params
    /*-------------------------------------------------------------------*/

    // send context to enclave
    meta.data->context.writeTo(ss);
    meta.data->secretKey.writeTo(ss);
    meta.data->publicKey.writeTo(ss);
    ContextStringTemp = ss.str();
    context_len = ContextStringTemp.size();
    hecontext = (uint8_t*)ContextStringTemp.c_str();
    // cout << "Host: check hecontext: " << hecontext << endl;
    cout << "Host: check context_len, containing context, sk and pk: " << context_len << endl;

    cout << "Host: transform HE params (context, sk and pk) into enclave:" << endl;
    result = enclave_init(enclave, &ret, hecontext, context_len, &rawCtxt, &rawCtxt_len, &num_rawCtxt);
    if (result != OE_OK){
        cerr << "Host: transform HE params failed. OE result = " << result << endl;
        ret = 1;
        goto exit;
    }
    if (ret != 0)
    {
        cerr << "Host: transform HE params failed. ret = " << ret << endl;
        goto exit;
    }

    cout << "Host: receive Root Tree Mats from enclave:" << endl;
    // receive root tree Mat after enclave initialization
    EnclaveCtxt = vector<Ctxt>(num_rawCtxt, Ctxt(meta.data->publicKey));
    ss.str(std::string());
    ss.clear();
    EnclaveTreeTemp = std::string(rawCtxt, rawCtxt + rawCtxt_len);
    ss << EnclaveTreeTemp;
    for (size_t i = 0; i < num_rawCtxt; ++i)
    {
        cout << "Host: read " << i << "-th Root Tree Ctxt from enclave!" << endl;
        EnclaveCtxt[i].Ctxt::read(ss);

        // Mat<long> test;
        // HEmatrix.decryptZmat(test, EnclaveCtxt[i]);
        // cout << "Decrypt root node EnclaveCtxt[" << i << "]: " << endl;
        // cout << test << endl;
    }

    totalLength += rawCtxt_len;
    cout << "Host: the length of Root Ctxts from Enclave: " << rawCtxt_len << endl;
    cout << "Host: Check the length of totalLength: " << totalLength << endl;
    
    /*---------------------------------------*/
    //  Load DataMatrix
    /*---------------------------------------*/
    cout << "Host: Load the raw data matrix." << endl;
    ProcessDMatrix(RawDmatrix, numMat, nrows, true);
    // cout << "Host: Check the first sub Mat in RawDmatrix: " << endl;
    // cout << RawDmatrix[0] << endl;

    /*---------------------------------------*/
    //  GenPoly
    /*---------------------------------------*/
    cout << "Host: Generate the Multi Poly." << endl;
    HEmatrix.genMultBPoly(Initpoly);
    // HEmatrix.genMultPoly(Initpoly);
    // cout << "Host: Generate the Shift Poly." << endl;
    // HEmatrix.genShiftPoly(shiftpoly);

    /*---------------------------------------*/
    //  Encryption
    /*---------------------------------------*/
    // This part should be the cost of client
    cout << "Host: encrypt the matrix and send them to the queue." << endl;
    start= chrono::steady_clock::now();
    for (size_t i = 0; i < numMat; ++i)
    {
        HEmatrix.genInitActxt(Actxts, RawDmatrix[i]);
        // cout << "Host: check Actxts size (should always be 16): " << Actxts.size() << endl;
        QCtxt.push(Actxts);
    }
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Host: Client Encryption Time= " << timeElapsed << " s" << endl;
    cout << "------------------------------------------------------------------------" << endl;
    // cout << "Host: check QCtxt size: " << QCtxt.size() << endl;
    ClientTimeTemp = timeElapsed;

    /*---------------------------------------*/
    //  Multiplication
    /*---------------------------------------*/
    // Confirm hte batch size
    num_HostCtxt = 64 / nrows;
    cout << "Host: check num_HostCtxt: " << num_HostCtxt << endl;

    start= chrono::steady_clock::now();
    while (!QCtxt.empty())
    {
        if (QCtxt.size() < num_HostCtxt)
        {
            queueCtxt_Idx = QCtxt.size(); // may < num_HostCtxt
        }else{
            queueCtxt_Idx = num_HostCtxt;
        }
        cout << "Host: check final queueCtxt_Idx: " << queueCtxt_Idx << endl;
        
        // ensure HostCtxt store the latest ctxts
        HostCtxt.clear();
        // HostCtxt = vector<Ctxt>(num_HostCtxt, Ctxt(meta.data->publicKey));
        for(size_t i = 0; i < queueCtxt_Idx; ++i)
        {
            HostCtxt.push_back(QCtxt.front());
            QCtxt.pop();
        }
        
        // cout << "Host: check HostCtxt size: " << HostCtxt.size() << endl;
        // cout << "Host: check EnclaveCtxt size: " << EnclaveCtxt.size() << endl;
        // cout << "Host: check resCtxt size: " << (HostCtxt.size() * EnclaveCtxt.size()) << endl;
        resCtxt = vector<Ctxt>((HostCtxt.size() * EnclaveCtxt.size()), Ctxt(meta.data->publicKey));

        // NTL_EXEC_RANGE(HostCtxt.size(), first, last);
        // for (size_t i = first; i < last; ++i)
        for (size_t i = 0; i < HostCtxt.size(); ++i)
        {
            // here if we use multi-thread, cause low performance, or try more threads
            // NTL_EXEC_RANGE(EnclaveCtxt.size(), first, last);
            // for (size_t j = first; j < last; ++j)
            for (size_t j = 0; j < EnclaveCtxt.size(); ++j)
            {
                // Ensure that HostCtxt[i] is not changed in last HEmatmul_preprocessing;
                tempActxts.clear();
                tempActxts.assign(HostCtxt[i].begin(), HostCtxt[i].end());
                // tempActxts[j].assign(HostCtxt[i].begin(), HostCtxt[i].end());
                cout << "Host: " << i << "-th HostCtxt, " << j << "-th EnclaveCtxt, Mat Multiplication: " << endl;
                // HEmatrix.HEmatmul_preprocessing(resCtxt[(i * EnclaveCtxt.size() + j)], tempActxts[j], EnclaveCtxt[j], Initpoly);
                HEmatrix.HEmatmul_preprocessing(resCtxt[(i * EnclaveCtxt.size() + j)], tempActxts, EnclaveCtxt[j], Initpoly);
                // HEmatrix.HEmatmul(resCtxt[(i * EnclaveCtxt.size() + j)], HostCtxt[i], EnclaveCtxt[j], Initpoly, shiftpoly);
            }
            // NTL_EXEC_RANGE_END;
        }
        // NTL_EXEC_RANGE_END;

        // clear the EnclaveCtxt, waiting for the ctxt from the enclave
        //EnclaveCtxt.clear();
        ret = CtxtExchange(resCtxt, EnclaveCtxt, meta, totalLength, queueCtxt_Idx);
        if (ret != 0)
        {
            cerr << "Host: CtxtExchange failed with " << ret << endl;
            goto exit;
        }
    }
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Host: Total Training Time = " << timeElapsed << " s" << endl;
    cout << "Host: Total Training Time (min) = " << (timeElapsed / 60) << " min" << endl;
    cout << "------------------------------------------------------------------------" << endl;

    cout << "Host: Total Communication Cost = " << ((double) totalLength / (double)(1024 * 1024)) << " MB" << endl;
    cout << "Host: Client Encryption Time= " << ClientTimeTemp << " s" << endl;

    start= chrono::steady_clock::now();
    num_EvaCtxt = 1;
    testTemp = "0101011001010100";
    EvaCtxt_len = testTemp.length();
    // cout << "Host: check EvaCtxt_len: " << EvaCtxt_len <<endl;
    EvaCtxt = (uint8_t*) testTemp.c_str();
    // cout << "Host: check EvaCtxt: " << EvaCtxt <<endl;
    // invoke ECALL for classification
    result = HT_Classify(enclave, &ret, EvaCtxt, EvaCtxt_len, num_EvaCtxt);
    if (result != OE_OK){
        cerr << "Host: calling into enclave for classification failed. OE result = " << result << endl;
        ret = 1;
        goto exit;
    }
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Host: Total Evaluation Time= " << timeElapsed << " s" << endl;
    cout << "------------------------------------------------------------------------" << endl;

exit:
    // cannot free, otherwise cause double free corruption
    // cout << "Host: free memory for context!" << endl;
    // free(hecontext);
    // cout << "Host: free memory for rawCtxt!" << endl;
    // free(rawCtxt);
    cout << "Host: called close_encryptor" << endl;
    result = close_encryptor(enclave);
    if (result != OE_OK)
    {
        ret = 1;
    }
    cout << "Host: terminate the enclave" << endl;
    oe_terminate_enclave(enclave);
    cout << "Host: done  " << ((ret == 0) ? "succeeded" : "failed") << endl;
    return ret;
}
