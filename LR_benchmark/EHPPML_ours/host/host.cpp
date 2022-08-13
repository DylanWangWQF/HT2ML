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
#include "LoadData.h"
#include "fosdsc_u.h"

oe_enclave_t* enclave = NULL;

int CallEnclaveForMatrix(Ctxt& result1, Ctxt& result2, Ctxt& enclaveCtxt, Meta& meta, uint64_t& totalLength)
{
    oe_result_t result;
    int ret = 0;
    stringstream css;

    // ctxt from the host
    string HostCtxtTemp = "";
    size_t ectxt_len = 0;
    uint8_t *ectxt = NULL;
    // ctxt from the enclave
    string EnclaveCtxtTemp = "";
    uint8_t *octxt = NULL;
    size_t octxt_len = 0;

    result1.writeTo(css);
    result2.writeTo(css);
    HostCtxtTemp = css.str();
    ectxt_len = HostCtxtTemp.length();
    ectxt = (uint8_t*) HostCtxtTemp.c_str();
    totalLength += ectxt_len;

    // cout << "Host: send/receive Ctxt to/from enclave:" << endl;
    result = MatrixOperation(enclave, &ret, ectxt, ectxt_len, &octxt, &octxt_len);
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

    css.str(std::string());
    css.clear();
    EnclaveCtxtTemp = std::string(octxt, octxt + octxt_len);
    css << EnclaveCtxtTemp;
    enclaveCtxt.Ctxt::read(css);
    totalLength += octxt_len;
    
exit:
    // cout << "Host: free memory for ectxt and octxt!" << endl;
    free(octxt);
    return ret;
}

int main(int argc, const char * argv[]) {
    SetNumThreads(16);
    /*---------------------------------------*/
    //  Load the dataset
    /*---------------------------------------*/
    cout << "Host: Loading dataset..." << endl;
    long nrows = 16;
    long ncols = nrows;
    Mat<long>* DataMat;
    Mat<long>* TranDataMat;
    Mat<long>* LabelMat;
    long num_DataMat;
    // string dataset = "../../scripts/2_dim_LR.dat";
    // string dataset = "../../scripts/4_dim_LR.dat";
    string dataset = "../../scripts/6_dim_LR.dat";
    // string dataset = "../../scripts/8_dim_LR.dat";
    // string dataset = "../../scripts/16_dim_LR.dat";
    ProcessDataMatrix(DataMat, TranDataMat, LabelMat, num_DataMat, nrows, dataset);

    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    // Setup the encalve params
    oe_result_t result;
    int ret = 0;
    const uint32_t flags = OE_ENCLAVE_FLAG_DEBUG_AUTO;

    // Setup the Matrix params
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols);
    Params param(/*m=*/9472, /*p=*/127, /*r=*/1, /*bits=*/119, /*c=*/2);
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

    // Setup ctxts
    zzX* Initpoly;
    vector<Ctxt> Actxts; // Actxts: used in Mul_preprocessing
    vector<vector<Ctxt>> EncTranDataMat;
    vector<vector<Ctxt>> EncTranDataMat_2nd;
    vector<Ctxt> EncDataMat(num_DataMat, Ctxt(meta.data->publicKey));
    vector<Ctxt> EncLabelMat(num_DataMat, Ctxt(meta.data->publicKey));
    vector<Ctxt> EncResultMat1(num_DataMat, Ctxt(meta.data->publicKey));
    vector<Ctxt> EncResultMat2(num_DataMat, Ctxt(meta.data->publicKey));
    Ctxt result1(meta.data->publicKey);
    Ctxt result2(meta.data->publicKey);
    Ctxt enclaveCtxt(meta.data->publicKey);

    // record the communication cost and run time
    uint64_t totalLength = 0;
    string ClientTemp = "";
    uint64_t client_totalLength = 0;

    auto start= std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = 0.0;
    
    auto host_end = std::chrono::steady_clock::now();
    auto host_diff = host_end - start;
    double host_timeElapsed = 0.0;

    auto client_start= std::chrono::steady_clock::now();
    auto client_end = std::chrono::steady_clock::now();
    auto client_diff = end - start;
    double client_timeElapsed = 0.0;

    auto server_start= std::chrono::steady_clock::now();
    auto server_end = std::chrono::steady_clock::now();
    auto server_diff = end - start;
    double server_timeElapsed = 0.0;

    /*---------------------------------------*/
    //  Create the enclave
    /*---------------------------------------*/
    // Command Format: ./host/fosdsc_host ./enclave/enclave.signed
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
    cout << "Host: check size of context and keys: " << ((double) context_len / (double)(1024 * 1024)) << " MB" << endl;

    cout << "Host: transform HE params (context, sk and pk) into enclave:" << endl;
    result = enclave_init(enclave, &ret, hecontext, context_len);
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

    /*---------------------------------------*/
    //  Pre-process polynominals
    /*---------------------------------------*/
    HEmatrix.genMultBPoly(Initpoly);

    /*---------------------------------------*/
    //  Encrypt the matrix
    /*---------------------------------------*/
    cout << "Encrypting the data and labels..." << endl;
    client_start= chrono::steady_clock::now();

    for (long i = 0; i < num_DataMat; ++i)
    {
        HEmatrix.encryptZmat(EncDataMat[i], DataMat[i]);
        HEmatrix.encryptZmat(EncLabelMat[i], LabelMat[i]);
        HEmatrix.genInitActxt(Actxts, TranDataMat[i]);
        EncTranDataMat.push_back(Actxts);
        EncTranDataMat_2nd.push_back(Actxts);
    }

    client_end = std::chrono::steady_clock::now();
    client_diff = client_end - client_start;
    client_timeElapsed = chrono::duration <double, milli> (client_diff).count()/1000.0;
    // record the size of generated ctxts
    ss.str(std::string());
    ss.clear();
    for (long i = 0; i < num_DataMat; i++)
    {
        EncDataMat[i].writeTo(ss);
        EncLabelMat[i].writeTo(ss);
        for (long  j = 0; j < EncTranDataMat[0].size(); j++)
        {
            EncTranDataMat[i][j].writeTo(ss);
        }
    }
    ClientTemp = ss.str();
    client_totalLength = ClientTemp.size();
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Host: Client Encryption Time = " << client_timeElapsed << " s" << endl;
    cout << "Host: Client Communication Cost = : " << ((double) client_totalLength / (double)(1024 * 1024)) << " MB" << endl;
    cout << "------------------------------------------------------------------------" << endl;

    /*---------------------------------------*/
    //  Homomorphically perform Linear Regression Training
    /*---------------------------------------*/
    // cout << endl << "TranDataMat * DataMat and Send/Receive ctxts to/from the enclave..." << endl;
    start= chrono::steady_clock::now();
    // 1. TranDataMat * DataMat
    for (long k = 0; k < num_DataMat; k++)
    {
        HEmatrix.HEmatmul_preprocessing(EncResultMat1[k], EncTranDataMat[k], EncDataMat[k], Initpoly);
        result1 += EncResultMat1[k];
        HEmatrix.HEmatmul_preprocessing(EncResultMat2[k], EncTranDataMat_2nd[k], EncLabelMat[k], Initpoly);
        result2 += EncResultMat2[k];
    }
    host_end = std::chrono::steady_clock::now();

    // 2. Send/Receive ctxts to/from the enclave
    // cout << endl << "Send/Receive ctxts to/from the enclave..." << endl;
    CallEnclaveForMatrix(result1, result2, enclaveCtxt, meta, totalLength);

    // Calculate host runtime
    host_diff = host_end - start;
    host_timeElapsed = chrono::duration <double, milli> (host_diff).count()/1000.0;

    // Calculate total runtime
    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;

    cout << "------------------------------------------------------------------------" << endl;
    cout << "Host: RunTime in Host = " << host_timeElapsed << " s" << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Host: LR Training Time= " << timeElapsed << " s" << endl;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Host: Communication Cost between Server and Enclave = : " << ((double) totalLength / (double)(1024 * 1024)) << " MB" << endl;
    cout << "------------------------------------------------------------------------" << endl;

exit:
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
