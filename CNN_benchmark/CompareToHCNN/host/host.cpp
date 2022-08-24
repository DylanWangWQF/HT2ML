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

#include "CKKSmatrix.h"
#include "LoadData.h"
#include "fosdsc_u.h"

oe_enclave_t* enclave = NULL;

int CallEnclaveForNonLin(vector<vector<Ctxt>>& resCtxt, vector<Ctxt>& enclaveCtxt, ckksMeta& meta, uint64_t& totalLength)
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

    for (size_t i = 0; i < 6; i++)
    {
        for (size_t j = 0; j < 9; j++)
        {
            resCtxt[i][j].writeTo(css);
        }
    }
    HostCtxtTemp = css.str();
    ectxt_len = HostCtxtTemp.length();
    ectxt = (uint8_t*) HostCtxtTemp.c_str();
    totalLength += ectxt_len;

    // cout << "Host: send/receive Ctxt to/from enclave:" << endl;
    result = nonLinearLayer(enclave, &ret, ectxt, ectxt_len, &octxt, &octxt_len);
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
    for (size_t i = 0; i < 14; ++i)
    {
        enclaveCtxt[i].Ctxt::read(css);
    }
    totalLength += octxt_len;
    
exit:
    // cout << "Host: free memory for ectxt and octxt!" << endl;
    free(octxt);
    return ret;
}

int main(int argc, const char * argv[]) {
    /*---------------------------------------*/
    //  Load the dataset
    //  Crop the images
    /*---------------------------------------*/
    cout << "Loading test imgs & labels..." << endl;
    int MNIST_HEIGHT = 28;
    int MNIST_WIDTH = 28;
    vector<vector<float>> test_imgs;
    vector<unsigned char> test_lbls;
    int num_imgs = 64;
    int num_channels = 6;
    int kernel_size = 5;
    int stride = 1;
    
    // Extract windows from 28 * 28 images (e.g., imgMat[0] ---> first image)
    int window_size = ceil((MNIST_HEIGHT - kernel_size) / stride) + 1; // 24
    int num_windows = window_size * window_size; // 24 * 24 = 576
    int num_kernels = kernel_size * kernel_size; // 25
    Mat<RR>* feaMat = new Mat<RR>[num_kernels]; // matrix for ct.I_{i,j}

    CropImages(test_imgs, test_lbls, feaMat, MNIST_HEIGHT, MNIST_WIDTH, num_imgs, kernel_size, stride, window_size);
    
    /*---------------------------------------*/
    //  Load the model
    /*---------------------------------------*/
    cout << "Loading the model params..." << endl;
    int kernel_dim1, kernel_dim2;
    int dense1_dim1, dense1_dim2;

    Mat<RR>* kernel_weights;
    Mat<RR>* dense1_weights;

    LoadModel(kernel_weights, dense1_weights, kernel_dim1, kernel_dim2, dense1_dim1, dense1_dim2, num_imgs);
    
    /*-----------------------------------------------Load images and models end!----------------------------------------------------------------------*/

    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    SetNumThreads(16);
    // Setup the encalve params
    oe_result_t result;
    int ret = 0;
    const uint32_t flags = OE_ENCLAVE_FLAG_DEBUG_AUTO;

    // Setup the Matrix params
    ckksMatpar ckksmatpar;
    long ncols = 64, nrows = 64, subdim = 16;
    readckksMatpar(ckksmatpar, nrows, ncols, subdim); // subdim used in below RecMul
    ckksParams param(/*m=*/16 * 1024, /*bits=*/179, /*precision=*/20, /*c=*/2);
    // ckksParams param(/*m=*/16 * 1024, /*bits=*/235, /*precision=*/20, /*c=*/4);
    ckksMeta meta;
    meta(param);
    cout << "HE Context contents: " << endl;
    meta.data->context.printout();
    CKKSmatrix CKKSmatrix(ckksmatpar, meta);

    //setup stringstream
    stringstream ss;

    // Setup the passed context, pk and sk
    string ContextStringTemp = "";
    size_t context_len = 0;
    uint8_t* hecontext = NULL;

    // Setup ctxts during inference
    vector<vector<Ctxt>> ct_I(num_kernels, vector<Ctxt>(9, Ctxt(meta.data->publicKey))); // 576 / 64 = 9, vector size = 25 * 9
    vector<Ctxt> ct_K(kernel_dim1 * kernel_dim2, Ctxt(meta.data->publicKey)); // 30 * 5
    // original-HE
    // vector<Ctxt> ct_W(14, Ctxt(meta.data->publicKey));
    // optimised-HE
    vector<vector<Ctxt>> ct_W;
    vector<vector<Ctxt>> ct_Ck(num_channels, vector<Ctxt>(9, Ctxt(meta.data->publicKey)));
    vector<Ctxt> ct_enclave(14, Ctxt(meta.data->publicKey));
    // original-HE
    // vector<vector<EncodedPtxt>> Initpoly;
    // vector<EncodedPtxt> shiftpoly;
    // optimised-HE
    vector<EncodedPtxt> Initpoly;
    vector<Ctxt> rmul_temp(14, Ctxt(meta.data->publicKey));
    Ctxt inference_result(meta.data->publicKey);

    // record the communication cost and run time
    uint64_t totalLength = 0;
    string ClientTemp = "";
    uint64_t client_totalLength = 0;
    string ServerTemp = "";
    uint64_t server_totalLength = 0;

    auto start= std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = 0.0;

    auto client_start= std::chrono::steady_clock::now();
    auto client_end = std::chrono::steady_clock::now();
    auto client_diff = end - start;
    double client_timeElapsed = 0.0;

    auto server_start= std::chrono::steady_clock::now();
    auto server_end = std::chrono::steady_clock::now();
    auto server_diff = end - start;
    double server_timeElapsed = 0.0;

    // Mat<RR> test;

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
    // original-HE
    // CKKSmatrix.genMultPoly(Initpoly);
    // CKKSmatrix.genShiftPoly(shiftpoly);
    // optimised-HE
    CKKSmatrix.genMultBPoly(Initpoly);

    /*---------------------------------------*/
    //  Encrypt the images
    /*---------------------------------------*/
    cout << "Encrypting the images (feaMat, size = 7 * 7)..." << endl;
    client_start= chrono::steady_clock::now();
    for (int i = 0; i < num_kernels; i++) // 25
    {
        for (int  j = 0; j < 9; j++)
        {
            Mat<RR> mat_temp;
            mat_temp.SetDims(64, 64);
            for (int m = 0; m < 64; m++)
            {
                for (int n = 0; n < 64; n++)
                {
                    mat_temp[m][n] = feaMat[i][(m + j * 64)][n];
                }
            }
            CKKSmatrix.encryptRmat(ct_I[i][j], mat_temp);
        }
    }
    client_end = std::chrono::steady_clock::now();
    client_diff = client_end - client_start;
    client_timeElapsed = chrono::duration <double, milli> (client_diff).count()/1000.0;
    ss.str(std::string());
    ss.clear();
    for (int i = 0; i < num_kernels; i++) // 25
    {
        for (int  j = 0; j < 9; j++)
        {
            ct_I[i][j].writeTo(ss);
        }
    }
    ClientTemp = ss.str();
    client_totalLength = ClientTemp.size();
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Host: Client Encryption Time = " << client_timeElapsed << " s" << endl;
    cout << "Host: Client Communication Cost = : " << ((double) client_totalLength / (double)(1024 * 1024)) << " MB" << endl;
    cout << "------------------------------------------------------------------------" << endl;
    
    /*---------------------------------------*/
    //  Encrypt the model params
    /*---------------------------------------*/
    // 1. ct.K_{i, j}_{k}
    // cout << "Preparing model params ciphertext, ct_K..." << endl;
    server_start= chrono::steady_clock::now();
    for (int i = 0; i < ct_K.size(); i++) // 30 * 5
    {
        CKKSmatrix.encryptRmat(ct_K[i], kernel_weights[i]);
    }

    // 2. ct.W_{k}, preprocessed ctxt
    // cout << "Preparing model params ciphertext, ct_W..." << endl;

    // original-HE
    // for (int i = 0; i < 14; i++) // 864/64=13.5=14
    // {
    //     CKKSmatrix.encryptRmat(ct_W[i], dense1_weights[i]);
    // }

    // optimised-HE
    for (int i = 0; i < 14; i++) // 864/64=13.5=14
    {
        vector<Ctxt> temp;
        CKKSmatrix.genInitRecActxt(temp, dense1_weights[i]);
        ct_W.push_back(temp);
    }
    
    server_end = std::chrono::steady_clock::now();
    server_diff = server_end - server_start;
    server_timeElapsed = chrono::duration <double, milli> (server_diff).count()/1000.0;

    ss.str(std::string());
    ss.clear();
    for (int i = 0; i < ct_K.size(); i++)
    {
        ct_K[i].writeTo(ss);
    }

    // original-HE
    // for (int i = 0; i < 14; i++)
    // {
    //     ct_W[i].writeTo(ss);  
    // }

    // optimised-HE
    for (int i = 0; i < 14; i++)
    {
        for (int j = 0; j < ct_W[i].size(); j++)
        {
            ct_W[i][j].writeTo(ss);
        }  
    }

    ServerTemp = ss.str();
    server_totalLength = ServerTemp.size();
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Host: Server Encryption Time = " << server_timeElapsed << " s" << endl;
    cout << "Host: Server Communication Cost = : " << ((double) server_totalLength / (double)(1024 * 1024)) << " MB" << endl;
    cout << "------------------------------------------------------------------------" << endl;

    /*---------------------------------------*/
    //  Homomorphically perform CNN inference
    /*---------------------------------------*/
    start= chrono::steady_clock::now();
    // 1. HE conv layer, ct.I_{i, j} and ct.K_{i, j}_{k}
    // cout << endl << "Inference-convolution layer..." << endl;
    for (int k = 0; k < num_channels; k++) // 6
    {
        for (int j = 0; j < 9; j++)
        {
            for (int i = 0; i < num_kernels; i++) // 25
            {
                Ctxt tmp = ct_I[i][j];
                tmp.multLowLvl(ct_K[k * num_kernels + i]);
                ct_Ck[k][j] += tmp;           
            }
            ct_Ck[k][j].reLinearize(); // return to a canonical state
        }
    }

    // 2. activation and pooling layer
    // cout << endl << "Inference-activation and pooling layer..." << endl;
    // should output 864 * 64 => 14 of 64 * 64
    CallEnclaveForNonLin(ct_Ck, ct_enclave, meta, totalLength);

    // CKKSmatrix.decryptRmat(test, ct_enclave[0]);
    // cout << endl << "Generated matrix ct_enclave[0]: " << endl << test << endl; 

    // 3. HE FC-2 layer, 10 * 64 X 64 * 64 = 10 * 64, ct.V and ct_F
    // cout << endl << "Inference-FC layer..." << endl;

    // original-HE
    // for (int i = 0; i < 14; i++)
    // {
    //     CKKSmatrix.HErmatmul(rmul_temp[i], ct_W[i], ct_enclave[i], Initpoly, shiftpoly);
    //     inference_result += rmul_temp[i];
    // }

    // optimised-HE
    for (int i = 0; i < 14; i++)
    {
        CKKSmatrix.HErmatmul_preprocessing(rmul_temp[i], ct_W[i], ct_enclave[i], Initpoly);
        inference_result += rmul_temp[i];
    }
    // inference_result.bumpNoiseBound(1e-7);

    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Host: Inference Time= " << timeElapsed << " s" << endl;
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
