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

int RefreshSingleCtxt(Ctxt& resCtxt, uint64_t& totalLength)
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

    resCtxt.writeTo(css);

    HostCtxtTemp = css.str();
    ectxt_len = HostCtxtTemp.length();
    ectxt = (uint8_t*) HostCtxtTemp.c_str();
    totalLength += ectxt_len;

    // cout << "Host: send/receive Ctxt to/from enclave:" << endl;
    result = singleCtxtTransform(enclave, &ret, ectxt, ectxt_len, &octxt, &octxt_len);
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
    resCtxt.Ctxt::read(css);
    totalLength += octxt_len;

exit:
    cout << "Host: free memory for ectxt and octxt!" << endl;
    free(octxt);
    return ret;
}

int RefreshMultipleCtxts(vector<Ctxt>& resCtxt, uint64_t& totalLength, ckksMeta& meta)
{
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

    for (size_t i = 0; i < num_ectxt; i++)
    {
        resCtxt[i].writeTo(css);
    }

    HostCtxtTemp = css.str();
    ectxt_len = HostCtxtTemp.length();
    ectxt = (uint8_t*) HostCtxtTemp.c_str();
    totalLength += ectxt_len;

    // cout << "Host: send/receive Ctxt to/from enclave:" << endl;
    result = multipleCtxtsTransform(enclave, &ret, ectxt, ectxt_len, num_ectxt, &octxt, &octxt_len);
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

    resCtxt.clear();
    resCtxt = vector<Ctxt>(num_ectxt, Ctxt(meta.data->publicKey));
    css.str(std::string());
    css.clear();
    EnclaveCtxtTemp = std::string(octxt, octxt + octxt_len);
    css << EnclaveCtxtTemp;
    for (size_t i = 0; i < num_ectxt; ++i)
    {
        resCtxt[i].Ctxt::read(css);
    }
    totalLength += octxt_len;

exit:
    cout << "Host: free memory for ectxt and octxt!" << endl;
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
    int num_channels = 4;
    int kernel_size = 7;
    int stride = 3;
    
    // Extract windows from 28 * 28 images (e.g., imgMat[0] ---> first image)
    int window_size = ceil((MNIST_HEIGHT - kernel_size) / stride) + 1; // 8
    int num_windows = window_size * window_size; // 64
    int num_kernels = kernel_size * kernel_size; // 49
    Mat<RR>* feaMat = new Mat<RR>[num_kernels]; // matrix for ct.I_{i,j}

    CropImages(test_imgs, test_lbls, feaMat, MNIST_HEIGHT, MNIST_WIDTH, num_imgs, kernel_size, stride, num_windows, window_size);
    
    /*---------------------------------------*/
    //  Load the model
    /*---------------------------------------*/
    cout << "Loading the model params..." << endl;
    int kernel_dim1, kernel_dim2;
    int dense1_dim1, dense1_dim2;
    int dense2_dim1, dense2_dim2;

    // 1. matrix for ct.K_{i, j}_{k}, num = 28 * 7
    Mat<RR>* kernel_weights;
    // 2. matrix for ct.W_{k}, num = 256/64 =4
    Mat<RR>* dense1_weights;
    // 3. matrix for ct.V, num = 1, pad zeros into matrix
    Mat<RR> dense2_weights;

    LoadModel(kernel_weights, dense1_weights, dense2_weights, kernel_dim1, kernel_dim2, dense1_dim1, dense1_dim2, dense2_dim1, dense2_dim2, num_windows, num_imgs);
    
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
    // ckksParams param(/*m=*/16 * 1024, /*bits=*/300, /*precision=*/20, /*c=*/2);
    // The below set can be used in E2DM baseline and comment the functions that invoke enclave.
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
    vector<Ctxt> ct_I(num_kernels, Ctxt(meta.data->publicKey));
    vector<Ctxt> ct_K(kernel_dim1 * kernel_dim2, Ctxt(meta.data->publicKey));
    vector<vector<Ctxt>> ct_W;
    vector<Ctxt> ct_V;
    vector<Ctxt> ct_Ck(num_channels, Ctxt(meta.data->publicKey));
    Ctxt ct_F(meta.data->publicKey);
    vector<EncodedPtxt> Initpoly;
    Ctxt inference_result(meta.data->publicKey);

    // record the communication cost and run time
    auto start= std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = 0.0;
    uint64_t totalLength = 0;

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
    // cout << "Host: check size of context and keys: " << ((double) context_len / (double)(1024 * 1024)) << " MB" << endl;

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
    CKKSmatrix.genMultBPoly(Initpoly);

    /*---------------------------------------*/
    //  Encrypt the images
    /*---------------------------------------*/
    cout << "Encrypting the images (feaMat, size = 7 * 7)..." << endl;
    for (int i = 0; i < num_kernels; i++) // 49
    {
        CKKSmatrix.encryptRmat(ct_I[i], feaMat[i]);
    }
    
    /*---------------------------------------*/
    //  Encrypt the model params
    /*---------------------------------------*/
    // 1. ct.K_{i, j}_{k}
    // cout << "Preparing model params ciphertext, ct_K..." << endl;
    for (int i = 0; i < ct_K.size(); i++) // 28 * 7
    {
        CKKSmatrix.encryptRmat(ct_K[i], kernel_weights[i]);
    }

    // 2. ct.W_{k}, preprocessed ctxt
    // cout << "Preparing model params ciphertext, ct_W..." << endl;
    for (int i = 0; i < dense1_dim2 / dense1_dim1; i++) // 256/64=4
    {
        vector<Ctxt> temp;
        CKKSmatrix.genInitActxt(temp, dense1_weights[i]);
        ct_W.push_back(temp);
    }

    // 3. ct.V, preprocessed ctxt
    // cout << "Preparing model params ciphertext, ct_V..." << endl;
    CKKSmatrix.genInitRecActxt(ct_V, dense2_weights);

    /*---------------------------------------*/
    //  Homomorphically perform CNN inference
    /*---------------------------------------*/
    start= chrono::steady_clock::now();
    // 1. HE conv layer, ct.I_{i, j} and ct.K_{i, j}_{k}
    // cout << endl << "Inference-convolution layer..." << endl;
    for (int k = 0; k < num_channels; k++) // 4
    {
        for (int i = 0; i < num_kernels; i++) // 49
        {
            Ctxt tmp = ct_I[i];
            tmp.multLowLvl(ct_K[k * num_kernels + i]);
            // tmp *= ct_K[k * num_kernels + i];
            ct_Ck[k] += tmp;
        }
        ct_Ck[k].reLinearize(); // return to a canonical state
    }

    // 2. HE square layer
    // cout << endl << "Inference-square1 layer..." << endl;
    for (int k = 0; k < num_channels; k++) // 4
    {
        ct_Ck[k].square();
        ct_Ck[k].bumpNoiseBound(1e-5);
    }

    // Todo: set a if() statement to check if invoking refreshCtxts()
    ret = RefreshMultipleCtxts(ct_Ck, totalLength, meta);
    if (ret != 0)
    {
        cerr << "Host: RefreshMultipleCtxts failed with " << ret << endl;
        goto exit;
    }

    // 3. HE FC-1 layer, 64 * 256 X 256 * 64 = 64 * 64, ct.W_{k} and ct_C_{k}
    // cout << endl << "Inference-FC-1 layer..." << endl;
    for (int k = 0; k < num_channels; k++) // 4
    {
        Ctxt temp(meta.data->publicKey);
        CKKSmatrix.HEmatmul_preprocessing(temp, ct_W[k], ct_Ck[k], Initpoly);
        ct_F += temp;
    }

    // cout << "ct_F.capacity=" << ct_F.capacity() << " ";
    // cout << "ct_F.isCorrect=" << ct_F.isCorrect() << " ";
    // cout << "ct_F.errorBound=" << ct_F.errorBound() << "\n";

    ret = RefreshSingleCtxt(ct_F, totalLength);
    if (ret != 0)
    {
        cerr << "Host: RefreshSingleCtxt failed with " << ret << endl;
        goto exit;
    }

    // 4. HE square layer
    // cout << endl << "Inference-square2 layer..." << endl;
    ct_F.square();

    // cout << endl << "ct_F.capacity=" << ct_F.capacity() << " ";
    // cout << "ct_F.isCorrect=" << ct_F.isCorrect() << " ";
    // cout << "ct_F.errorBound=" << ct_F.errorBound() << "\n";

    // 5. HE FC-2 layer, 10 * 64 X 64 * 64 = 10 * 64, ct.V and ct_F
    // cout << endl << "Inference-FC-2 layer..." << endl;
    CKKSmatrix.HErmatmul_preprocessing(inference_result, ct_V, ct_F, Initpoly);
    inference_result.bumpNoiseBound(1e-7);

    // cout << endl << "inference_result.isCorrect=" << inference_result.isCorrect() << endl;

    end = std::chrono::steady_clock::now();
    diff = end - start;
    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Host: HETEE Inference Time = " << timeElapsed << " s" << endl;
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
