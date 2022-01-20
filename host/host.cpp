#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <sys/time.h>
#include <chrono>

// libs
#include <NTL/BasicThreadPool.h>
#include <NTL/matrix.h>
#include <openenclave/host.h>
#include <helib/helib.h>

#include "data/DatasetOperations.h"
#include "HEMat/CKKSmatrix.h"
#include "utils/LoadData.h"
#include "utils/MatOpe.h"
#include "teehe_u.h"

oe_enclave_t* enclave = NULL;
// constexpr int MNIST_CHANNELS = 1;
constexpr int MNIST_HEIGHT = 28;
constexpr int MNIST_WIDTH = 28;
// constexpr int MNIST_LABELS = 10;

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

    cout << "Host: send/receive Ctxt to/from enclave:" << endl;
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

int RefreshMultipleCtxts(vector<Ctxt>& resCtxt, uint64_t& totalLength)
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
    size_t num_octxt = 0;

    for (size_t i = 0; i < num_ectxt; i++)
    {
        resCtxt[i].writeTo(css);
    }

    HostCtxtTemp = css.str();
    ectxt_len = HostCtxtTemp.length();
    ectxt = (uint8_t*) HostCtxtTemp.c_str();
    totalLength += ectxt_len;

    cout << "Host: send/receive Ctxt to/from enclave:" << endl;
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

    css.str(std::string());
    css.clear();
    EnclaveCtxtTemp = std::string(octxt, octxt + octxt_len);
    css << EnclaveCtxtTemp;
    for (size_t i = 0; i < num_octxt; ++i)
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
    SetNumThreads(4); // set threads for parallel operations when using NTL

    /*---------------------------------------*/
    //  Load the dataset
    /*---------------------------------------*/
    cout << "Loading test imgs & labels..." << endl;
    size_t test_img_limit = 0;
    vector<vector<float>> test_imgs;
    vector<unsigned char> test_lbls;
    string datasets_dir = "/home/dylan/code/TEE-HE-PPML/scripts";
    test_imgs = loadMnistTestImages(datasets_dir, test_img_limit);
    test_lbls = loadMnistTestLabels(datasets_dir, test_img_limit);
    // cout << "Dimension of imgs: (" << test_imgs.size() << ", " << test_imgs[0].size() << ")" << endl << endl;
    // cout << "Print first image: " << endl;
    // for (int i = 0; i < test_imgs[0].size(); i++) {
    //     cout <<  test_imgs[0][i] << " ";
    // }
    // cout << endl << "Number of labels for test: " << test_lbls.size() << endl << endl;

    /*---------------------------------------*/
    //  Crop the images
    /*---------------------------------------*/
    int num_imgs = 64;
    int num_channels = 4;
    int kernel_size = 7;
    int stride = 3;
    Mat<RR>* imgMat = new Mat<RR>[num_imgs]; // 64 * 28 * 28
    for (int k = 0; k < num_imgs; k++)
    {
        imgMat[k].SetDims(MNIST_HEIGHT, MNIST_WIDTH);
        int p_loc = 0;
        for(long i = 0; i < MNIST_HEIGHT ; i++)
        {
            for(long j = 0; j < MNIST_WIDTH; j++)
            {
                imgMat[k][i][j] = to_RR(test_imgs[k][p_loc]);
                // imgMat[k][i][j] = to_RR( (int)(100.0 * test_imgs[k][p_loc] + 0.5) / 100.0); // two pos after point
                p_loc++;
            }
        }
    }
    // cout << "Generated image matrix imgMat: " << endl << imgMat[0] << endl;

    // Extract windows from 28 * 28 images (e.g., imgMat[0] ---> first image)
    int window_size = ceil((MNIST_HEIGHT - kernel_size) / stride) + 1; // 8
    int num_windows = window_size * window_size; // 64
    int num_kernels = kernel_size * kernel_size; // 49
    Mat<RR>* feaMat = new Mat<RR>[num_kernels]; // matrix for ct.I_{i,j}

    int temp_i, temp_j, temp_x, temp_y, temp_idx2;
    int temp_idx1 = 0;
    for (temp_i = 0; temp_i < kernel_size; ++temp_i) // 7 * 7
    {
        for (temp_j = 0; temp_j < kernel_size; ++temp_j)
        {
            feaMat[temp_idx1].SetDims(num_windows, num_imgs); // 64 * 64
            for(temp_y = 0; temp_y < num_imgs; ++temp_y) // from 0-th image to 63-th image
            {
                temp_idx2 = 0;
                for(temp_x = 0; temp_x < num_windows ; ++temp_x) // from 0-th window (0, 0) to 63-th window (7, 7)
                {
                    feaMat[temp_idx1][temp_x][temp_y] = imgMat[temp_y][(temp_x / window_size) * stride + temp_i][(temp_x % window_size) * stride + temp_j];
                }
            }
            temp_idx1++;
        }
    }
    // cout << "Generated feature matrix feaMat: " << endl << feaMat[0] << endl;
    
    /*---------------------------------------*/
    //  Load the model
    /*---------------------------------------*/
    cout << "Loading the model params..." << endl;
    Mat<RR> raw_kernel_weights; // 28 * 7
    int kernel_dim1, kernel_dim2;
    Mat<RR> raw_dense1_weights; // 64 * 256
    int dense1_dim1, dense1_dim2;
    Mat<RR> raw_dense2_weights; // 10* 64
    int dense2_dim1, dense2_dim2;
    string datafile1 = "/home/dylan/code/TEE-HE-PPML/host/model/kernels_weights.dat";
    string datafile2 = "/home/dylan/code/TEE-HE-PPML/host/model/dense1_weights.dat";
    string datafile3 = "/home/dylan/code/TEE-HE-PPML/host/model/dense2_weights.dat";
    if (!LoadData(raw_kernel_weights, kernel_dim1, kernel_dim2, datafile1)) {
        return 0;
    }
    // cout << "Generated kernel weights matrix raw_kernel_weights: " << endl << raw_kernel_weights << endl;
    if (!LoadData(raw_dense1_weights, dense1_dim1, dense1_dim2, datafile2)) {
        return 0;
    }
    if (!LoadData(raw_dense2_weights, dense2_dim1, dense2_dim2, datafile3)) {
        return 0;
    }

    // Crop the weights matrix according to E2DM
    // 1. matrix for ct.K_{i, j}_{k}, num = 28 * 7
    Mat<RR>* kernel_weights = new Mat<RR>[kernel_dim1 * kernel_dim2];
    int temp_idx = 0;
    for (int m = 0; m < kernel_dim1; m++) // 28
    {
        for (int n = 0; n < kernel_dim2; n++) // 7
        {
            kernel_weights[temp_idx].SetDims(num_windows, num_imgs);
            for (int i = 0; i < num_windows; i++) // 64
            {
                for (int j = 0; j < num_imgs; j++) // 64
                {
                    kernel_weights[temp_idx][i][j] = raw_kernel_weights[m][n];
                    // kernel_weights[temp_idx][i][j] = to_RR((int)(100.0 * to_double(raw_kernel_weights[m][n]) + 0.5) / 100.0);
                }
            }
            temp_idx++;
        }
    }
    // cout << "Generated kernel_weights matrix kernel_weights[0]: " << endl << kernel_weights[0] << endl;

    // 2. matrix for ct.W_{k}, num = 256/64 =4
    Mat<RR>* dense1_weights = new Mat<RR>[dense1_dim2 / dense1_dim1];
    for (int k = 0; k < dense1_dim2 / dense1_dim1; k++)
    {
        dense1_weights[k].SetDims(dense1_dim1, dense1_dim1);
        for (int i = 0; i < dense1_dim1; i++) // 64
        {
            for (int j = 0; j < dense1_dim1; j++) // 256
            {
                dense1_weights[k][i][j] = raw_dense1_weights[i][k * dense1_dim1 + j];
                // dense1_weights[k][i][j] = to_RR((int)(100.0 * to_double(raw_dense1_weights[i][k * dense1_dim1 + j]) + 0.5) / 100.0);
            }
        }
    }
    // cout << "Generated dense1_weights matrix dense1_weights[0]: " << endl << dense1_weights[0] << endl;
    
    // 3. matrix for ct.V, num = 1, pad zeros into matrix
    Mat<RR> dense2_weights;
    dense2_weights.SetDims(16, dense2_dim2); // 16 * 64
    for (int i = 0; i < 16; i++) // pad 16-10=6 rows with zeros
    {
        for (int j = 0; j < dense2_dim2; j++) // 64
        {
            if (i < 10)
            {
                dense2_weights[i][j] = raw_dense2_weights[i][j];
            }
            else
            {
                dense2_weights[i][j] = to_RR(0);
            }
        }
    }
    // cout << "Generated dense2_weights matrix dense2_weights: " << endl << dense2_weights << endl;

    /*-----------------------------------------------Load images and models end!----------------------------------------------------------------------*/

    /*---------------------------------------*/
    //  Setup and Initialize HE params
    /*---------------------------------------*/
    // Setup the encalve params
    oe_result_t result;
    int ret = 0;
    const uint32_t flags = OE_ENCLAVE_FLAG_DEBUG_AUTO;

    ckksMatpar ckksmatpar;
    long ncols = 64, nrows = 64, subdim = 16;
    readckksMatpar(ckksmatpar, nrows, ncols, subdim); // subdim used in below RecMul
    ckksParams param(/*m=*/16 * 1024, /*bits=*/179, /*precision=*/20, /*c=*/2);
    ckksMeta meta;
    meta(param);
    cout << "HE Context contents: " << endl;
    meta.data->context.printout();
    CKKSmatrix CKKSmatrix(ckksmatpar, meta);

    //setup stringstream, record size and time
    stringstream ss;
    uint64_t totalLength = 0;

    // Setup the passed context, pk and sk, ctxts
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

    /*---------------------------------------*/
    //  Create the enclave
    //  Send context and keys to enclave
    /*---------------------------------------*/
    // Command Format: ./host/teehe_host ./enclave/enclave.signed
    cout << "Host: create enclave for image:" << argv[1] << endl;
    result = oe_create_teehe_enclave(argv[1], OE_ENCLAVE_TYPE_SGX, flags, NULL, 0, &enclave);
    if (result != OE_OK)
    {
        cerr << "oe_create_teehe_enclave() failed with " << argv[0] << " " << result << endl;
        ret = 1;
        goto exit;
    }

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
    cout << "Preparing model params ciphertext, ct_K..." << endl;
    for (int i = 0; i < ct_K.size(); i++) // 28 * 7
    {
        CKKSmatrix.encryptRmat(ct_K[i], kernel_weights[i]);
    }

    // 2. ct.W_{k}, preprocessed ctxt
    cout << "Preparing model params ciphertext, ct_W..." << endl;
    for (int i = 0; i < dense1_dim2 / dense1_dim1; i++) // 256/64=4
    {
        vector<Ctxt> temp;
        CKKSmatrix.genInitActxt(temp, dense1_weights[i]);
        ct_W.push_back(temp);
    }

    // 3. ct.V, preprocessed ctxt
    cout << "Preparing model params ciphertext, ct_V..." << endl;
    CKKSmatrix.genInitRecActxt(ct_V, dense2_weights);

    /*---------------------------------------*/
    //  Homomorphically perform CNN inference
    /*---------------------------------------*/
    // 1. HE conv layer, ct.I_{i, j} and ct.K_{i, j}_{k}
    cout << "Inference-convolution layer..." << endl;
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
    cout << "Inference-square1 layer..." << endl;
    for (int k = 0; k < num_channels; k++) // 4
    {
        ct_Ck[k].square();
    }

    // Todo: set a if() statement to check if invoking refreshCtxts()
    ret = RefreshMultipleCtxts(ct_Ck, totalLength);
    if (ret != 0)
    {
        cerr << "Host: RefreshMultipleCtxts failed with " << ret << endl;
        goto exit;
    }

    // 3. HE FC-1 layer, 64 * 256 X 256 * 64 = 64 * 64, ct.W_{k} and ct_C_{k}
    cout << "Inference-FC-1 layer..." << endl;
    for (int k = 0; k < num_channels; k++) // 4
    {
        Ctxt temp(meta.data->publicKey);
        CKKSmatrix.HEmatmul_preprocessing(temp, ct_W[k], ct_Ck[k], Initpoly);
        ct_F += temp;
    }

    ret = RefreshSingleCtxt(ct_F, totalLength);
    if (ret != 0)
    {
        cerr << "Host: RefreshSingleCtxt failed with " << ret << endl;
        goto exit;
    }

    // 4. HE square layer
    cout << "Inference-square2 layer..." << endl;
    ct_F.square();

    // 5. HE FC-2 layer, 10 * 64 X 64 * 64 = 10 * 64, ct.V and ct_F
    cout << "Inference-FC-2 layer..." << endl;
    CKKSmatrix.HErmatmul_preprocessing(inference_result, ct_V, ct_F, Initpoly);

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
