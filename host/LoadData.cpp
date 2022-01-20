#include <fstream>
#include <string>

#include "LoadData.h"
#include "data/DatasetOperations.h"

bool LoadData(Mat<RR> &rawData, int &weight_dim1, int &weight_dim2, string &filename){
    ifstream fin;
    fin.open(filename);
    if (!fin) {
        cout << "Unable to read data file." << endl;
        return false;
    }
    
    fin >> weight_dim1 >> weight_dim2;
    double data;

    rawData.SetDims(weight_dim1, weight_dim2);
    for (long i = 0; i < weight_dim1; i++) {
        for (long j = 0; j < weight_dim2; j++) {
            fin >> data;
            rawData[i][j] = to_RR(data);
        }
    }

    return true;
}

void CropImages(vector<vector<float>>& test_imgs, vector<unsigned char>& test_lbls, Mat<RR>*& feaMat, int MNIST_HEIGHT, int MNIST_WIDTH, int num_imgs, int kernel_size, int stride, int num_windows, int window_size)
{
    size_t test_img_limit = 0;
    string datasets_dir = "/home/dylan/mysamples/TEE-HE-PPML/scripts";
    test_imgs = loadMnistTestImages(datasets_dir, test_img_limit);
    test_lbls = loadMnistTestLabels(datasets_dir, test_img_limit);
    // cout << "Dimension of imgs: (" << test_imgs.size() << ", " << test_imgs[0].size() << ")" << endl << endl;
    // cout << "Print first image: " << endl;
    // for (int i = 0; i < test_imgs[0].size(); i++) {
    //     cout <<  test_imgs[0][i] << " ";
    // }
    // cout << endl << "Number of labels for test: " << test_lbls.size() << endl << endl;

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
    
}

void LoadModel(Mat<RR>*& kernel_weights, Mat<RR>*& dense1_weights, Mat<RR>& dense2_weights, int& kernel_dim1, int& kernel_dim2, int& dense1_dim1, int& dense1_dim2, int& dense2_dim1, int& dense2_dim2, int num_windows, int num_imgs)
{
    Mat<RR> raw_kernel_weights; // 28 * 7
    Mat<RR> raw_dense1_weights; // 64 * 256
    Mat<RR> raw_dense2_weights; // 10* 64

    string datafile1 = "/home/dylan/mysamples/TEE-HE-PPML/host/model/kernels_weights.dat";
    string datafile2 = "/home/dylan/mysamples/TEE-HE-PPML/host/model/dense1_weights.dat";
    string datafile3 = "/home/dylan/mysamples/TEE-HE-PPML/host/model/dense2_weights.dat";
    if (!LoadData(raw_kernel_weights, kernel_dim1, kernel_dim2, datafile1)) {
        return;
    }
    // cout << "Generated kernel weights matrix raw_kernel_weights: " << endl << raw_kernel_weights << endl;
    if (!LoadData(raw_dense1_weights, dense1_dim1, dense1_dim2, datafile2)) {
        return;
    }
    if (!LoadData(raw_dense2_weights, dense2_dim1, dense2_dim2, datafile3)) {
        return;
    }

    kernel_weights = new Mat<RR>[kernel_dim1 * kernel_dim2];
    dense1_weights = new Mat<RR>[dense1_dim2 / dense1_dim1];

    // Crop the weights matrix according to E2DM
    // 1. matrix for ct.K_{i, j}_{k}, num = 28 * 7
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

}