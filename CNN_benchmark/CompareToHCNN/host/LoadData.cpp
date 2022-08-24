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

void CropImages(vector<vector<float>>& test_imgs, vector<unsigned char>& test_lbls, Mat<RR>*& feaMat, int MNIST_HEIGHT, int MNIST_WIDTH, int num_imgs, int kernel_size, int stride, int window_size)
{
    size_t test_img_limit = 0;
    string datasets_dir = "/home/dylan/code/HETEE/CNN_benchmark/CompareToHCNN/scripts";
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
    for (temp_i = 0; temp_i < kernel_size; ++temp_i) // 5 * 5
    {
        for (temp_j = 0; temp_j < kernel_size; ++temp_j)
        {
            feaMat[temp_idx1].SetDims(576, num_imgs); // 576 * 64
            for(temp_y = 0; temp_y < num_imgs; ++temp_y) // from 0-th image to 63-th image
            {
                temp_idx2 = 0;
                for(temp_x = 0; temp_x < 576 ; ++temp_x) // from 0-th window (0, 0) to 576-th window (23, 23), window_size = 24
                {
                    feaMat[temp_idx1][temp_x][temp_y] = imgMat[temp_y][(temp_x / window_size) * stride + temp_i][(temp_x % window_size) * stride + temp_j];
                }
            }
            temp_idx1++;
        }
    }
    // cout << "Generated feature matrix feaMat: " << endl << feaMat[0][575][0] << endl;
    
}

void LoadModel(Mat<RR>*& kernel_weights, Mat<RR>*& dense1_weights, int& kernel_dim1, int& kernel_dim2, int& dense1_dim1, int& dense1_dim2, int num_imgs)
{
    Mat<RR> raw_kernel_weights; // 30 * 5
    Mat<RR> raw_dense1_weights; // 10 * 864

    string datafile1 = "/home/dylan/code/HETEE/CNN_benchmark/CompareToHCNN/host/model/kernels_weights.dat";
    string datafile2 = "/home/dylan/code/HETEE/CNN_benchmark/CompareToHCNN/host/model/dense1_weights.dat";
    if (!LoadData(raw_kernel_weights, kernel_dim1, kernel_dim2, datafile1)) {
        return;
    }
    // cout << "Generated kernel weights matrix raw_kernel_weights: " << endl << raw_kernel_weights << endl;
    if (!LoadData(raw_dense1_weights, dense1_dim1, dense1_dim2, datafile2)) {
        return;
    }

    kernel_weights = new Mat<RR>[kernel_dim1 * kernel_dim2];
    dense1_weights = new Mat<RR>[14]; // 864 / 64 = 13.5

    // Crop the weights matrix according to E2DM
    // 1. matrix for ct.K_{i, j}_{k}, num = 30 * 5
    int temp_idx = 0;
    for (int m = 0; m < kernel_dim1; m++) // 30
    {
        for (int n = 0; n < kernel_dim2; n++) // 5
        {
            kernel_weights[temp_idx].SetDims(64, num_imgs);
            for (int i = 0; i < 64; i++) // 64
            {
                for (int j = 0; j < num_imgs; j++) // 64
                {
                    kernel_weights[temp_idx][i][j] = raw_kernel_weights[m][n];
                }
            }
            temp_idx++;
        }
    }
    // cout << "Generated kernel_weights matrix kernel_weights[0]: " << endl << kernel_weights[0] << endl;

    // 2. matrix for ct.W_{k}, num = 14
    for (int k = 0; k < 14; k++)
    {
        // original-HE
        // dense1_weights[k].SetDims(64, 64);
        // optimised-HE
        dense1_weights[k].SetDims(16, 64);
        for (int i = 0; i < 16; i++)
        {
            for (int j = 0; j < 64; j++)
            {
                int idx = k * 64 + j;
                if ((i < 10) && (idx < 864))
                {
                    dense1_weights[k][i][j] = raw_dense1_weights[i][idx];
                }
                else
                {
                    dense1_weights[k][i][j] = to_RR(0);
                }
            }
        }
    }
    // cout << endl << "Generated dense1_weights matrix dense1_weights[0]: " << endl << dense1_weights[0] << endl;
    // cout << endl << "Generated dense1_weights matrix dense1_weights[1]: " << endl << dense1_weights[1] << endl;

}