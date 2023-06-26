#ifndef LoadData_h
#define LoadData_h

#include <vector>
#include <NTL/RR.h>
#include <NTL/mat_RR.h>
#include <NTL/matrix.h>

using namespace std;
using namespace NTL;

/*
 Load Raw data and Transform into Matrix
 */
bool LoadData(Mat<RR> &rawData, int &weight_dim1, int &weight_dim2, string &filename);

void CropImages(vector<vector<float>>& test_imgs, vector<unsigned char>& test_lbls, Mat<RR>*& feaMat, int MNIST_HEIGHT, int MNIST_WIDTH, int num_imgs, int kernel_size, int stride, int window_size);

void LoadModel(Mat<RR>*& kernel_weights, Mat<RR>*& dense1_weights, int& kernel_dim1, int& kernel_dim2, int& dense1_dim1, int& dense1_dim2, int num_imgs);

#endif /* LoadData_h */