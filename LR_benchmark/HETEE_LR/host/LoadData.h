//
//  LoadData.h
//  mHT
//
//  Created by Qifan Wang on 27/05/21.
//

#ifndef LoadData_h
#define LoadData_h

#include <vector>
#include <NTL/matrix.h>

using namespace std;
using namespace NTL;

// Load Raw data and Transform into Matrix
// bool LoadData(Mat<long> &rawData, Mat<long> &rawLabel, long &MatrixDim, string &filename);
// Initialize Mat with 0
void MatInit(Mat<long>& matInit);
// Process raw Dataset, sub Mat with dimension
void ProcessDataMatrix(Mat<long>*& Amat, Mat<long>*& ATranmat, Mat<long>*& Bmat, long& numMat, long& MatrixDim, string &filename);

#endif /* LoadData_h */