//
//  LoadData.h
//  mHT
//
//  Created by Qifan Wang on 27/05/21.
//

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

#endif /* LoadData_h */