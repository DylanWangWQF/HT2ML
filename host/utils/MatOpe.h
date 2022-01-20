#ifndef MatOpe_h
#define MatOpe_h

#include <vector>
#include <NTL/RR.h>
#include <NTL/mat_RR.h>
#include <NTL/matrix.h>

using namespace std;
using namespace NTL;

void ComponentProduct(Mat<RR> &mat1, Mat<RR> &mat2, Mat<RR> &mat3);

#endif /* MatOpe_h */