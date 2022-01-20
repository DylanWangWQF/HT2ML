#include "MatOpe.h"

void ComponentProduct(Mat<RR> &mat1, Mat<RR> &mat2, Mat<RR> &mat3)
{
    mat1.SetDims(mat2.NumRows(), mat2.NumCols());
    for (long i = 0; i < mat2.NumRows(); i++)
    {
        for (long j = 0; j < mat2.NumCols(); j++)
        {
            mat1[i][j] = mat2[i][j] * mat3[i][j];
        }
    }
    
}