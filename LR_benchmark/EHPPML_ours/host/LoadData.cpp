#include <fstream>
#include <string>

#include "LoadData.h"
#include <NTL/BasicThreadPool.h>

// bool LoadData(Mat<long> &rawData, Mat<long> &rawLabel, long &MatrixDim, string &filename){
//     ifstream fin;
//     fin.open(filename);
//     if (!fin) {
//         cout << "Unable to read data file." << endl;
//         return false;
//     }
    
//     long numAttr; //number of attributes
//     long n; //the number of data samples
//     long data;
//     fin >> numAttr >> n;
    
//     // MatrixDim should be 16, 32, 64
//     rawData.SetDims(n, MatrixDim); // data sample matrix
//     // MatInit(rawData);
//     rawLabel.SetDims(n, MatrixDim); // data sample matrix
//     // MatInit(rawLabel);
//     for (long i = 0; i < n; i++) {
//         for (long j = 0; j < numAttr; j++) {
//             fin >> rawData[i][j];
//         }
//         fin >> rawLabel[i][0];
//     }
//     return true;
// }

void ProcessDataMatrix(Mat<long>*& Amat, Mat<long>*& ATranmat, Mat<long>*& Bmat, long& numMat, long& MatrixDim, string& filename){
    ifstream fin;
    fin.open(filename);
    if (!fin) {
        cout << "Unable to read data file." << endl;
        return;
    }
    
    long numAttr; //number of attributes
    long n; //the number of data samples
    long data;
    fin >> numAttr >> n;
    
    // MatrixDim should be 16, 32, 64
    Mat<long> rawData;
    rawData.SetDims(n, MatrixDim); // data sample matrix
    // MatInit(rawData);
    Mat<long> rawLabel;
    rawLabel.SetDims(n, MatrixDim); // data sample matrix
    // MatInit(rawLabel);
    for (long i = 0; i < n; i++) {
        for (long j = 0; j < numAttr; j++) {
            fin >> rawData[i][j];
        }
        fin >> rawLabel[i][0];
    }

    long nrows = MatrixDim;
    long ncols = nrows;
    
    if (rawData.NumRows() % nrows != 0)
    {
        numMat = rawData.NumRows()/nrows + 1;
        // cout << "Number of SubMatrix: " << numMat << endl;
        Amat = new Mat<long>[numMat];
        Bmat = new Mat<long>[numMat];

        NTL_EXEC_RANGE(numMat - 1, first, last);
        for(long k = first; k < last; ++k){
            Amat[k].SetDims(nrows, ncols);
            Bmat[k].SetDims(nrows, ncols);
            //MatInit(Amat[k]);
            for(long i = 0; i < nrows; i++){
                for(long j = 0; j < ncols; j++){
                    Amat[k][i][j] = rawData[k * nrows + i][j];
                    Bmat[k][i][j] = rawLabel[k * nrows + i][j];
                }
            }
        }
        NTL_EXEC_RANGE_END;

        // The last matrix
        Amat[numMat - 1].SetDims(nrows, ncols);
        Bmat[numMat - 1].SetDims(nrows, ncols);
        //MatInit(Amat[numMat - 1]);
        for(long i = 0; i < nrows ; i++){
            for(long j = 0; j < ncols; j++){
                if (((numMat - 1) * nrows + i) < rawData.NumRows()){
                    Amat[numMat - 1][i][j] = rawData[(numMat - 1) * nrows + i][j];
                    Bmat[numMat - 1][i][j] = rawLabel[(numMat - 1) * nrows + i][j];
                }else{
                    Amat[numMat - 1][i][j] = 0;
                    Bmat[numMat - 1][i][j] = 0;
                }
            }
        }
    }else{
        //In this case, last matrix doesn't need to pad with zero
        numMat = rawData.NumRows()/nrows;
        Amat = new Mat<long>[numMat];
        Bmat = new Mat<long>[numMat];
        NTL_EXEC_RANGE(numMat, first, last);
        for(long k = first; k < last; ++k){
            Amat[k].SetDims(nrows, ncols);
            Bmat[k].SetDims(nrows, ncols);
            //MatInit(Amat[k]);
            for(long i = 0; i < nrows ; i++){
                for(long j = 0; j < ncols; j++){
                    Amat[k][i][j] = rawData[k * nrows + i][j];
                    Bmat[k][i][j] = rawLabel[k * nrows + i][j];
                }
            }
        }
        NTL_EXEC_RANGE_END;
    }

    // get the transpose of data sample matrix
    ATranmat = new Mat<long>[numMat];
    NTL_EXEC_RANGE(numMat, first, last);
    for(long k = first; k < last; ++k){
        ATranmat[k].SetDims(nrows, ncols);
        //MatInit(ATranmat[k]);
        for(long i = 0; i < nrows ; i++){
            for(long j = 0; j < ncols; j++){
                ATranmat[k][i][j] = Amat[k][j][i];
            }
        }
    }
    NTL_EXEC_RANGE_END;
}

void MatInit(Mat<long>& matInit)
{
    for (int i = 0; i < matInit.NumCols(); ++i){
        for (int j = 0; j < matInit.NumRows(); ++j){
            matInit[j][i] = 0;
        }
    }
}