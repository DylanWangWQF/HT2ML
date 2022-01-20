#include <fstream>
#include <string>

#include "LoadData.h"


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