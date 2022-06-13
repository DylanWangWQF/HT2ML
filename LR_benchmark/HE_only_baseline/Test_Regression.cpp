//
//  Test_Regression.cpp
//  ePPDSC
//
//  Created by Qifan Wang on 25/12/21.
//

#include <stdio.h>
#include "../src/helib/helib.h"
#include "Regression.h"
#include "Test_Regression.h"
#include <time.h>
#include <ctime>
#include <fstream>

void LRTest::RunRegressionTest(){
    // Load data
    Mat<long> rawData;
    vector<long> labels;
    long dim;
    string datafile = "/Users/qifanwang/code/MyExp/ePPDSC/ePPDSC/scripts/dataset.dat";
    
    if (!LoadData(rawData, labels, dim, datafile))
        return;
    
    //initialize HE parameters
    LRparams param(/*m =*/8471,
                      /*p =*/127,
                      /*r =*/1,
                      /*bits =*/200,
                      /*c=*/2,
                      /*gens =*/std::vector<long>{3744, 2366},
                      /*ords =*/std::vector<long>{42, -2},
                      /*mvec =*/std::vector<long>{43, 197});
    LRmeta meta;
    meta(param);
    meta.data->context.printout();
    Regression regression(meta);
    
    // Batch data samples
    cout << endl << "Begin to batch data!" << endl;
    vector<vector<vector<long>>> ptxtData;
    vector<vector<long>> ptxtLabels;
    BatchData(ptxtData, ptxtLabels, rawData, labels, param.p, meta.data->ea.size());
    
    // Encrypt data
    cout << endl << "Begin to encrypt data!" << endl;
    regression.EncryptData(ptxtData, ptxtLabels);
//    cout << "Check data size: " << regression.data.size() << endl;
//    cout << "Check labels size: " << regression.labels.size() << endl;
    
    // Linear Regression
    cout << endl << "Begin to perform linear regression!" << endl;
    vector<Ctxt> theta;
    Ctxt det(meta.data->publicKey);
    regression.Regress(theta, det);
    
    // Decrypt the result and check
    cout << endl << "Computed values: " << endl;
    for (unsigned i = 0; i < theta.size(); i++) {
        vector<long> theta_vec;
        meta.data->ea.decrypt(theta[i], meta.data->secretKey, theta_vec);
        cout << "  theta[" << i << "] = " << theta_vec[0] << endl;
    }
    
    vector<long> det_dec;
    meta.data->ea.decrypt(det, meta.data->secretKey, det_dec);
    cout << "  Determinant: " << det_dec[0] << endl << endl;
    
    
    // Plaintext training
    vector<long> theta_Pt;
    long det_Pt;
    RegressPT(theta_Pt, det_Pt, rawData, labels);
    cout << "Expected values: " << endl;
    for (unsigned i = 0; i < theta.size(); i++)
    {
        cout << "  theta[" << i << "] = " << theta_Pt[i] % (param.p) << endl;
    }
    cout << "  Determinant: " << det_Pt % (param.p) << endl;
    cout << endl;
}
