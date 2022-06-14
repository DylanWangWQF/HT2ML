//
//  Test_Regression.cpp
//  ePPDSC
//
//  Created by Qifan Wang on 25/12/21.
//

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sys/time.h>
#include <chrono>
#include <NTL/BasicThreadPool.h>

#include "Regression.h"

int main(int argc, const char * argv[]) {
    SetNumThreads(16);

    // Load data
    Mat<long> rawData;
    vector<long> labels;
    long dim;
    // string datafile = "/home/dylan/code/EHPPML/LR_benchmark/scripts/dataset.dat";
    string datafile = "/home/dylan/code/EHPPML/LR_benchmark/scripts/6_dim_LR.dat";
    
    if (!LoadData(rawData, labels, dim, datafile))
        return 0;
    
    //initialize HE parameters
    // LRparams param(/*m =*/8471,
    //                   /*p =*/127,
    //                   /*r =*/1,
    //                   /*bits =*/200,
    //                   /*c=*/2,
    //                   /*gens =*/std::vector<long>{3744, 2366},
    //                   /*ords =*/std::vector<long>{42, -2},
    //                   /*mvec =*/std::vector<long>{43, 197});
    LRparams param(/*m =*/17507,
                      /*p =*/127,
                      /*r =*/1,
                      /*bits =*/350,
                      /*c=*/2,
                      /*gens =*/std::vector<long>{9395, 2502, 12629},
                      /*ords =*/std::vector<long>{40, 6, -2},
                      /*mvec =*/std::vector<long>{41, 7, 61}); //nslots = 480

    LRmeta meta;
    meta(param);
    meta.data->context.printout();
    Regression regression(meta);
    
    // Batch data samples
    cout << endl << "Begin to batch data!" << endl;
    vector<vector<vector<long>>> ptxtData;
    vector<vector<long>> ptxtLabels;
    BatchData(ptxtData, ptxtLabels, rawData, labels, param.p, meta.data->ea.size());

    // cout << "Check the batch data:" << endl;
    // for (int i = 0; i < ptxtData.size(); i++)
    // {
    //     cout << "1st dimension of i = " << i << endl;
    //     for (int j = 0; j < ptxtData[0].size(); j++)
    //     {
    //         for (int k = 0; k < ptxtData[0][0].size(); k++)
    //         {
    //             cout << ptxtData[i][j][k] << " ";
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    // }
    
    // Encrypt data
    cout << endl << "Begin to encrypt data!" << endl;
    regression.EncryptData(ptxtData, ptxtLabels);
    // cout << "Check data size: " << regression.data.size() << endl;
    // cout << "Check labels size: " << regression.labels.size() << endl;
    
    // Linear Regression
    cout << endl << "Linear regression in HE (1-4): " << endl;
    auto start= std::chrono::steady_clock::now();

    vector<Ctxt> theta;
    Ctxt det(meta.data->publicKey);
    regression.Regress(theta, det);

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    double timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
    
    // Decrypt the result and check
    cout << endl << "4. Expected theta in HE version: " << endl;
    for (unsigned i = 0; i < theta.size(); i++) {
        vector<long> theta_vec;
        meta.data->ea.decrypt(theta[i], meta.data->secretKey, theta_vec);
        cout << "  theta[" << i << "] = " << theta_vec[0] << endl;
    }
    
    vector<long> det_dec;
    meta.data->ea.decrypt(det, meta.data->secretKey, det_dec);
    cout << "Determinant: " << det_dec[0] << endl << endl;

    cout << "--------------------------------------------------" << endl;
    
    // Plaintext training
    vector<long> theta_Pt;
    long det_Pt;
    cout << endl << "Linear regression in plaintext (1-4): " << endl;
    RegressPT(theta_Pt, det_Pt, rawData, labels);
    cout << endl << "4. Expected theta in plaintext version: " << endl;
    for (unsigned i = 0; i < theta.size(); i++)
    {
        cout << "  theta[" << i << "] = " << theta_Pt[i] % (param.p) << endl;
    }
    cout << "Determinant: " << det_Pt % (param.p) << endl;

    cout << "------------------------------------------------------------------------" << endl;
    cout << "Total HE Training Time = " << timeElapsed << " s" << endl;
    cout << "------------------------------------------------------------------------" << endl;

    return 0;
}