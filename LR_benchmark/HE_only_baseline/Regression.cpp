//
//  Regression.cpp
//  ePPDSC
//
//  Created by Qifan Wang on 25/12/21.
//

#include <stdio.h>
#include "Regression.h"

bool LoadData(Mat<long> &rawData, vector<long> &labels, long &dim, string &filename) {
    
    ifstream fin;
    fin.open(filename);
    if (!fin) {
        cout << "Unable to read data file." << endl;
        return false;
    }
  
//    rawData.kill();
//    labels.clear();

    long label, n;

    fin >> dim >> n;
    rawData.SetDims(n, dim);
    for (long i = 0; i < n; i++)
    {
        for (long j = 0; j < dim; j++)
        {
            fin >> rawData[i][j];
        }
        fin >> label;
        labels.push_back(label);
    }  
    return true;
}

void BatchData(vector<vector<vector<long>>> &ptxtData, vector<vector<long>> &ptxtLabels, const Mat<long> &rawData, const vector<long> &labels, long p, long nslots) {
    ptxtData.clear();
    ptxtLabels.clear();
    for (long i = 0; i < rawData.NumRows(); i += nslots) // number of data
    {
        vector<vector<long>> row;
        vector<long> curLabel;
        for (long j = 0; j < rawData.NumCols(); j++) // dimension
        {
            vector<long> columnBatch;
            for (long k = i; ((k < i + nslots) && (k < rawData.NumRows())); k++)
            {
                columnBatch.push_back(rawData[k][j] % p);
                if(j == 0) curLabel.push_back(labels[k] % p);
            }
            // pad zero to columnBatch
            if (columnBatch.size() < nslots)
            {
//                cout << j << "-th column's size before pad = " << columnBatch.size() << endl;
//                for (long dummy = 0; dummy < (nslots - columnBatch.size()); dummy++)
//                {
//                    columnBatch.push_back(0);
//                    if(j == 0) curLabel.push_back(0);
//                }
                columnBatch.resize(nslots, 0);
                if(j == 0) curLabel.resize(nslots, 0);
//                cout << j << "-th column's size after pad = " << columnBatch.size() << endl;
            }
            row.push_back(columnBatch);
        }
        ptxtData.push_back(row);
        ptxtLabels.push_back(curLabel);
    }
    return;
}

void Regression::EncryptData(const vector<vector<vector<long>>> &ptxtData, const vector<vector<long>> &ptxtLabels){
    for (unsigned i = 0; i < ptxtData.size(); i++)
    {
        vector<Ctxt> encExample(ptxtData[i].size(), Ctxt(meta.data->publicKey));
        for (unsigned j = 0; j < ptxtData[i].size(); j++) {
            meta.data->ea.encrypt(encExample[j], meta.data->publicKey, ptxtData[i][j]);
        }
        Ctxt encLabel(meta.data->publicKey);
        meta.data->ea.encrypt(encLabel, meta.data->publicKey, ptxtLabels[i]);
        
        data.push_back(encExample);
        labels.push_back(encLabel);
    }
}

//void Regression::Clear(){
//    data.Clear();
//    labels.clear();
//}

void Regression::Regress(vector<Ctxt> &theta, Ctxt &det){
    //---------------------------------
    //           X^T * y
    //---------------------------------
    vector<vector<Ctxt>> last;
    // cout << endl << "Multiply Vector!" << endl;
    MultiplyVector(last, data, labels, true, false);
    // cout << endl << "totalSums for [last]!" << endl;
    for (unsigned i = 0; i < last.size(); i++) {
        for (unsigned j = 0; j < last[0].size(); j++) {
            totalSums(last[i][j]);
        }
    }
    
    // Print and check the result of data * labels
    // cout << endl << "1. The result of data * labels for enc version!" << endl;
    // for (int i = 0; i < last.size(); i++) {
    //     for (int j = 0; j < last[0].size(); j++) {
    //         vector<long> temp_dec;
    //         meta.data->ea.decrypt(last[i][j], meta.data->secretKey, temp_dec);
    //         cout << temp_dec[0] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;
    
    //---------------------------------
    //           X^T * X
    //---------------------------------
    vector<vector<Ctxt>> res_MultByTrans;
    // cout << endl << "MultByTranspose!" << endl;
    MultByTranspose(res_MultByTrans, data, true);
    // cout << endl << "totalSums for [res_MultByTrans]!" << endl;
    for (unsigned i = 0; i < res_MultByTrans.size(); i++) {
        for (unsigned j = 0; j < res_MultByTrans[0].size(); j++) {
            totalSums(res_MultByTrans[i][j]);
        }
    }
    
    // Check the result of data^T * data
    // cout << endl << "2. The result of data^T * data for enc version!" << endl;
    // for (int i = 0; i < res_MultByTrans.size(); i++) {
    //     for (int j = 0; j < res_MultByTrans[0].size(); j++) {
    //         vector<long> temp_dec;
    //         meta.data->ea.decrypt(res_MultByTrans[i][j], meta.data->secretKey, temp_dec);
    //         cout << temp_dec[0] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;

    //---------------------------------
    //         (X^T * X)^(-1)
    //---------------------------------
    auto inverse_start= chrono::steady_clock::now();

    long numCols = data[0].size();
    if (numCols == 1) {
        det = res_MultByTrans[0][0];
        theta.assign(1, last[0][0]);
        return;
    }

    vector<vector<Ctxt>> adj;
    Invert(det, adj, res_MultByTrans, false);

    auto inverse_end = std::chrono::steady_clock::now();
    auto inverse_diff = inverse_end - inverse_start;
    auto inverse_timeElapsed = chrono::duration <double, milli> (inverse_diff).count()/1000.0;
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Matrix Inverse under Encryption Time = " << inverse_timeElapsed << " s" << endl;
    cout << "------------------------------------------------------------------------" << endl;
  
    // Check the result of data^T * data
    // cout << endl << "3. The result of Invert(data^T * data） for enc version!" << endl;
    // for (int i = 0; i < adj.size(); i++) {
    //     for (int j = 0; j < adj[0].size(); j++) {
    //         vector<long> temp_dec;
    //         meta.data->ea.decrypt(adj[i][j], meta.data->secretKey, temp_dec);
    //         cout << temp_dec[0] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;
    

    //---------------------------------
    //  [(X^T * X)^(-1)] * [X^T * y]
    //---------------------------------
    vector<vector<Ctxt>> final_result;
    // cout << endl << "MultiplyMatrix!" << endl;
    MultiplyMatrix(final_result, adj, last, false, false);

    // cout << endl << "Assign result to theta!" << endl;
    long final_res_cols = final_result.size();
    // theta.resize(final_res_cols, Ctxt(meta.data->publicKey));
    // for (unsigned i = 0; i < final_res_cols; i++)
    // {
    //     theta[i] = final_result[i][0];
    // }
    for (unsigned i = 0; i < final_res_cols; i++)
    {
        theta.push_back(final_result[i][0]);
    }
}

void RegressPT(vector<long> &theta, long &det, Mat<long> &data, vector<long> &labels){
    PtMatrix<long> A;
    for (long i = 0; i < data.NumRows(); i++) {
        vector<long> tempData;
        for (long j = 0; j < data.NumCols(); j++) {
            tempData.push_back(data[i][j]);
        }
        A.AddRow(tempData);
    }
    A.Transpose();

    PtMatrix<long> tmp = A * labels;
    
    // Check the result of A * labels
    cout << endl << "1. The result of A * labels!" << endl;
    for (int i = 0; i < tmp.NumRows(); i++) {
        for (int j = 0; j < tmp.NumCols(); j++) {
            cout << tmp(i,j) << " ";
        }
        cout << endl;
    }
    cout << endl;

    A.MultByTranspose();
    cout << endl << "2. The result of A^T * A!" << endl;
    for (int i = 0; i < A.NumRows(); i++) {
        for (int j = 0; j < A.NumCols(); j++) {
            cout << A(i,j) << " ";
        }
        cout << endl;
    }
    cout << endl;
    
    if (data.NumCols() == 1)
    {
        det = A(0,0);
        theta.assign(1, tmp(0,0));
        return;
    }

    A.Invert(det);
    cout << endl << "3. The result of Invert(A)!" << endl;
    for (int i = 0; i < A.NumRows(); i++) {
        for (int j = 0; j < A.NumCols(); j++) {
            cout << A(i,j) << " ";
        }
        cout << endl;
    }
    cout << endl;
    
    A *= tmp;

    theta.resize(A.NumRows());
    for (unsigned i = 0; i < A.NumRows(); i++)
    {
        theta[i] = A(i,0);
    }
}

void Regression::CheckWithThinBoot(Ctxt& c){
    if (c.bitCapacity() <= 50)
    {
        cout << endl << "Invoke bootstrapping!" << endl;
//        meta.data->publicKey.thinReCrypt(c);
        meta.data->publicKey.reCrypt(c);
    }
}

void Regression::MultiplyMatrix(vector<vector<Ctxt>> &output, vector<vector<Ctxt>> &input1, vector<vector<Ctxt>> &input2, bool isTranspose1, bool isTranspose2){
    long row = isTranspose1 ? input1[0].size() : input1.size();
    long col1 = isTranspose1 ? input1.size() : input1[0].size();
    long col = isTranspose2 ? input2.size() : input2[0].size();
    output = vector<vector<Ctxt>>(row, vector<Ctxt>(col, Ctxt(meta.data->publicKey)));
    
    for (long i = 0; i < row; i++)
    {
        for (long j = 0; j < col; j++)
        {
            output[i][j] = isTranspose1 ? input1[0][i] : input1[i][0];
            output[i][j] *= (isTranspose2 ? input2[j][0] : input2[0][j]);
            CheckWithThinBoot(output[i][j]);
            
            for (long k = 1; k < col1; k++)
            {
                Ctxt tmp = isTranspose1 ? input1[k][i] : input1[i][k];
                tmp *= (isTranspose2 ? input2[j][k] : input2[k][j]);
                CheckWithThinBoot(tmp);
                output[i][j] += tmp;
          }
        }
    }
}

// Functionlaity of Matrix * Vector
void Regression::MultiplyVector(vector<vector<Ctxt>> &output, vector<vector<Ctxt>> &input1, vector<Ctxt> &input2, bool isTranspose1, bool isTranspose2){
    long row1 = isTranspose1 ? input1[0].size() : input1.size();
    long col1 = isTranspose1 ? input1.size() : input1[0].size();
    output = vector<vector<Ctxt>>(row1, vector<Ctxt>(1, Ctxt(meta.data->publicKey)));
    
    for (long i = 0; i < row1; i++)
    {
        output[i][0] = isTranspose1 ? input1[0][i] : input1[i][0];
        output[i][0] *= input2[0];
        CheckWithThinBoot(output[i][0]);
        for (long j = 1; j < col1; j++)
        {
            Ctxt tmp = isTranspose1 ? input1[j][i] : input1[i][j];
            tmp *= input2[j];
            CheckWithThinBoot(tmp);
            output[i][0] += tmp;
        }
    }
}

void Regression::MultByTranspose(vector<vector<Ctxt>> &output, vector<vector<Ctxt>> &input, bool isTranspose){
    long row = isTranspose ? input[0].size() : input.size();
    long col = isTranspose ? input.size() : input[0].size();
    output = vector<vector<Ctxt>>(row, vector<Ctxt>(row, Ctxt(meta.data->publicKey)));
    
    for (unsigned i = 0; i < row; i++)
    {
        for (long j = i; j < row; j++)
        {
            output[i][j] = isTranspose ? input[0][i] : input[i][0];
            output[i][j] *= isTranspose ? input[0][j] : input[j][0];
            CheckWithThinBoot(output[i][j]);
            
            for (long k = 1; k < col; k++)
            {
                Ctxt tmp = isTranspose ? input[k][i] : input[i][k];
                tmp *= isTranspose ? input[k][j] : input[j][k];
                CheckWithThinBoot(tmp);
                output[i][j] += tmp;
            }
          
            if (i != j)
            {
                output[j][i] = output[i][j];
            }
        }
    }
}

void Regression::Invert(Ctxt &det, vector<vector<Ctxt>> &adj, vector<vector<Ctxt>> &input, bool isTranspose){
    long dim = isTranspose ? input[0].size() : input.size();
    adj = vector<vector<Ctxt>>(dim, vector<Ctxt>(dim, Ctxt(meta.data->publicKey)));
    
    vector<bool> usedRows(dim), usedCols(dim);
    for (long i = 0; i < dim; i++)
    {
        for (long j = 0; j < dim; j++)
        {
          usedRows[i] = usedCols[j] = true;
          Determinant(adj[j][i], usedRows, usedCols, dim-1, input, isTranspose);
          usedRows[i] = usedCols[j] = false;
          if ((i+j) % 2 == 1)
          {
              adj[j][i].multByConstant(long(-1));
          }
        }
    }
    
    // Since we already computed the adjugate matrix, we can
    // reuse the results to compute the determinant
    det = input[0][0];
    det *= adj[0][0];
    CheckWithThinBoot(det);
    for (unsigned i = 1; i < dim; i++)
    {
        Ctxt tmp = isTranspose ? input[i][0] : input[0][i];
        tmp *= adj[i][0];
        CheckWithThinBoot(tmp);
        det += tmp;
    }
}

void Regression::Determinant(Ctxt& det, vector<bool> &usedRows, vector<bool> &usedCols, long dim, vector<vector<Ctxt>> &input, bool isTranspose){
    long matDim = isTranspose ? input[0].size() : input.size();
    long row = 0;
    while (usedRows[row]) row++;
    
    bool negative = false;
    bool first = true;
    
    for (unsigned col = 0; col < matDim; col++)
    {
        if (usedCols[col]) continue;
     
        if (dim == 1)
        {
            det = isTranspose ? input[col][row] : input[row][col];
            return;
        }
     
        Ctxt tmp = isTranspose ? input[col][row] : input[row][col];
        if (negative) tmp.multByConstant(long(-1));
        negative = !negative;
        
        usedRows[row] = usedCols[col] = true;
        Ctxt tmp2(meta.data->publicKey);
        Determinant(tmp2, usedRows, usedCols, dim-1, input, isTranspose);
        usedRows[row] = usedCols[col] = false;
        
        tmp *= tmp2;
        CheckWithThinBoot(tmp);
        
        if (first)
        {
            det = tmp;
            first = false;
        } else
        {
            det += tmp;
        }
    }
}
