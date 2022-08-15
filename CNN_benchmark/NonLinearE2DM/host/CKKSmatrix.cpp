//
//  CKKSmatrix.cpp
//  ePPDSC
//
//  Created by Qifan Wang on 30/12/21.
//

#include <cmath>
#include <math.h>  // pow
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <string>
#include <vector>
#include <random>
#include <iomanip>

#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZX.h>
#include <NTL/mat_RR.h>
#include <NTL/vec_RR.h>

#include "CKKSmatrix.h"

void readckksMatpar(ckksMatpar& ckksmatpar, long nrows, long ncols, long subdim, long nbatching){
    ckksmatpar.nrows = nrows;
    ckksmatpar.ncols = ncols;
        
    long dim = nrows;
    if(dim < ncols) {
        dim = ncols;
    }
    
    ckksmatpar.dim = (1<< (long)ceil(log2(dim)));  //! power of two
    ckksmatpar.dim1 = ckksmatpar.dim - 1;
    ckksmatpar.logdim = (long) log2(ckksmatpar.dim);
    ckksmatpar.sqrdim = (long) ceil(sqrt(ckksmatpar.dim));
    ckksmatpar.subdim = subdim;  //! used for non-squared matrix multiplication
    ckksmatpar.nbatching = (1<< (long)ceil(log2(nbatching)));
    
    if(nbatching == 1){
        ckksmatpar.neededslots = ckksmatpar.dim * ckksmatpar.dim;  // just encode a single matrix
    }
    else{
        ckksmatpar.neededslots = (ckksmatpar.dim * ckksmatpar.dim) * ckksmatpar.nbatching;
    }
}

void printRvector(vec_RR& vec, long print_size){
    long len;
    
    if(print_size == 0){
        len = vec.length();
    }
    else{
        len = print_size;
    }
    
    cout << "   [" ;
    for(int i = 0; i < len; ++i){
        cout << " " << vec[i] << ((i != len - 1) ? "\t" : "]\n");
    }
}


//!@ Input: RR-matrix
//!@ Function: print the matrix
void printRmatrix(Mat<RR>& mat, const long print_size){
    long rlen, clen;
    
    if(print_size == 0){
        rlen = mat.NumRows();
        clen = mat.NumCols();
    }
    else{
        rlen = print_size;
        clen = print_size;
    }
    
    cout << "From printRmatrix!" << endl;
    for(int i = 0; i< rlen; ++i){
        cout << "   [";
        for(int j = 0; j < clen; ++j){
            cout << mat[i][j] << ((j != clen - 1) ? "\t" : "]\n");
        }
    }
}

//!@ Input: A and B
//!@ Function: return the maximum norm of the difference of two input matrices A and B
RR getError(mat_RR Amat, mat_RR Bmat, long nrows, long ncols){
    RR ret = to_RR("0");
    
    for(long i = 0; i < nrows; ++i){
        for(long j = 0; j < ncols; ++j){
            RR temp = abs(Amat[i][j]-Bmat[i][j]);
            if (ret < temp){
                ret = temp;
            }
            if(temp > 1e-2){
                cout << "From getError!" << endl;
                cout << "(" << i << "," << j  << ") = " << Amat[i][j] << ", " << Bmat[i][j]<< endl;
            }
        }
    }
    return ret;
}

void CKKSmatrix::encryptRmat(Ctxt& ctxt, mat_RR& mat){
    vector<complex<double>> cmsg(ckksmatpar.neededslots);
    
    NTL_EXEC_RANGE(ckksmatpar.nrows, first, last);
    for(long i = first; i < last; i++){
        for(long j = 0; j < ckksmatpar.ncols; j++){
            double dtemp;
            conv(dtemp, mat[i][j]);
            cmsg[i * ckksmatpar.dim + j].real(dtemp);
        }
    }
    NTL_EXEC_RANGE_END;
    
//    ckksmeta.data->ea.encrypt(ctxt, ckksmeta.data->publicKey, cmsg, 1.0);
    PtxtArray pa(ckksmeta.data->context, cmsg);
    pa.encrypt(ctxt);
}

void CKKSmatrix::decryptRmat(mat_RR& mat, Ctxt& ctxt){
    vector<complex<double>> cmsg;
    PtxtArray pa(ckksmeta.data->context);
    pa.decryptComplex(ctxt, ckksmeta.data->secretKey);
    pa.store(cmsg);
    mat.SetDims(ckksmatpar.dim, ckksmatpar.dim);
    
    long k = 0;
    for(long i = 0; i < ckksmatpar.dim; i++){
        for(long j = 0; j < ckksmatpar.dim; j++){
            mat[i][j] = to_RR(cmsg[k].real());
            k++;
        }
    }
}

void CKKSmatrix::encryptParallelRmat(Ctxt& ctxt, mat_RR*& mat, long nbatching){
//    vector<complex<double>> cmsg(ckksmatpar.neededslots);
    vector<complex<double>> cmsg(ckksmeta.data->ea.size());
    double* dtemp = new double[nbatching];
    
    NTL_EXEC_RANGE(nbatching, first, last);
    for(long k = first; k < last; ++k){
        // encode the kth matrix
        for(long i = 0; i < ckksmatpar.nrows; ++i){
            for(long j = 0; j < ckksmatpar.ncols; ++j){
                conv(dtemp[k], mat[k][i][j]);
                cmsg[(i * ckksmatpar.dim + j)* ckksmatpar.nbatching + k].real(dtemp[k]);
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    PtxtArray pa(ckksmeta.data->context, cmsg);
    pa.encrypt(ctxt);
    delete[] dtemp;
}

void CKKSmatrix::decryptParallelRmat(mat_RR*& mat, Ctxt& ctxt){
    vector<complex<double>> cmsg;
    PtxtArray pa(ckksmeta.data->context);
    pa.decryptComplex(ctxt, ckksmeta.data->secretKey);
    pa.store(cmsg);
    mat = new mat_RR[ckksmatpar.nbatching];
    
    for(long k = 0; k < ckksmatpar.nbatching; ++k){
        mat[k].SetDims(ckksmatpar.dim, ckksmatpar.dim);
        long l = 0;
        for(long i = 0; i < ckksmatpar.dim; ++i){
            for(long j = 0; j < ckksmatpar.dim; ++j){
                mat[k][i][j] = to_RR(cmsg[l * ckksmatpar.nbatching + k].real());
                l++;
            }
        }
    }
}

void CKKSmatrix::msgleftRotate(vector<complex<double>>& res, vector<complex<double>> vals, long dim, long nrot){
    long nshift = nrot % ckksmatpar.neededslots;
        
    long k = dim - nshift;
    for(long j = 0; j < k; ++j){
        res[j] = vals[j + nshift];
    }
    for(long j = k; j < dim; ++j){
        res[j] = vals[j - k];
    }
}

void CKKSmatrix::msgrightRotate(vector<complex<double>>& res, vector<complex<double>> vals, long dim, long nrot){
    long nshift = nrot % ckksmatpar.neededslots;
    long k = dim - nshift;
    
    //! vals[0],vals[1],....,vals[d-nrot-1])
    for(long j = 0; j < k; ++j){
        res[nshift + j] = vals[j];
    }
    //! (vals[d-nrot],vals[d-nrot+1],...,vals[d-1])
    for(long j = k; j < dim; ++j){
        res[j - k] = vals[j];
    }
}

void CKKSmatrix::msgleftRotateAndEqual(vector<complex<double>>& vals, long dim, long nrot){
    vector<complex<double>> res(dim);

    long nshift = nrot % ckksmatpar.neededslots;
    long k = dim - nshift;
    
    for(long j = 0; j < k; ++j){
        res[j] = vals[j + nshift];
    }
    for(long j = k; j < dim; ++j){
        res[j] = vals[j - k];
    }
    for(long j = 0; j < dim; ++j){
        vals[j] = res[j];
    }
}

void CKKSmatrix::msgrightRotateAndEqual(vector<complex<double>>& vals, long dim, long nrot){
    vector<complex<double>> res(dim);
        
    long nshift = nrot % ckksmatpar.neededslots;
    long k = dim - nrot;
    
    for(long j = 0; j < k; ++j){
        res[nshift + j] = vals[j];
    }
    for(long j = k; j < dim; ++j){
        res[j - k] = vals[j];
    }
    for(long j = 0; j < dim; ++j){
        vals[j] = res[j];
    }
}

void CKKSmatrix::genTransPoly(vector<EncodedPtxt>& transpoly){
    // define 1-dim's size
    vector<vector<complex<double>>> lvals(ckksmatpar.dim);
    vector<vector<complex<double>>> rvals(ckksmatpar.dim);
    long dsquare = ckksmatpar.neededslots - 1;
    transpoly = vector<EncodedPtxt>(2 * ckksmatpar.dim); //change
    long dimsqrdim = ckksmatpar.sqrdim * ckksmatpar.dim1;
    
    // k = i * sqrt(d) + j < d
    bool btmp;
    if((ckksmatpar.dim % ckksmatpar.sqrdim) == 0){
        btmp = false;   //! all the terms have the same numbers
    }
    else{
        btmp = true;
    }
    //! number of "i"
    long ibound = (long) ceil((double)ckksmatpar.dim/ckksmatpar.sqrdim);
    
    NTL_EXEC_RANGE(ibound, first, last);
    for(long i = first; i < last; ++i){
        long jbound = ckksmatpar.sqrdim;
        if ((btmp) && (i == ibound - 1)){   //! last term
            jbound = (ckksmatpar.dim % ckksmatpar.sqrdim);
        }
        
        for(long j = 0; j < jbound; ++j){
            long k = i * ckksmatpar.sqrdim + j;
            
            lvals[k] = vector<complex<double>>(ckksmatpar.neededslots);
            rvals[k] = vector<complex<double>>(ckksmatpar.neededslots);
            
            for(long l = 0; l < ckksmatpar.dim - k; ++l){
                long dtemp = l * (ckksmatpar.dim + 1) + k;
                lvals[k][dtemp].real(1.0);
                rvals[k][dsquare - dtemp].real(1.0);
            }
            
            msgrightRotateAndEqual(lvals[k], ckksmatpar.neededslots, i * dimsqrdim);
            msgleftRotateAndEqual(rvals[k], ckksmatpar.neededslots, i*dimsqrdim);
            
            PtxtArray pa1(ckksmeta.data->context, lvals[k]);
            pa1.encode(transpoly[k]);
            PtxtArray pa2(ckksmeta.data->context, rvals[k]);
            pa2.encode(transpoly[k + ckksmatpar.dim]);
        }
    }
    NTL_EXEC_RANGE_END;
}

void CKKSmatrix::transpose(Ctxt& res, Ctxt& ctxt, vector<EncodedPtxt>& transpoly){
    long dimsqrdim = ckksmatpar.sqrdim * ckksmatpar.dim1;
        
    bool btmp;
    if((ckksmatpar.dim % ckksmatpar.sqrdim) == 0){
        btmp = false;
    }
    else{
        btmp = true;
    }
    long ibound = (long) ceil((double) ckksmatpar.dim/ckksmatpar.sqrdim); //! number of "i"
    
    vector<Ctxt> Babyctxt1(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey));
    vector<Ctxt> Babyctxt2(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey));
    
    vector<vector<Ctxt>> ltemp(ckksmatpar.sqrdim, vector<Ctxt>(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey)));
    vector<vector<Ctxt>> rtemp(ckksmatpar.sqrdim, vector<Ctxt>(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey)));
    
    //! diagonal
//    res = scheme.multByPoly(ctxt, transpoly[0], ckksmatpar.cBits);
    res = ctxt;
    res.multByConstant(transpoly[0]);
    
    Babyctxt1[0] = ctxt;
    Babyctxt2[0] = ctxt;
    
    //! Babyctxt1[j] = rho(v; j * (d-1))
    //! res[j] = rho(-, (d-1)*srt(d) * i) * rho(v; j * (d-1))
    //ckksmatpar.sqrdim
    
    NTL_EXEC_RANGE(ckksmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);  // 1 <= j1 <= sqrdim - 1

        Babyctxt1[j1] = ctxt;
        rotate(Babyctxt1[j1], (-j1 * (ckksmatpar.dim1)));
        
        Babyctxt2[j1] = ctxt;
        rotate(Babyctxt2[j1], j1 * (ckksmatpar.dim1));
        
        ltemp[0][j1] = Babyctxt1[j1];
        ltemp[0][j1].multByConstant(transpoly[j1]);
        
        rtemp[0][j1] = Babyctxt2[j1];
        rtemp[0][j1].multByConstant(transpoly[j1 + ckksmatpar.dim]);

        ltemp[0][j1] += rtemp[0][j1];
    }
    NTL_EXEC_RANGE_END;

    for(long j = 1; j < ckksmatpar.sqrdim; ++j){
        res += ltemp[0][j];
    }

    //! 1 <= i < d
    //! [k]: rho(p[k];-dsqrt(d)*i1) * rho(v, j*(d-1)) -> rot by dsqrt(d)*i
    //! res[j] = rho( rho(-, (d-1)*srt(d) * i) * rho(v; j * (d-1)) ; i * (d-1)\sqrt(d) )
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long i = first; i < last; ++i){
        long jbound = ckksmatpar.sqrdim;
        if (btmp &&(i == (ibound - 2))){
            jbound = (ckksmatpar.dim % ckksmatpar.sqrdim);
        }
        
        long i1 = i + 1;
        long k = (i1) * ckksmatpar.sqrdim;

        //! 0 <= j < sqrdim
        NTL_EXEC_RANGE(jbound, first, last);
        for(long j = first; j < last; ++j){
            ltemp[i1][j] = Babyctxt1[j];
            ltemp[i1][j].multByConstant(transpoly[k + j]);
            rtemp[i1][j] = Babyctxt2[j];
            rtemp[i1][j].multByConstant(transpoly[k + j + ckksmatpar.dim]);
        }
        NTL_EXEC_RANGE_END;

        for(long j = 1; j < jbound; ++j){
            ltemp[i1][0] += ltemp[i1][j];
            rtemp[i1][0] += rtemp[i1][j];
        }

        rotate(ltemp[i1][0], (-(i1) * dimsqrdim));
        rotate(rtemp[i1][0], (i1) * dimsqrdim);

        ltemp[i1][0] += rtemp[i1][0];
    }
    NTL_EXEC_RANGE_END;
    
    for(long i = 1; i < ibound; ++i){
        res += ltemp[i][0];
    }
}

void CKKSmatrix::genTransPoly_Parallel(vector<EncodedPtxt>& transpoly){
    bool btmp;
    if((ckksmatpar.dim % ckksmatpar.sqrdim) == 0)
    {
        btmp = false;
    }
    else{
        btmp = true;
    }
    long ibound = (long) ceil((double)ckksmatpar.dim/ckksmatpar.sqrdim);
    
    vector<vector<complex<double>>> lvals(ckksmatpar.dim);
    vector<vector<complex<double>>> rvals(ckksmatpar.dim);
    
    transpoly = vector<EncodedPtxt>(2 * ckksmatpar.dim);
    
    long shiftunit = ckksmatpar.sqrdim * ckksmatpar.dim1 * ckksmatpar.nbatching;
    long dsquare = ((ckksmatpar.dim * ckksmatpar.dim) - 1) * ckksmatpar.nbatching;
    
    //! k = (i * sqrt(d) + j)  < d
    //! poly[k] = poly[i * sqrt(d) + j]  -> rho(poly[k]; i * (d-1)*sqrt(d) * l)
    for(long i = 0; i < ibound; i++){
        long jbound = ckksmatpar.sqrdim;
        if ((btmp) && (i == ibound - 1)){
            jbound = (ckksmatpar.dim % ckksmatpar.sqrdim);
        }
        for(long j = 0; j < jbound; j++){
            long k = (i * ckksmatpar.sqrdim + j);
        
            lvals[k] = vector<complex<double>>(ckksmatpar.neededslots);
            rvals[k] = vector<complex<double>>(ckksmatpar.neededslots);
            
            for(long l = 0; l < ckksmatpar.dim - k; l++){
                long dtemp= (l * (ckksmatpar.dim + 1) + k) * ckksmatpar.nbatching;
        
                for(long n = 0; n < ckksmatpar.nbatching; n++){
                    lvals[k][dtemp + n].real(1.0);
                    rvals[k][dsquare - dtemp + n].real(1.0);
                }
            }
        
            msgrightRotateAndEqual(lvals[k], ckksmatpar.neededslots, i * shiftunit);
            msgleftRotateAndEqual(rvals[k], ckksmatpar.neededslots, i * shiftunit);
            
            PtxtArray pa1(ckksmeta.data->context, lvals[k]);
            pa1.encode(transpoly[k]);
            PtxtArray pa2(ckksmeta.data->context, rvals[k]);
            pa2.encode(transpoly[k + ckksmatpar.dim]);
        }
    }
}

void CKKSmatrix::transpose_Parallel(Ctxt& res, Ctxt& ctxt, vector<EncodedPtxt>& transpoly){
    long shiftunit = ckksmatpar.sqrdim * ckksmatpar.dim1 * ckksmatpar.nbatching;  // (d-1) * sqrt(d) * l
    long unit  = ckksmatpar.dim1 * ckksmatpar.nbatching;
    
    bool btmp;
    if((ckksmatpar.dim % ckksmatpar.sqrdim) == 0){ btmp = false; }
    else{ btmp = true; }
    long ibound = (long) ceil((double) ckksmatpar.dim/ckksmatpar.sqrdim);
    
    vector<Ctxt> Babyctxt1(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey));
    vector<Ctxt> Babyctxt2(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey));
    
    vector<vector<Ctxt>> ltemp(ckksmatpar.sqrdim, vector<Ctxt>(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey)));
    vector<vector<Ctxt>> rtemp(ckksmatpar.sqrdim, vector<Ctxt>(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey)));
    
    //-------------------------
    // i = 0, sqrt(d) polynomials
    //-------------------------
    res = ctxt;
    res.multByConstant(transpoly[0]);

    Babyctxt1[0] = ctxt;
    Babyctxt2[0] = ctxt;
    
    NTL_EXEC_RANGE(ckksmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);
        
        Babyctxt1[j1] = ctxt;
        rotate(Babyctxt1[j1], (-j1 * unit));
        Babyctxt2[j1] = ctxt;
        rotate(Babyctxt2[j1], j1 * unit);
        
        ltemp[0][j1] = Babyctxt1[j1];
        ltemp[0][j1].multByConstant(transpoly[j1]);
        
        rtemp[0][j1] = Babyctxt2[j1];
        rtemp[0][j1].multByConstant(transpoly[j1 + ckksmatpar.dim]);

        ltemp[0][j1] += rtemp[0][j1];
    }
    NTL_EXEC_RANGE_END;
    
    for(long j = 1; j < ckksmatpar.sqrdim; ++j){
        res += ltemp[0][j];
    }
    
    //! 1 <= i < d
    //! [k]: rho(p[k];-dsqrt(d)*i1) * rho(v, j*(d-1)) -> rot by dsqrt(d)*i
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long i = first; i < last; ++i){
        long jbound = ckksmatpar.sqrdim;
        if ((btmp) && (i == ibound - 2)){
            jbound = (ckksmatpar.dim % ckksmatpar.sqrdim);
        }
        
        long i1 = i + 1;
        long k = (i1) * ckksmatpar.sqrdim;
        
        NTL_EXEC_RANGE(jbound, first, last);
        for(long j = first; j < last; ++j){
            ltemp[i1][j] = Babyctxt1[j];
            ltemp[i1][j].multByConstant(transpoly[k + j]);
            rtemp[i1][j] = Babyctxt2[j];
            rtemp[i1][j].multByConstant(transpoly[k + j + ckksmatpar.dim]);
        }
        NTL_EXEC_RANGE_END;
        
        for(long j = 1; j < jbound; ++j){
            ltemp[i1][0] += ltemp[i1][j];
            rtemp[i1][0] += rtemp[i1][j];
        }
        
        rotate(ltemp[i1][0], (-(i1) * shiftunit));
        rotate(rtemp[i1][0], (i1) * shiftunit);
        
        ltemp[i1][0] += rtemp[i1][0];
    }
    NTL_EXEC_RANGE_END;
    
    for(long i = 1; i < ibound; ++i){
        res += ltemp[i][0];
    }
}

void CKKSmatrix::genShiftPoly(vector<EncodedPtxt>& shiftpoly, long num){
    long length = ckksmatpar.dim1;
    if(num != 0){
        length = num;
    }
    
    vector<vector<complex<double>>> vals(length, vector<complex<double>>(ckksmatpar.neededslots));
    
    shiftpoly = vector<EncodedPtxt>(length);
    
    // i: shifted by (i+1)
    NTL_EXEC_RANGE(length, first, last);
    for(long i = first; i < last; ++i){
        for(long j = 0; j < ckksmatpar.dim; ++j){
            for(long k = 0; k < i + 1; ++k){
                vals[i][j * ckksmatpar.dim + k].real(1.0);
            }
        }
        PtxtArray pa(ckksmeta.data->context, vals[i]);
        pa.encode(shiftpoly[i]);
    }
    NTL_EXEC_RANGE_END;
}

void CKKSmatrix::shiftBycols(Ctxt& res, Ctxt& ctxt, long k, vector<EncodedPtxt>& shiftpoly){
    Ctxt ctemp(ctxt);
    ctemp.multByConstant(shiftpoly[k-1]);
    
    res = ctxt;
    res -= ctemp;
    rotate(ctemp, (ckksmatpar.dim - k));
    rotate(res, -k);
    res += ctemp;
}

void CKKSmatrix::genShiftPoly_Parallel(vector<EncodedPtxt>& shiftpoly){
    vector<vector<complex<double>>> vals(ckksmatpar.dim1, vector<complex<double>>(ckksmatpar.neededslots));
    shiftpoly = vector<EncodedPtxt>(ckksmatpar.dim1);
    
    NTL_EXEC_RANGE(ckksmatpar.dim1, first, last);
    for(long i = first; i < last; ++i){
        for(long j = 0; j < ckksmatpar.dim; ++j){
            for(long k = 0; k < i + 1; ++k){
                long dtemp = (j * ckksmatpar.dim + k) * ckksmatpar.nbatching;
                for(long n = 0; n < ckksmatpar.nbatching; ++n){
                    vals[i][dtemp + n] = 1;
                }
            }
        }
        PtxtArray pa(ckksmeta.data->context, vals[i]);
        pa.encode(shiftpoly[i]);
    }
    NTL_EXEC_RANGE_END;
}

void CKKSmatrix::shiftBycols_Parallel(Ctxt& res, Ctxt& ctxt, long k, vector<EncodedPtxt>& shiftpoly){
    Ctxt ctemp (ctxt);
    ctemp.multByConstant(shiftpoly[k-1]);
    res = ctxt;
    res -= ctemp;
    rotate(ctemp, (ckksmatpar.dim - k) * ckksmatpar.nbatching);
    rotate(res, -(k * ckksmatpar.nbatching));
    res += ctemp;
}

void CKKSmatrix::genMultPoly(vector<vector<EncodedPtxt>>& Initpoly){
    bool btmp;
    if((ckksmatpar.dim % ckksmatpar.sqrdim) == 0){
        btmp = false;
    }
    else{
        btmp = true;
    }
    long ibound = (long) ceil((double) ckksmatpar.dim/ckksmatpar.sqrdim);
    
    Initpoly = vector<vector<EncodedPtxt>>(3, vector<EncodedPtxt>(ckksmatpar.dim));
    
    vector<vector<complex<double>>> fvals1(ckksmatpar.dim);
    vector<vector<complex<double>>> fvals2(ckksmatpar.dim);
    vector<vector<complex<double>>> bvals(ckksmatpar.dim);

    NTL_EXEC_RANGE(ibound, first, last);
    for(long i = first; i < last; ++i){
        long jbound = ckksmatpar.sqrdim;
        if ((btmp)&&(i == ibound - 1)){
            jbound = (ckksmatpar.dim % ckksmatpar.sqrdim);
        }
        
        for(long j = 0; j < jbound; ++j){
            long k = i * ckksmatpar.sqrdim + j;
            fvals1[k] = vector<complex<double>>(ckksmatpar.neededslots);
            fvals2[k] = vector<complex<double>>(ckksmatpar.neededslots);
            bvals[k] = vector<complex<double>>(ckksmatpar.neededslots);
            
            for(long l = 0; l < ckksmatpar.dim - k; ++l){
                fvals1[k][k * ckksmatpar.dim + l].real(1.0);
            }
            msgleftRotate(fvals2[k], fvals1[k], ckksmatpar.neededslots, k * (2 * ckksmatpar.dim - 1));
            
            msgrightRotateAndEqual(fvals1[k], ckksmatpar.neededslots, i * ckksmatpar.sqrdim);
            
            PtxtArray pa1(ckksmeta.data->context, fvals1[k]);
            pa1.encode(Initpoly[0][k]);
        
            msgleftRotateAndEqual(fvals2[k], ckksmatpar.neededslots, i * ckksmatpar.sqrdim);
            
            PtxtArray pa2(ckksmeta.data->context, fvals2[k]);
            pa2.encode(Initpoly[1][k]);
            
            for(long l = 0; l < ckksmatpar.dim; ++l){
                bvals[k][l * ckksmatpar.dim + k].real(1.0);
            }
            
            msgrightRotateAndEqual(bvals[k], ckksmatpar.neededslots, i * ckksmatpar.sqrdim * ckksmatpar.dim);
            
            PtxtArray pa3(ckksmeta.data->context, bvals[k]);
            pa3.encode(Initpoly[2][k]);
        }
    }
    NTL_EXEC_RANGE_END;
}

void CKKSmatrix::genMultPoly_Parallel(vector<vector<EncodedPtxt>>& Initpoly){
    long btmp;
    if((ckksmatpar.dim % ckksmatpar.sqrdim) == 0){ btmp = false; }
    else{ btmp = true; }
    long ibound = (long) ceil((double) ckksmatpar.dim/ckksmatpar.sqrdim);

    
    Initpoly = vector<vector<EncodedPtxt>>(3, vector<EncodedPtxt>(ckksmatpar.dim));
    
    vector<vector<complex<double>>> fvals1(ckksmatpar.dim);
    vector<vector<complex<double>>> fvals2(ckksmatpar.dim);
    vector<vector<complex<double>>> bvals(ckksmatpar.dim);
    
    
    NTL_EXEC_RANGE(ibound, first, last);
    for(int i = first; i < last; ++i){
        long jbound = ckksmatpar.sqrdim;
        if ((btmp) && (i == (ibound - 1))){
            jbound = (ckksmatpar.dim % ckksmatpar.sqrdim);
        }
        
        for(long j = 0; j < jbound; ++j){
            long k = i * ckksmatpar.sqrdim + j;
            fvals1[k] = vector<complex<double>>(ckksmatpar.neededslots);
            fvals2[k] = vector<complex<double>>(ckksmatpar.neededslots);
            bvals[k] = vector<complex<double>>(ckksmatpar.neededslots);
            
            long start = k * ckksmatpar.dim * ckksmatpar.nbatching;
            long end = ( k * ckksmatpar.dim + ckksmatpar.dim - k) * ckksmatpar.nbatching;
            for(long l = start; l < end; ++l){
                fvals1[k][l].real(1.0);
            }
            
            msgleftRotate(fvals2[k], fvals1[k], ckksmatpar.neededslots, k * (2 * ckksmatpar.dim - 1) * ckksmatpar.nbatching);
            
            msgrightRotateAndEqual(fvals1[k], ckksmatpar.neededslots, i * ckksmatpar.sqrdim * ckksmatpar.nbatching);
            
            PtxtArray pa1(ckksmeta.data->context, fvals1[k]);
            pa1.encode(Initpoly[0][k]);
            
            msgleftRotateAndEqual(fvals2[k], ckksmatpar.neededslots, i * ckksmatpar.sqrdim * ckksmatpar.nbatching);
            
            PtxtArray pa2(ckksmeta.data->context, fvals2[k]);
            pa2.encode(Initpoly[1][k]);
            
            for(long l = 0; l < ckksmatpar.dim; ++l){
                long dtemp = (l * ckksmatpar.dim + k) * ckksmatpar.nbatching;
                for(long n = 0; n < ckksmatpar.nbatching; ++n){
                    bvals[k][dtemp + n].real(1.0);
                }
            }
            
            msgrightRotateAndEqual(bvals[k], ckksmatpar.neededslots, i * ckksmatpar.sqrdim * ckksmatpar.dim * ckksmatpar.nbatching);
            
            PtxtArray pa3(ckksmeta.data->context, bvals[k]);
            pa3.encode(Initpoly[2][k]);
        }
    }
    NTL_EXEC_RANGE_END;
}

void CKKSmatrix::genMultBPoly(vector<EncodedPtxt>& Initpoly){
    bool btmp;
    if((ckksmatpar.dim % ckksmatpar.sqrdim) == 0){ btmp = false; }
    else{ btmp = true; }
    long ibound = (long) ceil((double) ckksmatpar.dim/ckksmatpar.sqrdim);
    
    Initpoly = vector<EncodedPtxt>(ckksmatpar.dim);
    
    vector<vector<complex<double>>> bvals(ckksmatpar.dim);
    
    NTL_EXEC_RANGE(ibound, first, last);
    for(long i = first; i < last; ++i){
        long jbound = ckksmatpar.sqrdim;
        if ((btmp) && (i == ibound - 1)){
            jbound = (ckksmatpar.dim % ckksmatpar.sqrdim);
        }
        for(long j = 0; j < jbound; ++j){
            long k = i * ckksmatpar.sqrdim + j;
            bvals[k] = vector<complex<double>>(ckksmatpar.neededslots);
            
            for(long l = 0; l < ckksmatpar.dim; ++l){
                bvals[k][l * ckksmatpar.dim + k].real(1.0);
            }
            
            msgrightRotateAndEqual(bvals[k], ckksmatpar.neededslots, i * ckksmatpar.sqrdim * ckksmatpar.dim);
            PtxtArray pa(ckksmeta.data->context, bvals[k]);
            pa.encode(Initpoly[k]);
        }
    }
    NTL_EXEC_RANGE_END;
}

void CKKSmatrix::genInitCtxt(Ctxt& resA, Ctxt& resB, Ctxt& Actxt, Ctxt& Bctxt, vector<vector<EncodedPtxt>>& Initpoly){
    bool btmp;
    if((ckksmatpar.dim % ckksmatpar.sqrdim) == 0){
        btmp = false;
    }
    else{
        btmp = true;
    }
    
    long ibound = (long) ceil((double) ckksmatpar.dim/ckksmatpar.sqrdim);
    
    vector<vector<Ctxt>> Actemp1(ckksmatpar.sqrdim, vector<Ctxt>(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey)));
    vector<vector<Ctxt>> Actemp2(ckksmatpar.sqrdim, vector<Ctxt>(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey)));
    vector<vector<Ctxt>> Bctemp(ckksmatpar.sqrdim, vector<Ctxt>(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey)));
    
    vector<Ctxt> BaByctxt1(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    vector<Ctxt> BaByctxt2(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    vector<Ctxt> BaByctxtB(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    
    BaByctxt1[0] = Actxt;
    BaByctxt2[0] = Actxt;
    BaByctxtB[0] = Bctxt;
    
    //! i = 0:   Actxts[0] = v[0] + p1 * v[1] + ... + p[sqr(d)-1] *  v[sqr(d)-1]
    resA = BaByctxt1[0];
    resA.multByConstant(Initpoly[0][0]);
    resB = BaByctxtB[0];
    resB.multByConstant(Initpoly[2][0]);
    
    NTL_EXEC_RANGE(ckksmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);
        BaByctxt1[j1] = Actxt;
        rotate(BaByctxt1[j1], -j1);
        
        BaByctxt2[j1] = Actxt;
        rotate(BaByctxt2[j1], j1);
        
        BaByctxtB[j1] = Bctxt;
        rotate(BaByctxtB[j1], -((j1) * ckksmatpar.dim));
        
        Actemp1[0][j1] = BaByctxt1[j1];
        Actemp1[0][j1].multByConstant(Initpoly[0][j1]);
        Actemp2[0][j1] = BaByctxt2[j1];
        Actemp2[0][j1].multByConstant(Initpoly[1][j1]);
        Bctemp[0][j1] = BaByctxtB[j1];
        Bctemp[0][j1].multByConstant(Initpoly[2][j1]);
        
        Actemp1[0][j1] += Actemp2[0][j1];
    }
    NTL_EXEC_RANGE_END;

    for(long j = 1; j < ckksmatpar.sqrdim; ++j){
        resA += Actemp1[0][j];
        resB += Bctemp[0][j];
    }
    
    vector<Ctxt> Actxts1(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    vector<Ctxt> Actxts2(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    vector<Ctxt> Bctxts(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    
    NTL_EXEC_RANGE(ckksmatpar.dim - ckksmatpar.sqrdim, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k + ckksmatpar.sqrdim);
        long i = (long)(k1 / ckksmatpar.sqrdim);
        long j = (long)(k1 % ckksmatpar.sqrdim);
        
        Actxts1[k1] = BaByctxt1[j];
        Actxts1[k1].multByConstant(Initpoly[0][k1]);
        Actxts2[k1] = BaByctxt2[j];
        Actxts2[k1].multByConstant(Initpoly[1][k1]);
        Bctxts[k1] = BaByctxtB[j];
        Bctxts[k1].multByConstant(Initpoly[2][k1]);
    }
    NTL_EXEC_RANGE_END;
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k + 1) * ckksmatpar.sqrdim;
        long jbound = ckksmatpar.sqrdim;
        if((btmp) && (k == ibound - 2)){
            jbound = (ckksmatpar.dim % ckksmatpar.sqrdim);
        }
        for(long j = 1; j < jbound; ++j){
            Actxts1[k1] += Actxts1[k1+j];
            Actxts2[k1] += Actxts2[k1+j];
            Bctxts[k1] += Bctxts[k1+j];
        }
        
        long k2 = (k1 - (k1 % ckksmatpar.sqrdim));
        rotate(Actxts1[k1], -k2);
        rotate(Actxts2[k1], k2);
        rotate(Bctxts[k1], -(k2 * ckksmatpar.dim));
        
        Actxts1[k1] += Actxts2[k1];
    }
    NTL_EXEC_RANGE_END;
    
    for(long k = 1; k < ibound; ++k){
        long k1 = k * ckksmatpar.sqrdim;
        resA += Actxts1[k1];
        resB += Bctxts[k1];
    }
}

void CKKSmatrix::genInitCtxt_Parallel(Ctxt& resA, Ctxt& resB, Ctxt& Actxt, Ctxt& Bctxt, vector<vector<EncodedPtxt>>& Initpoly){
    bool btmp;
    if((ckksmatpar.dim % ckksmatpar.sqrdim) == 0){ btmp = false; }
    else{ btmp = true; }
    long ibound = (long) ceil((double)ckksmatpar.dim/ckksmatpar.sqrdim);

    
    vector<vector<Ctxt>> Actemp1(ckksmatpar.sqrdim, vector<Ctxt>(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey)));
    vector<vector<Ctxt>> Actemp2(ckksmatpar.sqrdim, vector<Ctxt>(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey)));
    vector<vector<Ctxt>> Bctemp(ckksmatpar.sqrdim, vector<Ctxt>(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey)));
    
    vector<Ctxt> BaByctxt1(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    vector<Ctxt> BaByctxt2(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    vector<Ctxt> BaByctxtB(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    
    BaByctxt1[0] = Actxt;
    BaByctxt2[0] = Actxt;
    BaByctxtB[0] = Bctxt;
    
    resA = BaByctxt1[0];
    resA.multByConstant(Initpoly[0][0]);
    resB = BaByctxtB[0];
    resB.multByConstant(Initpoly[2][0]);
    
    NTL_EXEC_RANGE(ckksmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);
        BaByctxt1[j1] = Actxt;
        rotate(BaByctxt1[j1], -(j1 * ckksmatpar.nbatching));
        BaByctxt2[j1] = Actxt;
        rotate(BaByctxt2[j1], j1 * ckksmatpar.nbatching);
        BaByctxtB[j1] = Bctxt;
        rotate(BaByctxtB[j1], -((j1) * ckksmatpar.dim * ckksmatpar.nbatching));
        
        Actemp1[0][j1] = BaByctxt1[j1];
        Actemp1[0][j1].multByConstant(Initpoly[0][j1]);
        Actemp2[0][j1] = BaByctxt2[j1];
        Actemp2[0][j1].multByConstant(Initpoly[1][j1]);
        Bctemp[0][j1] = BaByctxtB[j1];
        Bctemp[0][j1].multByConstant(Initpoly[2][j1]);
        
        Actemp1[0][j1] += Actemp2[0][j1];
    }
    NTL_EXEC_RANGE_END;
    
    for(long j = 1; j < ckksmatpar.sqrdim; ++j){
        resA += Actemp1[0][j];
        resB += Bctemp[0][j];
    }
 
    //---------------------------
    vector<Ctxt> Actxts1(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    vector<Ctxt> Actxts2(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    vector<Ctxt> Bctxts(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    
    NTL_EXEC_RANGE(ckksmatpar.dim - ibound, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k + ckksmatpar.sqrdim);
        long i = (long)(k1 / ckksmatpar.sqrdim);
        long j = (long)(k1 % ckksmatpar.sqrdim);
        
        Actxts1[k1] = BaByctxt1[j];
        Actxts1[k1].multByConstant(Initpoly[0][k1]);
        Actxts2[k1] = BaByctxt2[j];
        Actxts2[k1].multByConstant(Initpoly[1][k1]);
        Bctxts[k1] = BaByctxtB[j];
        Bctxts[k1].multByConstant(Initpoly[2][k1]);
    }
    NTL_EXEC_RANGE_END;
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k + 1) * ckksmatpar.sqrdim;
        long jbound = ckksmatpar.sqrdim;
        if((btmp) && (k == ibound - 2)){
            jbound = (ckksmatpar.dim % ckksmatpar.sqrdim);
        }
        for(long j = 1; j < jbound; ++j){
            Actxts1[k1] += Actxts1[k1+j];
            Actxts2[k1] += Actxts2[k1+j];
            Bctxts[k1] += Bctxts[k1+j];
        }
        
        long k2 = (k1 - (k1 % ckksmatpar.sqrdim)) * ckksmatpar.nbatching;
        rotate(Actxts1[k1], -k2);
        rotate(Actxts2[k1], k2);
        rotate(Bctxts[k1], -(k2 * ckksmatpar.dim));
        
        Actxts1[k1] += Actxts2[k1];
    }
    NTL_EXEC_RANGE_END;
    
    
    for(long k = 1; k < ibound; ++k){
        long k1 = k * ckksmatpar.sqrdim;
        resA += Actxts1[k1];
        resB += Bctxts[k1];
    }
}

void CKKSmatrix::genInitActxt(vector<Ctxt>& Actxts, Mat<RR>& mat){
    Mat<RR>* Amat = new Mat<RR>[ckksmatpar.dim];
    Actxts= vector<Ctxt>(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    vector<vector<complex<double>>> cmsg(ckksmatpar.dim, vector<complex<double>>(ckksmatpar.neededslots));

    NTL_EXEC_RANGE(ckksmatpar.dim, first, last);
    for(long k = first; k < last; ++k){
        Amat[k].SetDims(ckksmatpar.dim, ckksmatpar.dim);
        long dimk = ckksmatpar.dim - k;
        
        for(long i = 0; i < dimk; ++i){
            long nshift = k + i;
            long nshift2 = ckksmatpar.dim - nshift;
            
            for(long j = 0; j < nshift2; ++j){
                Amat[k][i][j] = mat[i][j + nshift];
            }
            
            for(long j = nshift2; j < ckksmatpar.dim; ++j){
                Amat[k][i][j] = mat[i][j - nshift2];
            }
        }
        
        //! i = d - k : Amat[k][i] <- mat[i]
        if(k!=0){
            for(long j = 0; j < ckksmatpar.dim; ++j){
                Amat[k][dimk][j] = mat[dimk][j];
            }
        }
        
        //! d - k + 1 <= i < d: shift by (k+i-d)-positions from mat[i]
        for(long i = dimk + 1; i < ckksmatpar.dim; ++i){
            long nshift =  i - dimk;
            long nshift2 = ckksmatpar.dim - nshift;
            
            for(long j = 0; j < nshift2; ++j){
                Amat[k][i][j] = mat[i][j + nshift];
            }
            for(long j = nshift2; j < ckksmatpar.dim; ++j){
                Amat[k][i][j] = mat[i][j - nshift2];
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! encryption of d Rmat
    NTL_EXEC_RANGE(ckksmatpar.dim, first, last);
    for(long k = first; k < last; ++k){
        for(long i = 0; i < ckksmatpar.nrows; ++i){
            for(long j = 0; j < ckksmatpar.ncols; ++j){
                double dtemp;
                conv(dtemp, Amat[k][i][j]);
                cmsg[k][i * ckksmatpar.dim + j].real(dtemp);
            }
        }
        PtxtArray pa(ckksmeta.data->context, cmsg[k]);
        pa.encrypt(Actxts[k]);
    }
    NTL_EXEC_RANGE_END;
}

void CKKSmatrix::genInitBctxt(Ctxt& resB, Ctxt& Bctxt, vector<EncodedPtxt>& Initpoly){
    bool btmp;
    if((ckksmatpar.dim % ckksmatpar.sqrdim) == 0){
        btmp = false;
    }
    else{
        btmp = true;
    }
    long ibound = (long) ceil((double)ckksmatpar.dim/ckksmatpar.sqrdim);
    
    
    vector<vector<Ctxt>> Bctemp(ckksmatpar.sqrdim, vector<Ctxt>(ckksmatpar.sqrdim, Ctxt(ckksmeta.data->publicKey)));
    
    //! 0. Store some ciphertexts (0,1,...,d-1), ( ,d+1,...2d-1)
    vector<Ctxt> BaByctxtB(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    
    BaByctxtB[0] = Bctxt;
    
    //! i = 0:   Actxts[0] = v[0] + p1 * v[1] + ... + p[sqr(d)-1] *  v[sqr(d)-1]
    resB = BaByctxtB[0];
    resB.multByConstant(Initpoly[0]);

    NTL_EXEC_RANGE(ckksmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);
        BaByctxtB[j1] = Bctxt;
        rotate(BaByctxtB[j1], -((j1) * ckksmatpar.dim));
        Bctemp[0][j1] = BaByctxtB[j1];
        Bctemp[0][j1].multByConstant(Initpoly[j1]);
    }
    NTL_EXEC_RANGE_END;
    
    for(long j = 1; j < ckksmatpar.sqrdim; ++j){
        resB += Bctemp[0][j];
    }
    
    vector<Ctxt> Bctxts(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    
    NTL_EXEC_RANGE(ckksmatpar.dim - ckksmatpar.sqrdim, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k + ckksmatpar.sqrdim);
        long i = (long)(k1 / ckksmatpar.sqrdim);
        long j = (long)(k1 % ckksmatpar.sqrdim);
        
        Bctxts[k1] = BaByctxtB[j];
        Bctxts[k1].multByConstant(Initpoly[k1]);
    }
    NTL_EXEC_RANGE_END;
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k+1) * ckksmatpar.sqrdim;
        long jbound = ckksmatpar.sqrdim;
        if((btmp) && (k == ibound - 2)){
            jbound = (ckksmatpar.dim % ckksmatpar.sqrdim);
        }
        
        for(long j = 1; j < jbound; ++j){
            Bctxts[k1] += Bctxts[k1+j];
        }
        
        long k2 = (k1 - (k1 % ckksmatpar.sqrdim));
        rotate(Bctxts[k1], -(k2 * ckksmatpar.dim));
    }
    NTL_EXEC_RANGE_END;
    
    for(long k = 1; k < ibound; ++k){
        long k1 = k * ckksmatpar.sqrdim;
        resB += Bctxts[k1];
    }
}

void CKKSmatrix::genInitRecActxt(vector<Ctxt>& Actxts, Mat<RR>& mat){
    Mat<RR> replicate_mat;
    replicate_mat.SetDims(ckksmatpar.dim, ckksmatpar.dim);
    
    long index_rows = ckksmatpar.dim/ckksmatpar.subdim;
    
    NTL_EXEC_RANGE(index_rows, first, last);
    for(long k = first; k < last; ++k){
        for(long i = 0; i < ckksmatpar.subdim; ++i){
            for(long j = 0; j < ckksmatpar.dim; ++j){
                replicate_mat[k * ckksmatpar.subdim + i][j] = mat[i][j];
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! generate the (linear transformed) matrices
    Mat<RR>* Amat = new Mat<RR>[ckksmatpar.subdim];
    Actxts = vector<Ctxt>(ckksmatpar.subdim, Ctxt(ckksmeta.data->publicKey));
    
    NTL_EXEC_RANGE(ckksmatpar.subdim, first, last);
    for(long k = first; k < last; ++k){
        Amat[k].SetDims(ckksmatpar.dim, ckksmatpar.dim);
        long dimk= ckksmatpar.dim - k;
        
        //! 0 <= i < d - k: shift by (k+i)-positions from mat[i]
        for(long i = 0; i < dimk; ++i){
            long nshift = k + i;
            long nshift2 = ckksmatpar.dim - nshift;
            
            for(long j = 0; j < nshift2; ++j){
                Amat[k][i][j] = replicate_mat[i][j + nshift];
            }
            
            for(long j = nshift2; j < ckksmatpar.dim; ++j){
                Amat[k][i][j] = replicate_mat[i][j - nshift2];
            }
        }
        
        //! i = d - k : Amat[k][i] <- mat[i]
        if(k!=0){
            for(long j = 0; j < ckksmatpar.dim; ++j){
                Amat[k][dimk][j] = replicate_mat[dimk][j];
            }
        }
        
        //! d - k + 1 <= i < d: shift by (k+i-d)-positions from mat[i]
        for(long i = dimk + 1; i < ckksmatpar.dim; ++i){
            long nshift =  i - dimk;
            long nshift2 = ckksmatpar.dim - nshift;
            
            for(long j = 0; j < nshift2; ++j){
                Amat[k][i][j] = replicate_mat[i][j + nshift];
            }
            for(long j = nshift2; j < ckksmatpar.dim; ++j){
                Amat[k][i][j] = replicate_mat[i][j - nshift2];
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! encryption of d Rmat
    vector<vector<complex<double>>> cmsg(ckksmatpar.subdim, vector<complex<double>>(ckksmatpar.neededslots));
    NTL_EXEC_RANGE(ckksmatpar.subdim, first, last);
    for(long k = first; k < last; ++k){
        for(int i = 0; i < ckksmatpar.nrows; ++i){
            for(long j = 0; j < ckksmatpar.ncols; ++j){
                double dtemp;
                conv(dtemp, Amat[k][i][j]);
                cmsg[k][i * ckksmatpar.dim + j].real(dtemp);
            }
        }
        PtxtArray pa(ckksmeta.data->context, cmsg[k]);
        pa.encrypt(Actxts[k]);
    }
    NTL_EXEC_RANGE_END;
}

void CKKSmatrix::HEmatmul_Hadamard(Ctxt& res, vector<Ctxt> Actxts, vector<Ctxt> Bctxts, long num){
    long num1 = num - 1;
    // we remove ModDown operations in HEAAN, including:
    // scheme.modDownToAndEqual
    // Fixme: check Actxts[0].modulus > Bctxts[0].modulus? In multiplyBy, findBaseSet is invoked.
    res = Actxts[0];
    res.multiplyBy(Bctxts[0]);
    
    NTL_EXEC_RANGE(num1, first, last);
    for(long i = first; i < last; ++i){
        long i1 = (i + 1);
         Actxts[i1].multiplyBy(Bctxts[i1]);
        // NB: reduce the reLinearize usage in multiplyBy()
        // Actxts[i1].multLowLvl(Bctxts[i1]);
    }
    NTL_EXEC_RANGE_END;
    
    for(int i = 1; i < num; ++i){
        res += Actxts[i];
    }
}

void CKKSmatrix::HEmatmul(Ctxt& res, Ctxt& Actxt, Ctxt& Bctxt, vector<vector<EncodedPtxt>>& Initpoly, vector<EncodedPtxt>& shiftpoly){
    vector<Ctxt> Actxts(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    vector<Ctxt> Bctxts(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
   
    //! 1. Generate the initial ciphertexts
    genInitCtxt(Actxts[0], Bctxts[0], Actxt, Bctxt, Initpoly);

    //! 2. Column shifting of Actxt[0], Row shifting of Bctxt[0]
    NTL_EXEC_RANGE(ckksmatpar.dim1, first, last);
    for(long i = first; i < last; ++i){
        long i1 = (i + 1);
        shiftBycols(Actxts[i1], Actxts[0], i1, shiftpoly);
        Bctxts[i1] = Bctxts[0];
        rotate(Bctxts[i1], -(ckksmatpar.dim * (i1)));
    }
    NTL_EXEC_RANGE_END;
    
    //! 3. Hadamard multiplication
    HEmatmul_Hadamard(res, Actxts, Bctxts, ckksmatpar.dim);
}

void CKKSmatrix::HEmatmul_Parallel(Ctxt& res, Ctxt& Actxt, Ctxt& Bctxt, vector<vector<EncodedPtxt>>& Initpoly, vector<EncodedPtxt>& shiftpoly){
    vector<Ctxt> Actxts(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    vector<Ctxt> Bctxts(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));
    
    //! 1. Generate the initial ciphertexts
    genInitCtxt_Parallel(Actxts[0], Bctxts[0], Actxt, Bctxt, Initpoly);

    //! 2. Column shifting of Actxt[0], Row shifting of Bctxt[0]
    long unit = ckksmatpar.dim  * ckksmatpar.nbatching;
    NTL_EXEC_RANGE(ckksmatpar.dim1, first, last);
    for(long i = first; i < last; ++i){
        long i1 = (i + 1);
        shiftBycols_Parallel(Actxts[i1], Actxts[0], i1, shiftpoly);
        Bctxts[i1] = Bctxts[0];
        rotate(Bctxts[i1], -(unit * (i1)));
    }
    NTL_EXEC_RANGE_END;
    
    //! 3. Hadamard mult : Actxts[0] * Bctxts[0] + ... + Actxts[d-1] * Bctxts[d-1]
    HEmatmul_Hadamard(res, Actxts, Bctxts, ckksmatpar.dim);
}

void CKKSmatrix::HErmatmul(Ctxt& res, Ctxt& Actxt, Ctxt& Bctxt, vector<vector<EncodedPtxt>>& Initpoly, vector<EncodedPtxt>& shiftpoly){
    vector<Ctxt> Actxts(ckksmatpar.subdim, Ctxt(ckksmeta.data->publicKey));
    vector<Ctxt> Bctxts(ckksmatpar.subdim, Ctxt(ckksmeta.data->publicKey));
    
    //! 1. Generate the initial ciphertexts
    genInitCtxt(Actxts[0], Bctxts[0], Actxt, Bctxt, Initpoly);
    
    //! 2. Column shifting of Actxt[0], Row shifting of Bctxt[0]
    NTL_EXEC_RANGE(ckksmatpar.subdim - 1, first, last);
    for(long i = first; i < last; ++i){
        long i1 = (i + 1);
        shiftBycols(Actxts[i1], Actxts[0], i1, shiftpoly);
        Bctxts[i1] = Bctxts[0];
        //equal to shiftByrows with i1
        rotate(Bctxts[i1], -(ckksmatpar.dim * (i1)));
    }
    NTL_EXEC_RANGE_END;
 
    //! 3. Hadamard mult : Actxts[0] * Bctxts[0] + ... + Actxts[d-1] * Bctxts[d-1]
    HEmatmul_Hadamard(res, Actxts, Bctxts, ckksmatpar.subdim);
    
    //! 4. shift and aggregate the results
    long index = (long) log2(ckksmatpar.dim/ckksmatpar.subdim);
    
    for(long i = 0; i < index; ++i){
        Ctxt ctemp(res);
        rotate(ctemp, -(ckksmatpar.dim * ckksmatpar.subdim * (1<<i)));
        res += ctemp;
    }
}

void CKKSmatrix::HEmatmul_preprocessing(Ctxt& res, vector<Ctxt>& Actxts, Ctxt& Bctxt, vector<EncodedPtxt>& Initpoly){
    //! 1. Generate the initial ciphertexts
    vector<Ctxt> Bctxts(ckksmatpar.dim, Ctxt(ckksmeta.data->publicKey));;
    genInitBctxt(Bctxts[0], Bctxt, Initpoly);
    
    //! 2. Row shifting of Bctxt[0]
    NTL_EXEC_RANGE(ckksmatpar.dim1, first, last);
    for(long i = first; i < last; ++i){
        long i1 = (i + 1);
        Bctxts[i1] = Bctxts[0];
        rotate(Bctxts[i1], -(ckksmatpar.dim * (i1)));
    }
    NTL_EXEC_RANGE_END;
    
    //! 3. Hadamard mult : Actxts[0] * Bctxts[0] + ... + Actxts[d-1] * Bctxts[d-1]
    NTL_EXEC_RANGE(ckksmatpar.dim, first, last);
    for(long i = first; i < last; ++i){
        Actxts[i].multiplyBy(Bctxts[i]);
    }
    NTL_EXEC_RANGE_END;
    
    res = Actxts[0];
    for(int i = 1; i < ckksmatpar.dim; ++i){
        res += Actxts[i];
    }
}

void CKKSmatrix::HErmatmul_preprocessing(Ctxt& res, vector<Ctxt>& Actxts, Ctxt& Bctxt, vector<EncodedPtxt>& Initpoly){
    //! 1. Generate the initial ciphertexts
    vector<Ctxt> Bctxts(ckksmatpar.subdim, Ctxt(ckksmeta.data->publicKey));
    genInitBctxt(Bctxts[0], Bctxt, Initpoly);
    
    //! 2. Row shifting of Bctxt[0]
    NTL_EXEC_RANGE(ckksmatpar.subdim - 1, first, last);
    //Fixme: int i or long i
    for(long i = first; i < last; ++i){
        long i1 = (i + 1);
        Bctxts[i1] = Bctxts[0];
        rotate(Bctxts[i1], -(ckksmatpar.dim * (i1)));
    }
    NTL_EXEC_RANGE_END;
    
    //! 3. Hadamard mult : Actxts[0] * Bctxts[0] + ... + Actxts[d-1] * Bctxts[d-1]
    NTL_EXEC_RANGE(ckksmatpar.subdim, first, last);
    for(long i = first; i < last; ++i){
        Actxts[i].multiplyBy(Bctxts[i]);
    }
    NTL_EXEC_RANGE_END;
    
    //! aggregate the results
    res = Actxts[0];
    for(int i = 1; i < ckksmatpar.subdim; ++i){
        res += Actxts[i];
    }
    
    //! 4. shift and aggregate the results
    long index = (long) log2(ckksmatpar.dim/ckksmatpar.subdim);
    for(long i = 0; i < index; ++i){
        Ctxt ctemp = res;
        rotate(ctemp, -(ckksmatpar.dim * ckksmatpar.subdim * (1<<i)));
        res += ctemp;
    }
}
