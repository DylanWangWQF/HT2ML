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
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <helib/helib.h>

#include "BGVmatrix.h"


void readHEMatpar(HEMatpar& HEmatpar, long nrows, long ncols, long subdim, long nbatching){
    HEmatpar.nrows = nrows;
    HEmatpar.ncols = ncols;
        
    long dim = nrows;
    if(dim < ncols) {
        dim = ncols;
    }
    
    HEmatpar.dim = (1<< (long)ceil(log2(dim)));  //! power of two
    HEmatpar.dim1 = HEmatpar.dim - 1;
    HEmatpar.logdim = (long) log2(HEmatpar.dim);
    HEmatpar.sqrdim = (long) ceil(sqrt(HEmatpar.dim));
    HEmatpar.subdim = subdim;  //! used for non-squared matrix multiplication
    HEmatpar.nbatching = (1<< (long)ceil(log2(nbatching)));
    
    if(nbatching == 1){
        HEmatpar.neededslots = HEmatpar.dim * HEmatpar.dim;  // just encode a single matrix
    }
    else{
        HEmatpar.neededslots = (HEmatpar.dim * HEmatpar.dim) * HEmatpar.nbatching;
    }
}

/*
 Note: In HElib, it only supports vector<long> in functions like: encrypt(), encode(). So we cannot use the long* as the input params.
 */

void HEmatrix::encryptZmat(Ctxt& ctxt, Mat<long>& mat){
    
    if (meta.data->ea.size() < HEmatpar.neededslots) {
        cout << "Number of slots in Ptxt less than needed slots in matrix, please re-FindM!" << endl;
        return;
    }
    
    // note: re-write the encryption of matrix for BGV
    //Ptxt<BGV> ptxt(meta.data->context);
    //Fixme: ea.size() != HEmatpar.neededslots
    vector<long> cmsg(meta.data->ea.size());
    
    NTL_EXEC_RANGE(HEmatpar.nrows, first, last);
    //Fixme: int i maybe enough as in HEmat
    for(long i = first; i < last; ++i){
        for(long j = 0; j < HEmatpar.ncols; ++j){
//            long temp;
//            conv(temp, mat[i][j]);
            cmsg[i * HEmatpar.dim + j] = mat[i][j];
            //ptxt[i * HEmatpar.dim  + j] = temp;
        }
    }
    NTL_EXEC_RANGE_END;

    //meta.data->publicKey.Encrypt(ctxt, ptxt);
    meta.data->ea.encrypt(ctxt, meta.data->publicKey, cmsg);
    //vector<long>().swap(cmsg);
    //cmsg.clear();
    //cmsg.shrink_to_fit();
}

void HEmatrix::decryptZmat(Mat<long>& mat, Ctxt& ctxt){
    //Ptxt<BGV> ptxt_result(meta.data->context);
    //meta.data->secretKey.Decrypt(ptxt_result, ctxt);
    vector<long> cmsg;
    meta.data->ea.decrypt(ctxt, meta.data->secretKey, cmsg);
    mat.SetDims(HEmatpar.dim, HEmatpar.dim);
    
    long k = 0;
    for(long i = 0; i < HEmatpar.dim; ++i){
        for(long j = 0; j < HEmatpar.dim; ++j){
            mat[i][j] = cmsg[k];
            k++;
        }
    }
}

void HEmatrix::encryptParallelZmat(Ctxt& ctxt, Mat<long>*& mat, long nbatching){
    if (meta.data->ea.size() < HEmatpar.neededslots) {
        cout << "Number of slots in Ptxt less than needed slots in matrix, please re-FindM!" << endl;
        return;
    }
    //Ptxt<BGV> ptxt(meta.data->context);
    //Note: Consider replace meta.data->ea.size() with HEmatpar.neededslots
    vector<long> cmsg(meta.data->ea.size());
//    long* ltemp = new long[nbatching];
    
    NTL_EXEC_RANGE(nbatching, first, last);
    for(long k = first; k < last; ++k){
        // encode the kth matrix
        for(long i = 0; i < HEmatpar.nrows; ++i){
            for(long j = 0; j < HEmatpar.ncols; ++j){
//                conv(ltemp[k], mat[k][i][j]);
                cmsg[(i * HEmatpar.dim + j) * HEmatpar.nbatching + k] = mat[k][i][j];
                //cmsg[(i * HEmatpar.dim + j) + (k * HEmatpar.dim * HEmatpar.dim)] = ltemp[k];
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    //meta.data->publicKey.Encrypt(ctxt, ptxt);
    meta.data->ea.encrypt(ctxt, meta.data->publicKey, cmsg);
    
    //cmsg.clear();
    //cmsg.shrink_to_fit();
//    delete[] ltemp;
}

void HEmatrix::decryptParallelZmat(Mat<long>*& mat, Ctxt& ctxt){
    //Ptxt<BGV> ptxt_result(meta.data->context);
    //meta.data->secretKey.Decrypt(ptxt_result, ctxt);
    vector<long> cmsg;
    meta.data->ea.decrypt(ctxt, meta.data->secretKey, cmsg);
    mat = new Mat<long>[HEmatpar.nbatching];
    
    for(long k = 0; k < HEmatpar.nbatching; ++k){
        mat[k].SetDims(HEmatpar.dim, HEmatpar.dim);
        long l = 0;
        for(long i = 0; i < HEmatpar.dim; ++i){
            for(long j = 0; j < HEmatpar.dim; ++j){
                mat[k][i][j] = cmsg[l * HEmatpar.nbatching + k];
                //mat[k][i][j] = cmsg[l + (k * HEmatpar.dim * HEmatpar.dim)];
                l++;
            }
        }
    }
    //cmsg.clear();
    //cmsg.shrink_to_fit();
}

//Fixme:ea.size() or HEmatpar.neededslots
void HEmatrix::msgleftRotate(vector<long>& res, vector<long> vals, long dim, long nrot){
    long nshift = (nrot) % HEmatpar.neededslots;
    //long nshift = (nrot) % meta.data->ea.size();
    
    long k = dim - nshift;
    for(long j = 0; j < k; ++j){
        res[j] = vals[j + nshift];
    }
    for(long j = k; j < dim; ++j){
        res[j] = vals[j - k];
    }
}

//void HEmatrix::msgleftRotate_Mul(vector<long>& res, vector<long> vals, long dim, long nrot){
//    //dim = HEmatpar.neededslots
//    long nshift = (nrot) % HEmatpar.neededslots;
//    //long nshift = (nrot) % meta.data->ea.size();
//
//    long k = dim - nshift;
//    for(long j = 0; j < k; ++j){
//        res[j] = vals[j + nshift];
//    }
//    for(long j = k; j < dim; ++j){
//        res[j] = vals[j - k];
//    }
//    for (long i = dim; i < meta.data->ea.size(); i++) {
//        res[i] = vals[i];
//    }
//}

void HEmatrix::msgrightRotate(vector<long>& res, vector<long> vals, long dim, long nrot){
    long nshift = (nrot) % HEmatpar.neededslots;
    //long nshift = (nrot) % meta.data->ea.size();
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

//void HEmatrix::msgrightRotate_Mul(vector<long>& res, vector<long> vals, long dim, long nrot){
//    //dim = HEmatpar.neededslots
//    long nshift = (nrot) % HEmatpar.neededslots;
//    //long nshift = (nrot) % meta.data->ea.size();
//
//    long k = dim - nshift;
//    for(long j = 0; j < k; ++j){
//        res[nshift + j] = vals[j];
//    }
//    for(long j = k; j < dim; ++j){
//        res[j - k] = vals[j];
//    }
//    for (long i = dim; i < meta.data->ea.size(); i++) {
//        res[i] = vals[i];
//    }
//}

void HEmatrix::msgleftRotateAndEqual(vector<long>& vals, long dim, long nrot){
    vector<long> res(dim);

    long nshift = (nrot) % HEmatpar.neededslots;
    //long nshift = (nrot) % meta.data->ea.size();
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

void HEmatrix::msgrightRotateAndEqual(vector<long>& vals, long dim, long nrot){
    vector<long> res(dim);
    
    long nshift = (nrot) % HEmatpar.neededslots;
    //long nshift = (nrot) % meta.data->ea.size();
    //Fixauthor, should be dim - nshift?
    long k = dim - nrot;
    //long k = dim - nshift;
    
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

//------------------------------------------------
//! Transposition
//------------------------------------------------

//! Output: transpoly ("2 * dim"), generate the polynomial for transpose
//! Originally we need to generate
//! rho(p[k]; -(d-1)*sqrt(r)i) where k = sqrt(r) * i + j, r = nrows, 0 <= i, j < sqrt(r)
// (0, 1,   ..., d-1),  (-, d+1, ..., 2d-1)
void HEmatrix::genTransPoly(zzX*& transpoly){
    vector<vector<long>> lvals(HEmatpar.dim); // define 1-dim's size
    vector<vector<long>> rvals(HEmatpar.dim);
    long dsquare = HEmatpar.neededslots - 1;
    transpoly = new zzX[2 * HEmatpar.dim];
    long dimsqrdim = HEmatpar.sqrdim * HEmatpar.dim1; //! (d-1) * sqrt(d)
    
    // k = i * sqrt(d) + j < d
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){
        btmp = false;   //! all the terms have the same numbers
    }
    else{
        btmp = true;
    }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
    
    NTL_EXEC_RANGE(ibound, first, last);
    for(long i = first; i < last; ++i){
        long jbound = HEmatpar.sqrdim;
        if ((btmp) && (i == ibound - 1)){   //! last term
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        
        for(long j = 0; j < jbound; ++j){
            long k = i * HEmatpar.sqrdim + j;
            
            //Fix1
//            lvals[k] = vector<long>(HEmatpar.neededslots);
//            rvals[k] = vector<long>(HEmatpar.neededslots);
            lvals[k] = vector<long>(meta.data->ea.size());
            rvals[k] = vector<long>(meta.data->ea.size());
            
            for(long l = 0; l < HEmatpar.dim - k; ++l){
                long ltemp = l * (HEmatpar.dim + 1) + k;
                lvals[k][ltemp] = 1;
                rvals[k][dsquare - ltemp] = 1;
            }
            
            //Fix2
            msgrightRotateAndEqual(lvals[k], HEmatpar.neededslots, i * dimsqrdim);
            msgleftRotateAndEqual(rvals[k], HEmatpar.neededslots, i * dimsqrdim);
            
            meta.data->ea.encode(transpoly[k], lvals[k]);
            meta.data->ea.encode(transpoly[k + HEmatpar.dim], rvals[k]);
        }
    }
    NTL_EXEC_RANGE_END;
}

void HEmatrix::genTransPoly_Parallel(zzX*& transpoly){
    
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){
        btmp = false;
    }    //! all the terms have the same numbers
    else{
        btmp = true;
    }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
    
    vector<vector<long>> lvals(HEmatpar.dim); // define 1-dim's size
    vector<vector<long>> rvals(HEmatpar.dim);
    
    transpoly = new zzX[2 * HEmatpar.dim];
    
    long shiftunit = HEmatpar.sqrdim * HEmatpar.dim1 * HEmatpar.nbatching ;   //! (d-1) * sqrt(d) * l
    long dsquare = ((HEmatpar.dim * HEmatpar.dim) - 1) * HEmatpar.nbatching ;  //! (d^2-1) * l
    
    //! k = (i * sqrt(d) + j)  < d
    //! poly[k] = poly[i * sqrt(d) + j]  -> rho(poly[k]; i * (d-1)*sqrt(d) * l)
    for(int i = 0; i < ibound; ++i){
        long jbound = HEmatpar.sqrdim;
        if ((btmp)&&(i == ibound - 1)){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        for(long j = 0; j < jbound; ++j){
            long k = (i * HEmatpar.sqrdim + j);
        
//            lvals[k] = vector<long>(HEmatpar.neededslots);
//            rvals[k] = vector<long>(HEmatpar.neededslots);
            lvals[k] = vector<long>(meta.data->ea.size());
            rvals[k] = vector<long>(meta.data->ea.size());
            
            for(long l = 0; l < HEmatpar.dim - k; ++l){
                long dtemp= (l * (HEmatpar.dim + 1) + k) * HEmatpar.nbatching;  //! starting index
        
                for(long n = 0; n < HEmatpar.nbatching; ++n){
                    lvals[k][dtemp + n] = 1;
                    rvals[k][dsquare - dtemp + n] = 1;
                }
            }
        
            //! Lrho(P[k], - i * (d-1) * sqrt(d) * (nbathcing))
            //! Rrho(P[k],- i * (d-1) * sqrt(d) * (nbathcing))
            msgrightRotateAndEqual(lvals[k], HEmatpar.neededslots, i * shiftunit);
            msgleftRotateAndEqual(rvals[k], HEmatpar.neededslots, i * shiftunit);
            
            meta.data->ea.encode(transpoly[k], lvals[k]);
            meta.data->ea.encode(transpoly[k + HEmatpar.dim], rvals[k]);
        }
    }
}

void HEmatrix::transpose(Ctxt& res, Ctxt& ctxt, zzX*& transpoly){
    
//    res = ctxt;
//    res.multByConstant(transpoly[HEmatpar.dim + 5]);
    
//    vector<Ctxt> Babyctxt1(HEmatpar.sqrdim, Ctxt(meta.data->publicKey));
//    Babyctxt1[1] = ctxt;
//    rotate(Babyctxt1[1], (-1 * (HEmatpar.dim1)));
//    vector<Ctxt> Babyctxt2(HEmatpar.sqrdim, Ctxt(meta.data->publicKey));
//    Babyctxt2[1] = ctxt;
//    rotate(Babyctxt2[1], (1 * (HEmatpar.dim1)));
//    res = Babyctxt2[1];
    
    long dimsqrdim = HEmatpar.sqrdim * HEmatpar.dim1;   // (d-1) * sqrt(d)
    
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){
        btmp = false;   //! all the terms have the same numbers
    }
    else{
        btmp = true;
    }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
    
    vector<Ctxt> Babyctxt1(HEmatpar.sqrdim, Ctxt(meta.data->publicKey));
    vector<Ctxt> Babyctxt2(HEmatpar.sqrdim, Ctxt(meta.data->publicKey));
    
    vector<vector<Ctxt>> ltemp(HEmatpar.sqrdim, vector<Ctxt>(HEmatpar.sqrdim, Ctxt(meta.data->publicKey)));
    vector<vector<Ctxt>> rtemp(HEmatpar.sqrdim, vector<Ctxt>(HEmatpar.sqrdim, Ctxt(meta.data->publicKey)));
    
    //! diagonal
    //ctxt.multByConstant(transpoly[0]); this changes the ctxt
    res = ctxt;
    HELIB_NTIMER_START(multiply_transpoly0);
    res.multByConstant(transpoly[0]);
    HELIB_NTIMER_STOP(multiply_transpoly0);
    printNamedTimer(cout, "multiply_transpoly0");
    
    Babyctxt1[0] = ctxt;
    Babyctxt2[0] = ctxt;
    
    //! Babyctxt1[j] = rho(v; j * (d-1))
    //! res[j] = rho(-, (d-1)*srt(d) * i) * rho(v; j * (d-1))
    NTL_EXEC_RANGE(HEmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);  // 1 <= j1 <= sqrdim - 1
        
        // prepare baby ctxt
        //-:leftrotate,+:rightrotate
        cout << j << "-th rotate start----------------------------" << endl;
        HELIB_NTIMER_START(rotateBabyctxt1);
        Babyctxt1[j1] = ctxt;
        rotate(Babyctxt1[j1], (-j1 * (HEmatpar.dim1)));
        HELIB_NTIMER_STOP(rotateBabyctxt1);
        printNamedTimer(cout, "rotateBabyctxt1");
        
        HELIB_NTIMER_START(rotateBabyctxt2);
        Babyctxt2[j1] = ctxt;
        rotate(Babyctxt2[j1], j1 * (HEmatpar.dim1));
        HELIB_NTIMER_STOP(rotateBabyctxt2);
        printNamedTimer(cout, "rotateBabyctxt2");
        cout << j << "-th rotate end----------------------------" << endl;
        
        cout << endl;

        HELIB_NTIMER_START(multiply_ltemp);
        ltemp[0][j1] = Babyctxt1[j1];
        ltemp[0][j1].multByConstant(transpoly[j1]);
        HELIB_NTIMER_STOP(multiply_ltemp);
        printNamedTimer(cout, "multiply_ltemp");
        
        HELIB_NTIMER_START(multiply_rtemp);
        rtemp[0][j1] = Babyctxt2[j1];
        rtemp[0][j1].multByConstant(transpoly[j1 + HEmatpar.dim]);
        HELIB_NTIMER_STOP(multiply_rtemp);
        printNamedTimer(cout, "multiply_rtemp");
        
        ltemp[0][j1] += rtemp[0][j1];
        
        cout << endl;
    }
    NTL_EXEC_RANGE_END;

    for(long j = 1; j < HEmatpar.sqrdim; ++j){
        res += ltemp[0][j];
    }

    //! 1 <= i < d
    //! [k]: rho(p[k];-dsqrt(d)*i1) * rho(v, j*(d-1)) -> rot by dsqrt(d)*i
    //! res[j] = rho( rho(-, (d-1)*srt(d) * i) * rho(v; j * (d-1)) ; i * (d-1)\sqrt(d) )
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long i = first; i < last; ++i){
        long jbound = HEmatpar.sqrdim;
        if (btmp && (i == (ibound - 2))){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        
        long i1 = i + 1;
        long k = (i1) * HEmatpar.sqrdim;

        //! 0 <= j < sqrdim
        NTL_EXEC_RANGE(jbound, first, last);
        for(long j = first; j < last; ++j){
            ltemp[i1][j] = Babyctxt1[j];
            ltemp[i1][j].multByConstant(transpoly[k + j]);
            rtemp[i1][j] = Babyctxt2[j];
            rtemp[i1][j].multByConstant(transpoly[k + j + HEmatpar.dim]);
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

void HEmatrix::transpose_Parallel(Ctxt& res, Ctxt& ctxt, zzX*& transpoly){
    long shiftunit = HEmatpar.sqrdim * HEmatpar.dim1 * HEmatpar.nbatching;  // (d-1) * sqrt(d) * l
    long unit = (HEmatpar.dim1) * HEmatpar.nbatching;
    
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){
        btmp = false;
    }    //! all the terms have the same numbers
    else{
        btmp = true;
    }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
    
    vector<Ctxt> Babyctxt1(HEmatpar.sqrdim, Ctxt(meta.data->publicKey));
    vector<Ctxt> Babyctxt2(HEmatpar.sqrdim, Ctxt(meta.data->publicKey));
    
    vector<vector<Ctxt>> ltemp(HEmatpar.sqrdim, vector<Ctxt>(HEmatpar.sqrdim, Ctxt(meta.data->publicKey)));
    vector<vector<Ctxt>> rtemp(HEmatpar.sqrdim, vector<Ctxt>(HEmatpar.sqrdim, Ctxt(meta.data->publicKey)));
    
    //-------------------------
    // i = 0, sqrt(d) polynomials
    //-------------------------
    res = ctxt;
    res.multByConstant(transpoly[0]);

    Babyctxt1[0] = ctxt;
    Babyctxt2[0] = ctxt;
    
    NTL_EXEC_RANGE(HEmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);  //! 1 <= j1 <= sqrdim - 1
        
        // rho(v; j * (d-1) * l)
        Babyctxt1[j1] = ctxt;
        rotate(Babyctxt1[j1], (-(j1 * unit)));
        Babyctxt2[j1] = ctxt;
        rotate(Babyctxt2[j1], (j1 * unit));
        
        ltemp[0][j1] = Babyctxt1[j1];
        ltemp[0][j1].multByConstant(transpoly[j1]);
        rtemp[0][j1] = Babyctxt2[j1];
        rtemp[0][j1].multByConstant(transpoly[j1 + HEmatpar.dim]);
        
        ltemp[0][j1] += rtemp[0][j1];
    }
    NTL_EXEC_RANGE_END;
    
    for(long j = 1; j < HEmatpar.sqrdim; ++j){
        res += ltemp[0][j];
    }
    
    //! 1 <= i < d
    //! [k]: rho(p[k];-dsqrt(d)*i1) * rho(v, j*(d-1)) -> rot by dsqrt(d)*i
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long i = first; i < last; ++i){
        long jbound = HEmatpar.sqrdim;
        if ((btmp)&&(i == ibound - 2)){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        
        long i1 = i + 1;
        long k = (i1) * HEmatpar.sqrdim;
        
        NTL_EXEC_RANGE(jbound, first, last);
        for(long j = first; j < last; ++j){
            ltemp[i1][j] = Babyctxt1[j];
            ltemp[i1][j].multByConstant(transpoly[k + j]);
            rtemp[i1][j] = Babyctxt2[j];
            rtemp[i1][j].multByConstant(transpoly[k + j + HEmatpar.dim]);
        }
        NTL_EXEC_RANGE_END;
        
        for(long j = 1; j < jbound; ++j){
            ltemp[i1][0] += ltemp[i1][j];
            rtemp[i1][0] += rtemp[i1][j];
        }
        
        rotate(ltemp[i1][0], -(i1 * shiftunit));
        rotate(rtemp[i1][0], i1 * shiftunit);

        ltemp[i1][0] += rtemp[i1][0];
    }
    NTL_EXEC_RANGE_END;
    
    for(long i = 1; i < ibound; ++i){
        res += ltemp[i][0];
    }
}

//------------------------------------------------
//! Shift
//------------------------------------------------

// num = 0: dimension d-1
// subdim - 1: dimension (subdim-1)
void HEmatrix::genShiftPoly(zzX*& shiftpoly, long num){
    long length = HEmatpar.dim1;
    if(num != 0){
        length = num;
    }
    
    vector<vector<long>> vals(length, vector<long>(HEmatpar.neededslots));
//    vector<vector<long>> vals(length, vector<long>(meta.data->ea.size()));
    shiftpoly = new zzX[length];
    
    // i: shifted by (i+1)
    NTL_EXEC_RANGE(length, first, last);
    for(long i = first; i < last; ++i){
        for(long j = 0; j < HEmatpar.dim; ++j){
            for(long k = 0; k < i + 1; ++k){
                vals[i][j * HEmatpar.dim + k] = 1;
            }
        }
        meta.data->ea.encode(shiftpoly[i], vals[i]);
    }
    NTL_EXEC_RANGE_END;
//    cout << "shift_vals_vector: " << endl;
//    cout << vals << endl;
}

void HEmatrix::genShiftPoly_Parallel(zzX*& shiftpoly){
    vector<vector<long>> vals(HEmatpar.dim1, vector<long>(HEmatpar.neededslots));
//    vector<vector<long>> vals(HEmatpar.dim1, vector<long>(meta.data->ea.size()));
    shiftpoly = new zzX[HEmatpar.dim1];
    
    NTL_EXEC_RANGE(HEmatpar.dim1, first, last);
    for(long i = first; i < last; ++i){
//        vals[i] = vector<long>(HEmatpar.neededslots);
        
        for(long j = 0; j < HEmatpar.dim; ++j){
            for(long k = 0; k < i + 1; ++k){
                long dtemp = (j * HEmatpar.dim + k) * HEmatpar.nbatching;
                
                for(long n = 0; n < HEmatpar.nbatching; ++n){
                    vals[i][dtemp + n] = 1;
                }
            }
        }
        meta.data->ea.encode(shiftpoly[i], vals[i]);
    }
    NTL_EXEC_RANGE_END;
}

void HEmatrix::shiftBycols(Ctxt& res, Ctxt& ctxt,  long k , zzX*& shiftpoly){
    //! ctemp = Enc(m[1],.., m[k],0,,,0 | ....)
    Ctxt ctemp(ctxt);
    ctemp.multByConstant(shiftpoly[k-1]);
    
    //! ctemp = Enc(0,,,0, m[k+1],...,m[d] | ....)
    
    //TBD1: do we need to ModDownToSet/BringToSet for res?
    //TBD2: or ModDown is needed for accelerating?
    res = ctxt;
    res -= ctemp;
    rotate(ctemp, (HEmatpar.dim - k));
    rotate(res, -k);
    res += ctemp;
}

void HEmatrix::shiftBycols_Parallel(Ctxt& res, Ctxt& ctxt, long k, zzX*& shiftpoly){
    //! ctemp = Enc(m[1],.., m[k],0,,,0 | ....)
    Ctxt ctemp (ctxt);
    ctemp.multByConstant(shiftpoly[k-1]);
    res = ctxt;
    res -= ctemp;
    rotate(ctemp, (HEmatpar.dim - k) * HEmatpar.nbatching);
    rotate(res, -(k * HEmatpar.nbatching));
    res += ctemp;
}

//------------------------------------------------
//! Matrix multiplication
//------------------------------------------------

//!@ Output: Initpoly[0] (constant left-polynomials for Amat): 0,1,....,d-1
//!@              Initpoly[1] (constant right-polynomials for Amat):  -1,...,-(d-1)
//!!                 Initpoly[2] (constant polynomials for Bmat)

void HEmatrix::genMultPoly(zzX**& Initpoly){
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){
        btmp = false;
    }
    else{
        btmp = true;
    }
    long ibound = (long) ceil((double) HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
    
    Initpoly = new zzX*[3];
    
    Initpoly[0] =  new zzX[HEmatpar.dim];
    Initpoly[1] =  new zzX[HEmatpar.dim];
    Initpoly[2] =  new zzX[HEmatpar.dim];
    
    vector<vector<long>> fvals1(HEmatpar.dim);
    vector<vector<long>> fvals2(HEmatpar.dim);
    vector<vector<long>> bvals(HEmatpar.dim);

    NTL_EXEC_RANGE(ibound, first, last);
    for(long i = first; i < last; ++i){
        long jbound = HEmatpar.sqrdim;
        if ((btmp)&&(i == ibound - 1)){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        
        for(long j = 0; j < jbound; ++j){
            long k = i * HEmatpar.sqrdim + j;
            
//            fvals1[k] = vector<long>(meta.data->ea.size());
//            fvals2[k] = vector<long>(meta.data->ea.size());
//            bvals[k] = vector<long>(meta.data->ea.size());
            fvals1[k] = vector<long>(HEmatpar.neededslots);
            fvals2[k] = vector<long>(HEmatpar.neededslots);
            bvals[k] = vector<long>(HEmatpar.neededslots);
            
            for(long l = 0; l < HEmatpar.dim - k; ++l){
                fvals1[k][k * HEmatpar.dim + l] = 1;
            }
            
            msgleftRotate(fvals2[k], fvals1[k], HEmatpar.neededslots, k * (2 * HEmatpar.dim - 1));
            
            msgrightRotateAndEqual(fvals1[k], HEmatpar.neededslots, i * HEmatpar.sqrdim);
            meta.data->ea.encode(Initpoly[0][k], fvals1[k]);
        
            msgleftRotateAndEqual(fvals2[k], HEmatpar.neededslots, i * HEmatpar.sqrdim);
            meta.data->ea.encode(Initpoly[1][k], fvals2[k]);
            
            for(long l = 0; l < HEmatpar.dim; ++l){
                bvals[k][l * HEmatpar.dim + k] = 1;
            }
            
            msgrightRotateAndEqual(bvals[k], HEmatpar.neededslots, i * HEmatpar.sqrdim * HEmatpar.dim);
            meta.data->ea.encode(Initpoly[2][k], bvals[k]);
        }
    }
    NTL_EXEC_RANGE_END;
}

void HEmatrix::genMultPoly_Parallel(zzX**& Initpoly){
    long btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){ btmp = false; }    //! all the terms have the same numbers
    else{ btmp = true; }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"

    
    Initpoly = new zzX*[3];
    
    Initpoly[0] =  new zzX[HEmatpar.dim];
    Initpoly[1] =  new zzX[HEmatpar.dim];
    Initpoly[2] =  new zzX[HEmatpar.dim];
    
    vector<vector<long>> fvals1(HEmatpar.dim);
    vector<vector<long>> fvals2(HEmatpar.dim);
    vector<vector<long>> bvals(HEmatpar.dim);
    
    
    NTL_EXEC_RANGE(ibound, first, last);
    for(long i = first; i < last; ++i){
        long jbound = HEmatpar.sqrdim;
        if ((btmp) && (i == (ibound - 1))){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        
        for(long j = 0; j < jbound; ++j){
            long k = i * HEmatpar.sqrdim + j;
            
//            fvals1[k] = vector<long>(HEmatpar.neededslots);
//            fvals2[k] = vector<long>(HEmatpar.neededslots);
//            bvals[k] = vector<long>(HEmatpar.neededslots);
            fvals1[k] = vector<long>(meta.data->ea.size());
            fvals2[k] = vector<long>(meta.data->ea.size());
            bvals[k] = vector<long>(meta.data->ea.size());
            
            //! original: k * d <= index <= k * d + d - k -1
            //! parallel: (k * d) * l + 0 <= ... <= (k * d) * l + (l-1)
            //!           (k * d + d - k -1) * l + 0 <= ... <= (k * d + d - k -1) * l + (l-1) <  (k * d + d - k -1) * l + l
            //!           from k * d * l <= index <  (k * d + d - k) * l
            long start = k * HEmatpar.dim * HEmatpar.nbatching;
            long end = ( k * HEmatpar.dim + HEmatpar.dim - k) * HEmatpar.nbatching;
            for(long l = start; l < end; ++l){
                fvals1[k][l] = 1;
            }
            
            msgleftRotate(fvals2[k], fvals1[k], HEmatpar.neededslots, k * (2 * HEmatpar.dim - 1) * HEmatpar.nbatching);
            
            msgrightRotateAndEqual(fvals1[k], HEmatpar.neededslots, i*HEmatpar.sqrdim * HEmatpar.nbatching);
            meta.data->ea.encode(Initpoly[0][k], fvals1[k]);
            
            msgleftRotateAndEqual(fvals2[k], HEmatpar.neededslots, i*HEmatpar.sqrdim * HEmatpar.nbatching);
            meta.data->ea.encode(Initpoly[1][k], fvals2[k]);
            
            for(long l = 0; l < HEmatpar.dim; ++l){
                long dtemp = (l * HEmatpar.dim + k) * HEmatpar.nbatching;
                
                for(long n = 0; n < HEmatpar.nbatching; ++n){
                    bvals[k][dtemp + n] = 1;
                }
            }
            msgrightRotateAndEqual(bvals[k], HEmatpar.neededslots, i * HEmatpar.sqrdim * HEmatpar.dim * HEmatpar.nbatching);
            meta.data->ea.encode(Initpoly[2][k], bvals[k]);
        }
    }
    NTL_EXEC_RANGE_END;
}

void HEmatrix::genInitCtxt(Ctxt& resA, Ctxt& resB, Ctxt& Actxt, Ctxt& Bctxt, zzX**& poly){
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){
        btmp = false;
    }    //! all the terms have the same numbers
    else{
        btmp = true;
    }
    
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
    
    vector<vector<Ctxt>> Actemp1(HEmatpar.sqrdim, vector<Ctxt>(HEmatpar.sqrdim, Ctxt(meta.data->publicKey)));
    vector<vector<Ctxt>> Actemp2(HEmatpar.sqrdim, vector<Ctxt>(HEmatpar.sqrdim, Ctxt(meta.data->publicKey)));
    vector<vector<Ctxt>> Bctemp(HEmatpar.sqrdim, vector<Ctxt>(HEmatpar.sqrdim, Ctxt(meta.data->publicKey)));
    
    
    //! 0. Store some ciphertexts (0,1,...,d-1), (,d+1,...2d-1)
    //! v, lrho(v;1), lrho(v;2), ..., lrho(v;d-1)
    //! -, rrho(v;1), rrho(v;2), ..., rrho(v;d-1)
    
    vector<Ctxt> BaByctxt1(HEmatpar.dim, Ctxt(meta.data->publicKey));
    vector<Ctxt> BaByctxt2(HEmatpar.dim, Ctxt(meta.data->publicKey));
    vector<Ctxt> BaByctxtB(HEmatpar.dim, Ctxt(meta.data->publicKey));
    
    BaByctxt1[0] = Actxt;
    BaByctxt2[0] = Actxt;
    BaByctxtB[0] = Bctxt;
    
    //! i = 0:   Actxts[0] = v[0] + p1 * v[1] + ... + p[sqr(d)-1] *  v[sqr(d)-1]
    resA = BaByctxt1[0];
    resA.multByConstant(poly[0][0]);
    resB = BaByctxtB[0];
    resB.multByConstant(poly[2][0]);
    
    NTL_EXEC_RANGE(HEmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);
        BaByctxt1[j1] = Actxt;
        rotate(BaByctxt1[j1], -j1);
//        Ctxt ctemp1(BaByctxt1[j1]);
//        shift(ctemp1, -(meta.data->ea.size() - HEmatpar.neededslots));
//        BaByctxt1[j1] += ctemp1;
        
        BaByctxt2[j1] = Actxt;
        rotate(BaByctxt2[j1], j1);
//        rotate(BaByctxt2[j1], (j1 + meta.data->ea.size() - HEmatpar.neededslots));
//        Ctxt ctemp2(BaByctxt2[j1]);
//        shift(ctemp2, (meta.data->ea.size() - HEmatpar.neededslots));
//        BaByctxt2[j1] += ctemp2;
        
        BaByctxtB[j1] = Bctxt;
        rotate(BaByctxtB[j1], -((j1) * HEmatpar.dim));
        
        Actemp1[0][j1] = BaByctxt1[j1];
        Actemp1[0][j1].multByConstant(poly[0][j1]);
        Actemp2[0][j1] = BaByctxt2[j1];
        Actemp2[0][j1].multByConstant(poly[1][j1]);
        Bctemp[0][j1] = BaByctxtB[j1];
        Bctemp[0][j1].multByConstant(poly[2][j1]);
        
        Actemp1[0][j1] += Actemp2[0][j1];
    }
    NTL_EXEC_RANGE_END;

    for(long j = 1; j < HEmatpar.sqrdim; ++j){
        resA += Actemp1[0][j];
        resB += Bctemp[0][j];
    }
    
    vector<Ctxt> Actxts1(HEmatpar.dim, Ctxt(meta.data->publicKey));
    vector<Ctxt> Actxts2(HEmatpar.dim, Ctxt(meta.data->publicKey));
    vector<Ctxt> Bctxts(HEmatpar.dim, Ctxt(meta.data->publicKey));
    
    NTL_EXEC_RANGE(HEmatpar.dim - HEmatpar.sqrdim, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k + HEmatpar.sqrdim);
        long i = (long)(k1 / HEmatpar.sqrdim);
        long j = (long)(k1 % HEmatpar.sqrdim);
        
        Actxts1[k1] = BaByctxt1[j];
        Actxts1[k1].multByConstant(poly[0][k1]);
        Actxts2[k1] = BaByctxt2[j];
        Actxts2[k1].multByConstant(poly[1][k1]);
        Bctxts[k1] = BaByctxtB[j];
        Bctxts[k1].multByConstant(poly[2][k1]);
    }
    NTL_EXEC_RANGE_END;
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k + 1) * HEmatpar.sqrdim;
        long jbound = HEmatpar.sqrdim;
        if((btmp) && (k == ibound - 2)){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        for(long j = 1; j < jbound; ++j){
            Actxts1[k1] += Actxts1[k1+j];
            Actxts2[k1] += Actxts2[k1+j];
            Bctxts[k1] += Bctxts[k1+j];
        }
        
        long k2 = (k1 - (k1 % HEmatpar.sqrdim));
        rotate(Actxts1[k1], -k2);
        rotate(Actxts2[k1], k2);
        rotate(Bctxts[k1], -(k2 * HEmatpar.dim));
        
        Actxts1[k1] += Actxts2[k1];
    }
    NTL_EXEC_RANGE_END;
    
    for(long k = 1; k < ibound; ++k){
        long k1 = k * HEmatpar.sqrdim;
        resA += Actxts1[k1];
        resB += Bctxts[k1];
    }
}

void HEmatrix::genInitCtxt_Parallel(Ctxt& resA, Ctxt& resB, Ctxt& Actxt, Ctxt& Bctxt, zzX**& poly){
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){ btmp = false; }    //! all the terms have the same numbers
    else{ btmp = true; }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"

    
    vector<vector<Ctxt>> Actemp1(HEmatpar.sqrdim, vector<Ctxt>(HEmatpar.sqrdim, Ctxt(meta.data->publicKey)));
    vector<vector<Ctxt>> Actemp2(HEmatpar.sqrdim, vector<Ctxt>(HEmatpar.sqrdim, Ctxt(meta.data->publicKey)));
    vector<vector<Ctxt>> Bctemp(HEmatpar.sqrdim, vector<Ctxt>(HEmatpar.sqrdim, Ctxt(meta.data->publicKey)));

    //-----------------------------------------------------------
    // 0. Store some ciphertexts (0,1,...,d-1), ( ,d+1,...2d-1)
    // v, lrho(v;1), lrho(v;2), ..., lrho(v;d-1)
    // -, rrho(v;1), rrho(v;2), ..., rrho(v;d-1)
    
    vector<Ctxt> BaByctxt1(HEmatpar.dim, Ctxt(meta.data->publicKey));
    vector<Ctxt> BaByctxt2(HEmatpar.dim, Ctxt(meta.data->publicKey));
    vector<Ctxt> BaByctxtB(HEmatpar.dim, Ctxt(meta.data->publicKey));
    
    BaByctxt1[0] = Actxt;
    BaByctxt2[0] = Actxt;
    BaByctxtB[0] = Bctxt;
    
    resA = BaByctxt1[0];
    resA.multByConstant(poly[0][0]);
    resB = BaByctxtB[0];
    resB.multByConstant(poly[2][0]);
    
    NTL_EXEC_RANGE(HEmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);
        BaByctxt1[j1] = Actxt;
        rotate(BaByctxt1[j1], -(j1 * HEmatpar.nbatching));
        BaByctxt2[j1] = Actxt;
        rotate(BaByctxt2[j1], j1 * HEmatpar.nbatching);
        BaByctxtB[j1] = Bctxt;
        rotate(BaByctxtB[j1], -((j1) * HEmatpar.dim * HEmatpar.nbatching));
        
        Actemp1[0][j1] = BaByctxt1[j1];
        Actemp1[0][j1].multByConstant(poly[0][j1]);
        Actemp2[0][j1] = BaByctxt2[j1];
        Actemp2[0][j1].multByConstant(poly[1][j1]);
        Bctemp[0][j1] = BaByctxtB[j1];
        Bctemp[0][j1].multByConstant(poly[2][j1]);
        
        Actemp1[0][j1] += Actemp2[0][j1];
    }
    NTL_EXEC_RANGE_END;
    
    for(long j = 1; j < HEmatpar.sqrdim; ++j){
        resA += Actemp1[0][j];
        resB += Bctemp[0][j];
    }
 
    //---------------------------
    vector<Ctxt> Actxts1(HEmatpar.dim, Ctxt(meta.data->publicKey));
    vector<Ctxt> Actxts2(HEmatpar.dim, Ctxt(meta.data->publicKey));
    vector<Ctxt> Bctxts(HEmatpar.dim, Ctxt(meta.data->publicKey));
    
    NTL_EXEC_RANGE(HEmatpar.dim - ibound, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k + HEmatpar.sqrdim);
        long i = (long)(k1 / HEmatpar.sqrdim);
        long j = (long)(k1 % HEmatpar.sqrdim);
        
        Actxts1[k1] = BaByctxt1[j];
        Actxts1[k1].multByConstant(poly[0][k1]);
        Actxts2[k1] = BaByctxt2[j];
        Actxts2[k1].multByConstant(poly[1][k1]);
        Bctxts[k1] = BaByctxtB[j];
        Bctxts[k1].multByConstant(poly[2][k1]);
    }
    NTL_EXEC_RANGE_END;
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k + 1) * HEmatpar.sqrdim;
        long jbound = HEmatpar.sqrdim;
        if((btmp) && (k == ibound - 2)){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        for(long j = 1; j < jbound; ++j){
            Actxts1[k1] += Actxts1[k1+j];
            Actxts2[k1] += Actxts2[k1+j];
            Bctxts[k1] += Bctxts[k1+j];
        }
        
        long k2 = (k1 - (k1 % HEmatpar.sqrdim)) * HEmatpar.nbatching;
        rotate(Actxts1[k1], -k2);
        rotate(Actxts2[k1], k2);
        rotate(Bctxts[k1], -(k2 * HEmatpar.dim));
        
        Actxts1[k1] += Actxts2[k1];
    }
    NTL_EXEC_RANGE_END;
    
    
    for(long k = 1; k < ibound; ++k){
        long k1 = k * HEmatpar.sqrdim;
        resA += Actxts1[k1];
        resB += Bctxts[k1];
    }
}

void HEmatrix::HEmatmul_Hadamard(Ctxt& res, vector<Ctxt>& Actxts, vector<Ctxt>& Bctxts, long num){
    
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
    }
    NTL_EXEC_RANGE_END;
    
    for(int i = 1; i < num; ++i){
        res += Actxts[i];
    }
}

void HEmatrix::HEmatmul(Ctxt& res, Ctxt& Actxt, Ctxt& Bctxt, zzX**& Initpoly, zzX*& shiftpoly){
    vector<Ctxt> Actxts(HEmatpar.dim, Ctxt(meta.data->publicKey));
    vector<Ctxt> Bctxts(HEmatpar.dim, Ctxt(meta.data->publicKey));
   
    //! 1. Generate the initial ciphertexts
    genInitCtxt(Actxts[0], Bctxts[0], Actxt, Bctxt, Initpoly);

    //! 2. Column shifting of Actxt[0], Row shifting of Bctxt[0]
    NTL_EXEC_RANGE(HEmatpar.dim1, first, last);
    for(long i = first; i < last; ++i){
        long i1 = (i + 1);
        shiftBycols(Actxts[i1], Actxts[0], i1, shiftpoly);
        Bctxts[i1] = Bctxts[0];
        rotate(Bctxts[i1], -(HEmatpar.dim * (i1)));
    }
    NTL_EXEC_RANGE_END;
    
    //! 3. Hadamard multiplication
    HEmatmul_Hadamard(res, Actxts, Bctxts, HEmatpar.dim);
}

void HEmatrix::HEmatmul_Parallel(Ctxt& res, Ctxt& Actxt, Ctxt& Bctxt, zzX**& Initpoly, zzX*& shiftpoly){
    vector<Ctxt> Actxts(HEmatpar.dim, Ctxt(meta.data->publicKey));
    vector<Ctxt> Bctxts(HEmatpar.dim, Ctxt(meta.data->publicKey));
    
    //! 1. Generate the initial ciphertexts
    genInitCtxt_Parallel(Actxts[0], Bctxts[0], Actxt, Bctxt, Initpoly);

    //! 2. Column shifting of Actxt[0], Row shifting of Bctxt[0]
    long unit = HEmatpar.dim  * HEmatpar.nbatching;
    NTL_EXEC_RANGE(HEmatpar.dim1, first, last);
    for(long i = first; i < last; ++i){
        long i1 = (i + 1);
        shiftBycols_Parallel(Actxts[i1], Actxts[0], i1, shiftpoly);
        Bctxts[i1] = Bctxts[0];
        rotate(Bctxts[i1], -(unit * (i1)));
    }
    NTL_EXEC_RANGE_END;
    
    //! 3. Hadamard mult : Actxts[0] * Bctxts[0] + ... + Actxts[d-1] * Bctxts[d-1]
    HEmatmul_Hadamard(res, Actxts, Bctxts, HEmatpar.dim);
}

void HEmatrix::HErmatmul(Ctxt& res, Ctxt& Actxt, Ctxt& Bctxt, zzX**& Initpoly, zzX*& shiftpoly){
    vector<Ctxt> Actxts(HEmatpar.subdim, Ctxt(meta.data->publicKey));
    vector<Ctxt> Bctxts(HEmatpar.subdim, Ctxt(meta.data->publicKey));
    
    //! 1. Generate the initial ciphertexts
    genInitCtxt(Actxts[0], Bctxts[0], Actxt, Bctxt, Initpoly);
    
    //! 2. Column shifting of Actxt[0], Row shifting of Bctxt[0]
    NTL_EXEC_RANGE(HEmatpar.subdim - 1, first, last);
    for(long i = first; i < last; ++i){
        long i1 = (i + 1);
        shiftBycols(Actxts[i1], Actxts[0], i1, shiftpoly);
        Bctxts[i1] = Bctxts[0];
        //equal to shiftByrows with i1
        rotate(Bctxts[i1], -(HEmatpar.dim * (i1)));
    }
    NTL_EXEC_RANGE_END;
 
    //! 3. Hadamard mult : Actxts[0] * Bctxts[0] + ... + Actxts[d-1] * Bctxts[d-1]
    HEmatmul_Hadamard(res, Actxts, Bctxts, HEmatpar.subdim);
    
    //! 4. shift and aggregate the results
    long index = (long) log2(HEmatpar.dim/HEmatpar.subdim);
    
    for(long i = 0; i < index; ++i){
        //Fixme: can we use vector<Ctxt> instead of Ctxt outside the for loop?
        Ctxt ctemp(res);
        rotate(ctemp, -(HEmatpar.dim * HEmatpar.subdim * (1<<i)));
        res += ctemp;
    }
}

void HEmatrix::HErmatmul_Parallel(Ctxt& res, Ctxt& Actxt, Ctxt& Bctxt, zzX**& Initpoly, zzX*& shiftpoly){
    vector<Ctxt> Actxts(HEmatpar.subdim, Ctxt(meta.data->publicKey));
    vector<Ctxt> Bctxts(HEmatpar.subdim, Ctxt(meta.data->publicKey));
    
    genInitCtxt_Parallel(Actxts[0], Bctxts[0], Actxt, Bctxt, Initpoly);
    
    //Fixme: HEmatpar.dim or HEmatpar.subdim
    long unit = HEmatpar.subdim * HEmatpar.nbatching;
    NTL_EXEC_RANGE(HEmatpar.subdim - 1, first, last);
    for(long i = first; i < last; ++i){
        long i1 = (i + 1);
        shiftBycols_Parallel(Actxts[i1], Actxts[0], i1, shiftpoly);
        Bctxts[i1] = Bctxts[0];
        rotate(Bctxts[i1], -(unit * (i1)));
    }
    NTL_EXEC_RANGE_END;
    
    HEmatmul_Hadamard(res, Actxts, Bctxts, HEmatpar.subdim);
    
    long index = (long) (log2(HEmatpar.dim/HEmatpar.subdim));
    for(long i = 0; i < index; ++i){
        Ctxt ctemp(res);
        //The difference is here, multiply batching
        rotate(ctemp, -(HEmatpar.dim * HEmatpar.subdim * (1<<i) * HEmatpar.nbatching));
        res += ctemp;
    }
}

void HEmatrix::genMultBPoly(zzX*& Initpoly){
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){ btmp = false; }    //! all the terms have the same numbers
    else{ btmp = true; }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
    
    Initpoly =  new zzX[HEmatpar.dim];
    vector<vector<long>> bvals(HEmatpar.dim, vector<long>(HEmatpar.neededslots));
    
    NTL_EXEC_RANGE(ibound, first, last);
    for(long i = first; i < last; ++i){
        long jbound = HEmatpar.sqrdim;
        if ((btmp)&&(i == ibound - 1)){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        for(long j = 0; j < jbound; ++j){
            long k = i * HEmatpar.sqrdim + j;
            
            for(long l = 0; l < HEmatpar.dim; ++l){
                bvals[k][l * HEmatpar.dim + k] = 1;
            }
            msgrightRotateAndEqual(bvals[k], HEmatpar.neededslots, i * HEmatpar.sqrdim * HEmatpar.dim);
            meta.data->ea.encode(Initpoly[k], bvals[k]);
        }
    }
    NTL_EXEC_RANGE_END;
}

//New function
//void HEmatrix::genMultBPoly_Parallel(zzX*& Initpoly){
//    bool btmp;
//    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){ btmp = false; }    //! all the terms have the same numbers
//    else{ btmp = true; }
//    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim); //! number of "i"
//
//    Initpoly =  new zzX[HEmatpar.dim];
//    vector<vector<long>> bvals(HEmatpar.dim, vector<long>(HEmatpar.neededslots));
//
//    NTL_EXEC_RANGE(ibound, first, last);
//    for(long i = first; i < last; ++i){
//        long jbound = HEmatpar.sqrdim;
//        if ((btmp)&&(i == ibound - 1)){
//            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
//        }
//        for(long j = 0; j < jbound; ++j){
//            long k = i * HEmatpar.sqrdim + j;
//            long start = k * HEmatpar.dim * HEmatpar.nbatching;
//            long end = ( k * HEmatpar.dim + HEmatpar.dim - k) * HEmatpar.nbatching;
//            for(long l = start; l < end; ++l){
//                bvals[k][l] = 1;
//            }
//            msgrightRotateAndEqual(bvals[k], HEmatpar.neededslots, i * HEmatpar.sqrdim * HEmatpar.dim * HEmatpar.nbatching);
//            meta.data->ea.encode(Initpoly[k], bvals[k]);
//        }
//    }
//    NTL_EXEC_RANGE_END;
//}

void HEmatrix::genInitActxt(vector<Ctxt>& Actxts, Mat<long>& mat){
//    Mat<ZZ>* Amat = new Mat<ZZ>[HEmatpar.dim];
    Mat<long>* Amat = new Mat<long>[HEmatpar.dim];
    //Fixme: Is this definition correct?
    Actxts= vector<Ctxt>(HEmatpar.dim, Ctxt(meta.data->publicKey));
    vector<vector<long>> cmsg(HEmatpar.dim, vector<long>(HEmatpar.neededslots));

    NTL_EXEC_RANGE(HEmatpar.dim, first, last);
    for(long k = first; k < last; ++k){
        Amat[k].SetDims(HEmatpar.dim, HEmatpar.dim);
        long dimk = HEmatpar.dim - k;
    
        //! 0 <= i < d - k: shift by (k+i)-positions from mat[i]
        for(long i = 0; i < dimk; ++i){
            long nshift = k + i;
            long nshift2 = HEmatpar.dim - nshift;
            
            for(long j = 0; j < nshift2; ++j){
                Amat[k][i][j] = mat[i][j + nshift];
            }
            
            for(long j = nshift2; j < HEmatpar.dim; ++j){
                Amat[k][i][j] = mat[i][j - nshift2];
            }
        }
        
        //! i = d - k : Amat[k][i] <- mat[i]
        if(k!=0){
            for(long j = 0; j < HEmatpar.dim; ++j){
                Amat[k][dimk][j] = mat[dimk][j];
            }
        }
        
        //! d - k + 1 <= i < d: shift by (k+i-d)-positions from mat[i]
        for(long i = dimk + 1; i < HEmatpar.dim; ++i){
            long nshift =  i - dimk;
            long nshift2 = HEmatpar.dim - nshift;
            
            for(long j = 0; j < nshift2; ++j){
                Amat[k][i][j] = mat[i][j + nshift];
            }
            for(long j = nshift2; j < HEmatpar.dim; ++j){
                Amat[k][i][j] = mat[i][j - nshift2];
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! encryption of d Rmat
    NTL_EXEC_RANGE(HEmatpar.dim, first, last);
    for(long k = first; k < last; ++k){
        for(long i = 0; i < HEmatpar.nrows; ++i){
            for(long j = 0; j < HEmatpar.ncols; ++j){
                long dtemp;
                conv(dtemp, Amat[k][i][j]);
                cmsg[k][i * HEmatpar.dim + j] = dtemp;
            }
        }
        meta.data->ea.encrypt(Actxts[k], meta.data->publicKey, cmsg[k]);
    }
    NTL_EXEC_RANGE_END;
}

void HEmatrix::genInitBctxt(Ctxt& resB, Ctxt& Bctxt, zzX*& poly){
    bool btmp;
    if((HEmatpar.dim % HEmatpar.sqrdim) == 0){
        btmp = false;
    }
    else{
        btmp = true;
    }
    long ibound = (long) ceil((double)HEmatpar.dim/HEmatpar.sqrdim);
    
    
    vector<vector<Ctxt>> Bctemp(HEmatpar.sqrdim, vector<Ctxt>(HEmatpar.sqrdim, Ctxt(meta.data->publicKey)));
    
    //! 0. Store some ciphertexts (0,1,...,d-1), ( ,d+1,...2d-1)
    vector<Ctxt> BaByctxtB(HEmatpar.dim, Ctxt(meta.data->publicKey));
    
    BaByctxtB[0] = Bctxt;
    
    //! i = 0:   Actxts[0] = v[0] + p1 * v[1] + ... + p[sqr(d)-1] *  v[sqr(d)-1]
    resB = BaByctxtB[0];
    resB.multByConstant(poly[0]);

    NTL_EXEC_RANGE(HEmatpar.sqrdim - 1, first, last);
    for(long j = first; j < last; ++j){
        long j1 = (j + 1);
        BaByctxtB[j1] = Bctxt;
        rotate(BaByctxtB[j1], -((j1) * HEmatpar.dim));
        Bctemp[0][j1] = BaByctxtB[j1];
        Bctemp[0][j1].multByConstant(poly[j1]);
    }
    NTL_EXEC_RANGE_END;
    
    for(long j = 1; j < HEmatpar.sqrdim; ++j){
        resB += Bctemp[0][j];
    }
    
    vector<Ctxt> Bctxts(HEmatpar.dim, Ctxt(meta.data->publicKey));
    
    NTL_EXEC_RANGE(HEmatpar.dim - HEmatpar.sqrdim, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k + HEmatpar.sqrdim);
        long i = (long)(k1 / HEmatpar.sqrdim);
        long j = (long)(k1 % HEmatpar.sqrdim);
        
        Bctxts[k1] = BaByctxtB[j];
        Bctxts[k1].multByConstant(poly[k1]);
    }
    NTL_EXEC_RANGE_END;
    
    NTL_EXEC_RANGE(ibound - 1, first, last);
    for(long k = first; k < last; ++k){
        long k1 = (k+1) * HEmatpar.sqrdim;
        long jbound = HEmatpar.sqrdim;
        if((btmp) && (k == ibound - 2)){
            jbound = (HEmatpar.dim % HEmatpar.sqrdim);
        }
        
        for(long j = 1; j < jbound; ++j){
            Bctxts[k1] += Bctxts[k1+j];
        }
        
        long k2 = (k1 - (k1 % HEmatpar.sqrdim));
        rotate(Bctxts[k1], -(k2 * HEmatpar.dim));
    }
    NTL_EXEC_RANGE_END;
    
    for(long k = 1; k < ibound; ++k){
        long k1 = k * HEmatpar.sqrdim;
        resB += Bctxts[k1];
    }
}

//Fixme: consider matrix<long>
void HEmatrix::genInitRecActxt(vector<Ctxt>& Actxts, Mat<long>& mat){
    // rep_mat: (mat; mat; ... ;mat) square mat
    Mat<long> replicate_mat;
    replicate_mat.SetDims(HEmatpar.dim, HEmatpar.dim);
    
    long index_rows = HEmatpar.dim/HEmatpar.subdim;
    
    NTL_EXEC_RANGE(index_rows, first, last);
    for(long k = first; k < last; ++k){
        for(long i = 0; i < HEmatpar.subdim; ++i){
            for(long j = 0; j < HEmatpar.dim; ++j){
                replicate_mat[k*HEmatpar.subdim + i][j] = mat[i][j];
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! generate the (linear transformed) matrices
    Mat<long>* Amat = new Mat<long>[HEmatpar.subdim];
    //Fixme: correct?
    Actxts = vector<Ctxt>(HEmatpar.subdim, Ctxt(meta.data->publicKey));
    
    NTL_EXEC_RANGE(HEmatpar.subdim, first, last);
    for(long k = first; k < last; ++k){
        Amat[k].SetDims(HEmatpar.dim, HEmatpar.dim);
        long dimk= HEmatpar.dim - k;
        
        //! 0 <= i < d - k: shift by (k+i)-positions from mat[i]
        for(long i = 0; i < dimk; ++i){
            long nshift = k + i;
            long nshift2 = HEmatpar.dim - nshift;
            
            for(long j = 0; j < nshift2; ++j){
                Amat[k][i][j] = replicate_mat[i][j + nshift];
            }
            
            for(long j = nshift2; j < HEmatpar.dim; ++j){
                Amat[k][i][j] = replicate_mat[i][j - nshift2];
            }
        }
        
        //! i = d - k : Amat[k][i] <- mat[i]
        if(k!=0){
            for(long j = 0; j < HEmatpar.dim; ++j){
                Amat[k][dimk][j] = replicate_mat[dimk][j];
            }
        }
        
        //! d - k + 1 <= i < d: shift by (k+i-d)-positions from mat[i]
        for(long i = dimk + 1; i < HEmatpar.dim; ++i){
            long nshift =  i - dimk;
            long nshift2 = HEmatpar.dim - nshift;
            
            for(long j = 0; j < nshift2; ++j){
                Amat[k][i][j] = replicate_mat[i][j + nshift];
            }
            for(long j = nshift2; j < HEmatpar.dim; ++j){
                Amat[k][i][j] = replicate_mat[i][j - nshift2];
            }
        }
    }
    NTL_EXEC_RANGE_END;
    
    //! encryption of d Rmat
    vector<vector<long>> cmsg(HEmatpar.subdim, vector<long>(HEmatpar.neededslots));
    NTL_EXEC_RANGE(HEmatpar.subdim, first, last);
    for(long k = first; k < last; ++k){
        for(int i = 0; i < HEmatpar.nrows; ++i){
            for(long j = 0; j < HEmatpar.ncols; ++j){
//                long dtemp;
//                conv(dtemp, Amat[k][i][j]);
                cmsg[k][i * HEmatpar.dim + j] = Amat[k][i][j];
            }
        }
        meta.data->ea.encrypt(Actxts[k], meta.data->publicKey, cmsg[k]);
    }
    NTL_EXEC_RANGE_END;
}

void HEmatrix::HEmatmul_preprocessing(Ctxt& res, vector<Ctxt>& Actxts, Ctxt& Bctxt, zzX*& Initpoly){

    //! 1. Generate the initial ciphertexts
    vector<Ctxt> Bctxts(HEmatpar.dim, Ctxt(meta.data->publicKey));;
    genInitBctxt(Bctxts[0], Bctxt, Initpoly);
    
    //! 2. Row shifting of Bctxt[0]
    NTL_EXEC_RANGE(HEmatpar.dim1, first, last);
    for(long i = first; i < last; ++i){
        long i1 = (i + 1);
        Bctxts[i1] = Bctxts[0];
        rotate(Bctxts[i1], -(HEmatpar.dim * (i1)));
    }
    NTL_EXEC_RANGE_END;
    
    //! 3. Hadamard mult : Actxts[0] * Bctxts[0] + ... + Actxts[d-1] * Bctxts[d-1]
    NTL_EXEC_RANGE(HEmatpar.dim, first, last);
    for(long i = first; i < last; ++i){
        Actxts[i].multiplyBy(Bctxts[i]);
    }
    NTL_EXEC_RANGE_END;
    
    res = Actxts[0];
    for(int i = 1; i < HEmatpar.dim; ++i){
        res += Actxts[i];
    }
}

void HEmatrix::HErmatmul_preprocessing(Ctxt& res, vector<Ctxt>& Actxts, Ctxt& Bctxt, zzX*& Initpoly){
    
    //! 1. Generate the initial ciphertexts
    vector<Ctxt> Bctxts(HEmatpar.subdim, Ctxt(meta.data->publicKey));
    genInitBctxt(Bctxts[0], Bctxt, Initpoly);
    
    //! 2. Row shifting of Bctxt[0]
    NTL_EXEC_RANGE(HEmatpar.subdim - 1, first, last);
    //Fixme: int i or long i
    for(long i = first; i < last; ++i){
        long i1 = (i + 1);
        Bctxts[i1] = Bctxts[0];
        rotate(Bctxts[i1], -(HEmatpar.dim * (i1)));
    }
    NTL_EXEC_RANGE_END;
    
    //! 3. Hadamard mult : Actxts[0] * Bctxts[0] + ... + Actxts[d-1] * Bctxts[d-1]
    NTL_EXEC_RANGE(HEmatpar.subdim, first, last);
    for(long i = first; i < last; ++i){
        Actxts[i].multiplyBy(Bctxts[i]);
    }
    NTL_EXEC_RANGE_END;
    
    //! aggregate the results
    res = Actxts[0];
    for(int i = 1; i < HEmatpar.subdim; ++i){
        res += Actxts[i];
    }
    
    //! 4. shift and aggregate the results
    long index = (long) log2(HEmatpar.dim/HEmatpar.subdim);
    for(long i = 0; i < index; ++i){
        Ctxt ctemp = res;
        rotate(ctemp, -(HEmatpar.dim * HEmatpar.subdim * (1<<i)));
        res += ctemp;
    }
}
