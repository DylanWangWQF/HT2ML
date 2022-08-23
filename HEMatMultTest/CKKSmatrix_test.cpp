//
//  CKKSmatrix_test.cpp
//  ePPDSC
//
//  Created by Qifan Wang on 2/01/22.
//

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <sys/time.h>
#include <fstream>

#include <NTL/RR.h>
#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include "NTL/RR.h"
#include "NTL/vec_RR.h"
#include "NTL/mat_RR.h"
#include <NTL/BasicThreadPool.h>
#include <NTL/mat_ZZ.h>

#include "../src/helib/helib.h"
#include "CKKSmatrix.h"
#include "CKKSmatrix_test.h"

using namespace std;
using namespace helib;
using namespace NTL;
using namespace chrono;

void CKKSmatrixTest::testMult(long nrows){
    ckksParams param(/*m=*/16 * 1024, /*bits=*/150, /*precision=*/20, /*c=*/2);
    ckksMatpar ckksmatpar;
    long ncols = nrows;
    readckksMatpar(ckksmatpar, nrows, ncols);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, po2dim, neededslots) = (" << nrows << "," << ckksmatpar.dim  << "," << ckksmatpar.sqrdim << "," << ckksmatpar.neededslots << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    HELIB_NTIMER_START(ckksmatrix_test_setup);
    ckksMeta meta;
    meta(param);
    cout << "Context contents: " << endl;
    meta.data->context.printout();
    CKKSmatrix CKKSmatrix(ckksmatpar, meta);
    HELIB_NTIMER_STOP(ckksmatrix_test_setup);
    printNamedTimer(cout, "ckksmatrix_test_setup");
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
    Mat<RR> Amat;
    Mat<RR> Bmat;
    
    Amat.SetDims(nrows, ncols);
    Bmat.SetDims(nrows, ncols);
    
    for(long i = 0; i < nrows ; i++){
        for(long j = 0; j < ncols; j++){
            Amat[i][j]= to_RR((((i * 2 + ncols * j) % 3)/10.0));
            Bmat[i][j]= to_RR((((i * ncols + j ) %3 )/10.0));
        }
    }
//    cout << "Generated matrix Amat: " << endl << Amat << endl;
//    cout << "Generated matrix Bmat: " << endl << Bmat << endl;
    
    /*---------------------------------------*/
    //  Encryption
    /*---------------------------------------*/
    HELIB_NTIMER_START(enc_ckksmatrix);
    Ctxt Actxt(meta.data->publicKey);
    CKKSmatrix.encryptRmat(Actxt, Amat);
    Ctxt Bctxt(meta.data->publicKey);
    CKKSmatrix.encryptRmat(Bctxt, Bmat);
    HELIB_NTIMER_STOP(enc_ckksmatrix);
    printNamedTimer(cout, "enc_ckksmatrix");
    
    /*---------------------------------------*/
    //  GenPoly
    /*---------------------------------------*/
    vector<vector<EncodedPtxt>> Initpoly;
    HELIB_NTIMER_START(GenMulPoly_testMul);
    CKKSmatrix.genMultPoly(Initpoly);
    HELIB_NTIMER_STOP(GenMulPoly_testMul);
    printNamedTimer(cout, "GenMulPoly_testMul");
    
    vector<EncodedPtxt> shiftpoly;
    HELIB_NTIMER_START(GenShiftPoly_testMul);
    CKKSmatrix.genShiftPoly(shiftpoly);
    HELIB_NTIMER_STOP(GenShiftPoly_testMul);
    printNamedTimer(cout, "GenShiftPoly_testMul");
    
    /*---------------------------------------*/
    //  Mult
    /*---------------------------------------*/
    
    HELIB_NTIMER_START(Mul_testMul);
    Ctxt Cctxt(meta.data->publicKey);
    cout << "Cctxt.capacity=" << Cctxt.capacity() << " ";
    cout << "Cctxt.errorBound=" << Cctxt.errorBound() << "\n";
    CKKSmatrix.HEmatmul(Cctxt, Actxt, Bctxt, Initpoly, shiftpoly);
    cout << "Cctxt.capacity=" << Cctxt.capacity() << " ";
    cout << "Cctxt.errorBound=" << Cctxt.errorBound() << "\n";
    HELIB_NTIMER_STOP(Mul_testMul);
    printNamedTimer(cout, "Mul_testMul");
    
    /*---------------------------------------*/
    //  Decryption
    /*---------------------------------------*/
    Mat<RR> CKKSresmat;
    CKKSmatrix.decryptRmat(CKKSresmat, Cctxt);
    
    /*---------------------------------------*/
    //  Error
    /*---------------------------------------*/
    Mat<RR> resmat;
    mul(resmat, Amat, Bmat);
    
//    RR error = getError(resmat, CKKSresmat, nrows, ncols);
//    cout << "Error (mul): " << error << endl;
//
//    cout << "------------------" << endl;
//    cout << "Plaintext" << endl;
//    printRmatrix(resmat, nrows);
//    cout << "------------------" << endl;
//    cout << "Encryption" << endl;
//    printRmatrix(CKKSresmat, nrows);
//    cout << "------------------" << endl;
}

void CKKSmatrixTest::testRMult(long nrows, long subdim){
    ckksParams param(/*m=*/1024, /*bits=*/358, /*precision=*/20, /*c=*/2);
    ckksMatpar ckksmatpar;
    long ncols = nrows;
    readckksMatpar(ckksmatpar, nrows, ncols, subdim);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, po2dim, neededslots) = (" << nrows << "," << ckksmatpar.dim  << "," << ckksmatpar.sqrdim << "," << ckksmatpar.neededslots << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    HELIB_NTIMER_START(ckksmatrix_test_setup);
    ckksMeta meta;
    meta(param);
    cout << "Context contents: " << endl;
    meta.data->context.printout();
    CKKSmatrix CKKSmatrix(ckksmatpar, meta);
    HELIB_NTIMER_STOP(ckksmatrix_test_setup);
    printNamedTimer(cout, "ckksmatrix_test_setup");
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
    Mat<RR> rAmat;
    rAmat.SetDims(subdim, ncols);
    
    Mat<RR> Bmat;
    Bmat.SetDims(nrows, ncols);
    
    for(long i = 0; i < subdim; i++){
        for(long j = 0; j < ncols; j++){
            rAmat[i][j]= to_RR((((i + ncols * j)%3)/10.0));
        }
    }
    
    for(long i = 0; i < nrows ; i++){
        for(long j = 0; j < ncols; j++){
            Bmat[i][j]= to_RR((((i*ncols + j )%3)/10.0));
        }
    }
    
    // replicate
    Mat<RR> rAmat1;
    rAmat1.SetDims(nrows, ncols);
    for(long i = 0; i < nrows/subdim; ++i){
        for(long j = 0; j < subdim; ++j){
            for(long k = 0; k< ncols; k++){
                rAmat1[i*subdim + j][k] = rAmat[j][k];
            }
        }
    }
    cout << "Generated matrix rAmat1: " << endl << rAmat1 << endl;
    cout << "Generated matrix Bmat: " << endl << Bmat << endl;
    
    /*---------------------------------------*/
    //  Encryption
    /*---------------------------------------*/
    HELIB_NTIMER_START(enc_ckksmatrix);
    Ctxt Actxt(meta.data->publicKey);
    CKKSmatrix.encryptRmat(Actxt, rAmat1);
    Ctxt Bctxt(meta.data->publicKey);
    CKKSmatrix.encryptRmat(Bctxt, Bmat);
    HELIB_NTIMER_STOP(enc_ckksmatrix);
    printNamedTimer(cout, "enc_ckksmatrix");
    
    /*---------------------------------------*/
    //  GenPoly
    /*---------------------------------------*/
    vector<vector<EncodedPtxt>> Initpoly;
    HELIB_NTIMER_START(GenMulPoly_testMul);
    CKKSmatrix.genMultPoly(Initpoly);
    HELIB_NTIMER_STOP(GenMulPoly_testMul);
    printNamedTimer(cout, "GenMulPoly_testMul");
    
    vector<EncodedPtxt> shiftpoly;
    HELIB_NTIMER_START(GenShiftPoly_testMul);
    CKKSmatrix.genShiftPoly(shiftpoly);
    HELIB_NTIMER_STOP(GenShiftPoly_testMul);
    printNamedTimer(cout, "GenShiftPoly_testMul");
    
    /*---------------------------------------*/
    //  Mult
    /*---------------------------------------*/
    HELIB_NTIMER_START(Mul_testMul);
    Ctxt Cctxt(meta.data->publicKey);
    CKKSmatrix.HErmatmul(Cctxt, Actxt, Bctxt, Initpoly, shiftpoly);
    HELIB_NTIMER_STOP(Mul_testMul);
    printNamedTimer(cout, "Mul_testMul");
    
    /*---------------------------------------*/
    //  Decryption
    /*---------------------------------------*/
    Mat<RR> CKKSresmat;
    CKKSmatrix.decryptRmat(CKKSresmat, Cctxt);
    
    /*---------------------------------------*/
    //  Error
    /*---------------------------------------*/
    Mat<RR> resmat;
    mul(resmat, rAmat1, Bmat);
    
    RR error = getError(resmat, CKKSresmat, subdim, ncols);
    cout << "Error (mul): " << error << endl;
    
    cout << "------------------" << endl;
    cout << "Plaintext" << endl;
    printRmatrix(resmat, nrows);
    cout << "------------------" << endl;
    cout << "Encryption" << endl;
    printRmatrix(CKKSresmat, nrows);
    cout << "------------------" << endl;
}

void CKKSmatrixTest::testMult_preprocessing(long nrows){
    ckksParams param(/*m=*/1024, /*bits=*/358, /*precision=*/20, /*c=*/2);
    ckksMatpar ckksmatpar;
    long ncols = nrows;
    readckksMatpar(ckksmatpar, nrows, ncols);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, po2dim, neededslots) = (" << nrows << "," << ckksmatpar.dim  << "," << ckksmatpar.sqrdim << "," << ckksmatpar.neededslots << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    HELIB_NTIMER_START(ckksmatrix_test_setup);
    ckksMeta meta;
    meta(param);
    cout << "Context contents: " << endl;
    meta.data->context.printout();
    CKKSmatrix CKKSmatrix(ckksmatpar, meta);
    HELIB_NTIMER_STOP(ckksmatrix_test_setup);
    printNamedTimer(cout, "ckksmatrix_test_setup");
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
    Mat<RR> Amat;
    Mat<RR> Bmat;
    
    Amat.SetDims(nrows, ncols);
    Bmat.SetDims(nrows, ncols);
    
    for(long i = 0; i < nrows ; i++){
        for(long j = 0; j < ncols; j++){
            Amat[i][j]= to_RR((((i * 2 + ncols * j) % 3)/10.0));
            Bmat[i][j]= to_RR((((i * ncols + j ) %3 )/10.0));
        }
    }
    cout << "Generated matrix Amat: " << endl << Amat << endl;
    cout << "Generated matrix Bmat: " << endl << Bmat << endl;
    
    /*---------------------------------------*/
    //  Encryption
    /*---------------------------------------*/
    HELIB_NTIMER_START(enc_ckksmatrix);
    vector<Ctxt> Actxts;
    CKKSmatrix.genInitActxt(Actxts, Amat);
    Ctxt Bctxt(meta.data->publicKey);
    CKKSmatrix.encryptRmat(Bctxt, Bmat);
    HELIB_NTIMER_STOP(enc_ckksmatrix);
    printNamedTimer(cout, "enc_ckksmatrix");
    
    /*---------------------------------------*/
    //  GenPoly
    /*---------------------------------------*/
//    vector<vector<EncodedPtxt>> Initpoly;
//    HELIB_NTIMER_START(GenMulPoly_testMul);
//    CKKSmatrix.genMultPoly(Initpoly);
//    HELIB_NTIMER_STOP(GenMulPoly_testMul);
//    printNamedTimer(cout, "GenMulPoly_testMul");
//
//    vector<EncodedPtxt> shiftpoly;
//    HELIB_NTIMER_START(GenShiftPoly_testMul);
//    CKKSmatrix.genShiftPoly(shiftpoly);
//    HELIB_NTIMER_STOP(GenShiftPoly_testMul);
//    printNamedTimer(cout, "GenShiftPoly_testMul");
    vector<EncodedPtxt> Initpoly;
    CKKSmatrix.genMultBPoly(Initpoly);
    
    /*---------------------------------------*/
    //  Mult
    /*---------------------------------------*/
    HELIB_NTIMER_START(Mul_testMul);
    Ctxt Cctxt(meta.data->publicKey);
    CKKSmatrix.HEmatmul_preprocessing(Cctxt, Actxts, Bctxt, Initpoly);
    HELIB_NTIMER_STOP(Mul_testMul);
    printNamedTimer(cout, "Mul_testMul");
    
    /*---------------------------------------*/
    //  Decryption
    /*---------------------------------------*/
    Mat<RR> CKKSresmat;
    CKKSmatrix.decryptRmat(CKKSresmat, Cctxt);
    
    /*---------------------------------------*/
    //  Error
    /*---------------------------------------*/
    Mat<RR> resmat;
    mul(resmat, Amat, Bmat);
    
    RR error = getError(resmat, CKKSresmat, nrows, ncols);
    cout << "Error (mul): " << error << endl;
    
    cout << "------------------" << endl;
    cout << "Plaintext" << endl;
    printRmatrix(resmat, nrows);
    cout << "------------------" << endl;
    cout << "Encryption" << endl;
    printRmatrix(CKKSresmat, nrows);
    cout << "------------------" << endl;
}

void CKKSmatrixTest::testRMult_preprocessing(long nrows, long subdim){
    ckksParams param(/*m=*/1024, /*bits=*/358, /*precision=*/20, /*c=*/2);
    ckksMatpar ckksmatpar;
    long ncols = nrows;
    readckksMatpar(ckksmatpar, nrows, ncols, subdim);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, po2dim, neededslots) = (" << nrows << "," << ckksmatpar.dim  << "," << ckksmatpar.sqrdim << "," << ckksmatpar.neededslots << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    HELIB_NTIMER_START(ckksmatrix_test_setup);
    ckksMeta meta;
    meta(param);
    cout << "Context contents: " << endl;
    meta.data->context.printout();
    CKKSmatrix CKKSmatrix(ckksmatpar, meta);
    HELIB_NTIMER_STOP(ckksmatrix_test_setup);
    printNamedTimer(cout, "ckksmatrix_test_setup");
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
    Mat<RR> rAmat;
    rAmat.SetDims(subdim, ncols);
    
    Mat<RR> Bmat;
    Bmat.SetDims(nrows, ncols);
    
    for(long i = 0; i < subdim; i++){
        for(long j = 0; j < ncols; j++){
            rAmat[i][j]= to_RR((((i + ncols * j)%3)/10.0));
        }
    }
    
    for(long i = 0; i < nrows ; i++){
        for(long j = 0; j < ncols; j++){
            Bmat[i][j]= to_RR((((i*ncols + j )%3)/10.0));
        }
    }
    
    // replicate
    Mat<RR> rAmat1;
    rAmat1.SetDims(nrows, ncols);
    for(long i = 0; i < nrows/subdim; ++i){
        for(long j = 0; j < subdim; ++j){
            for(long k = 0; k< ncols; k++){
                rAmat1[i*subdim + j][k] = rAmat[j][k];
            }
        }
    }
    cout << "Generated matrix rAmat1: " << endl << rAmat1 << endl;
    cout << "Generated matrix Bmat: " << endl << Bmat << endl;
    
    /*---------------------------------------*/
    //  Encryption
    /*---------------------------------------*/
    HELIB_NTIMER_START(enc_ckksmatrix);
    vector<Ctxt> Actxts;
    CKKSmatrix.genInitRecActxt(Actxts, rAmat);
    Ctxt Bctxt(meta.data->publicKey);
    CKKSmatrix.encryptRmat(Bctxt, Bmat);
    HELIB_NTIMER_STOP(enc_ckksmatrix);
    printNamedTimer(cout, "enc_ckksmatrix");
    
    /*---------------------------------------*/
    //  GenPoly
    /*---------------------------------------*/
//    vector<vector<EncodedPtxt>> Initpoly;
//    HELIB_NTIMER_START(GenMulPoly_testMul);
//    CKKSmatrix.genMultPoly(Initpoly);
//    HELIB_NTIMER_STOP(GenMulPoly_testMul);
//    printNamedTimer(cout, "GenMulPoly_testMul");
//
//    vector<EncodedPtxt> shiftpoly;
//    HELIB_NTIMER_START(GenShiftPoly_testMul);
//    CKKSmatrix.genShiftPoly(shiftpoly);
//    HELIB_NTIMER_STOP(GenShiftPoly_testMul);
//    printNamedTimer(cout, "GenShiftPoly_testMul");
    vector<EncodedPtxt> Initpoly;
    CKKSmatrix.genMultBPoly(Initpoly);
    
    /*---------------------------------------*/
    //  Mult
    /*---------------------------------------*/
    HELIB_NTIMER_START(Mul_testMul);
    Ctxt Cctxt(meta.data->publicKey);
    CKKSmatrix.HErmatmul_preprocessing(Cctxt, Actxts, Bctxt, Initpoly);
    HELIB_NTIMER_STOP(Mul_testMul);
    printNamedTimer(cout, "Mul_testMul");
    
    /*---------------------------------------*/
    //  Decryption
    /*---------------------------------------*/
    Mat<RR> CKKSresmat;
    CKKSmatrix.decryptRmat(CKKSresmat, Cctxt);
    
    /*---------------------------------------*/
    //  Error
    /*---------------------------------------*/
    Mat<RR> resmat;
    mul(resmat, rAmat, Bmat);
    
    RR error = getError(resmat, CKKSresmat, subdim, ncols);
    cout << "Error (mul): " << error << endl;
    
    cout << "------------------" << endl;
//    cout << "Plaintext" << endl;
    printRmatrix(resmat, subdim);
    cout << "------------------" << endl;
    cout << "Encryption" << endl;
    printRmatrix(CKKSresmat, nrows);
    cout << "------------------" << endl;
}
