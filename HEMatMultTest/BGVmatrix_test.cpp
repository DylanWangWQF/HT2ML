#include <cmath>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <sys/time.h>
#include <fstream>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/matrix.h>
#include <helib/helib.h>

#include "BGVmatrix.h"
#include "BGVmatrix_test.h"

using namespace std;
using namespace helib;
using namespace chrono;

struct testS {
    SecKey secretKey;
    const PubKey& publicKey;
    const EncryptedArray& ea;

    testS(const Context& context) :
          secretKey(context),
          publicKey((secretKey.GenSecKey(),
                     addSome1DMatrices(secretKey),
                     secretKey)),
          ea(context.getEA()){}
};

void HEmatrixTest::FindCorrectM(){
//    cout << "Finding m starts: " << endl;
//    int p, j;
//    for(p = 300; p <= 1000; p++)
//    {
//        for(j = 2 ; j < p; j++)
//        {
//            if(p % j == 0) break;
//        }
//        if(j >= p){
//            cout << "Current prime is: " << p << endl;
//            long m = FindM(/*k=*/80, /*nBits=*/500, /*c=*/2, /*p=*/p, /*d=*/1, /*s=*/240, /*chosen_m=*/0, /*verbose=*/true);
//            cout << "Selected m is: " << m << endl;
//            // p=2; ph(m)=4096, m=4369; k-bit security too low, for bits=200, k=30;
//            Params param(/*m=*/m, /*p=*/p, /*r=*/1, /*bits=*/500, /*c=*/2);
//            Meta meta;
//            meta(param);
//            if (meta.data->ea.size() == 256 || meta.data->ea.size() ==  4096) {
//                cout << "Print context start: " << endl;
//                meta.data->context.printout();
//                cout << "Print context end: " << endl;
//            }
//        }
//    }
    
//    long m = FindM(/*k=*/80, /*nBits=*/200, /*c=*/2, /*p=*/7, /*d=*/1, /*s=*/256, /*chosen_m=*/0, /*verbose=*/true);
//    cout << "Selected m is: " << m << endl;
    
//    HELIB_NTIMER_START(FindSlots);
//    Params param(/*m=*/16384, /*p=*/8191, /*r=*/1, /*bits=*/120, /*c=*/2);
//    Params param(/*m=*/8704, /*p=*/31, /*r=*/1, /*bits=*/119, /*c=*/2); 73 security
//    Params param(/*m=*/10240, /*p=*/127, /*r=*/1, /*bits=*/120, /*c=*/2);
//    Params param(/*m=*/20480, /*p=*/127, /*r=*/1, /*bits=*/200, /*c=*/2);
//    Params param(/*m=*/18944, /*p=*/127, /*r=*/1, /*bits=*/150, /*c=*/2);
//    Params param(/*m=*/20191, /*p=*/47, /*r=*/1, /*bits=*/200, /*c=*/2);
//    Params param(/*m=*/30720, /*p=*/257, /*r=*/1, /*bits=*/180, /*c=*/2);
//    Params param(/*m=*/15360, /*p=*/257, /*r=*/1, /*bits=*/180, /*c=*/2);
//    Params param(/*m=*/17408, /*p=*/127, /*r=*/1, /*bits=*/180, /*c=*/2);
//    Params param(/*m=*/13056, /*p=*/191, /*r=*/1, /*bits=*/180, /*c=*/2);
//    Params param(/*m=*/19968, /*p=*/257, /*r=*/1, /*bits=*/180, /*c=*/2);
//    Params param(/*m=*/30720, /*p=*/5119, /*r=*/1, /*bits=*/180, /*c=*/2);
//    Params param(/*m=*/32768, /*p=*/8191, /*r=*/1, /*bits=*/180, /*c=*/2);
//    Params param(/*m=*/13568, /*p=*/1151, /*r=*/1, /*bits=*/120, /*c=*/2);
//    Params param(/*m=*/9472, /*p=*/127, /*r=*/1, /*bits=*/119, /*c=*/2);
//    Params param(/*m=*/40960, /*p=*/40961, /*r=*/1, /*bits=*/180, /*c=*/2);
//    Meta meta;
//    meta(param);
//    meta.data->context.printout();
//    HELIB_NTIMER_STOP(FindSlots);
//    printNamedTimer(cout, "FindSlots");
    
    //proper params
    //security=82.2439-->/*m=*/10240, /*p=*/127, /*r=*/1, /*bits=*/110, /*c=*/2
    //security=80-->/*m=*/20480, /*p=*/127, /*r=*/1, /*bits=*/200, /*c=*/2
    //security=187.479-->/*m=*/40960, /*p=*/6143, /*r=*/1, /*bits=*/200, /*c=*/2
    //security=118.653-->/*m=*/18944, /*p=*/127, /*r=*/1, /*bits=*/150, /*c=*/2
    
    Params param(/*m=*/19887, /*p=*/127, /*r=*/1, /*bits=*/500, /*c=*/2);
    Meta meta;
    meta(param);
    meta.data->context.printout();
    
//    vector<long> vec1(meta.data->ea.size());
//    for(int i = 0; i < vec1.size() ; i++)
//    {
//        vec1[i] = 1;
//    }
//
//    Ctxt ctxt1(meta.data->publicKey);
//    meta.data->ea.encrypt(ctxt1, meta.data->publicKey, vec1);
//    Ctxt ctxt2(meta.data->publicKey);
//    meta.data->ea.encrypt(ctxt2, meta.data->publicKey, vec1);
//
//    stringstream ss;
//    ctxt1.writeTo(ss);
//    string ctxtr = ss.str();
//    cout << "size of ctxt1 before addition = " << (double) ctxtr.size()/(1024 * 1024) << " MB" << endl;
    
    // ctxt addition for 100 times, check the result ciphertext size
    
}

void HEmatrixTest::testMult(long nrows) {

    long ncols = nrows;
    
//    Params param(/*m=*/12800, /*p=*/127, /*r=*/1, /*bits=*/119, /*c=*/2);
//    Params param(/*m=*/16384, /*p=*/6143, /*r=*/1, /*bits=*/180, /*c=*/2);
//    Params param(/*m=*/32768, /*p=*/8191, /*r=*/1, /*bits=*/180, /*c=*/2);
//    Params param(/*m=*/20480, /*p=*/8191, /*r=*/1, /*bits=*/180, /*c=*/2);
//    Params param(/*m=*/17408, /*p=*/127, /*r=*/1, /*bits=*/180, /*c=*/2);
//    Params param(/*m=*/13056, /*p=*/191, /*r=*/1, /*bits=*/180, /*c=*/2);
//    Params param(/*m=*/12800, /*p=*/127, /*r=*/1, /*bits=*/119, /*c=*/2);
//    Params param(/*m=*/13568, /*p=*/1151, /*r=*/1, /*bits=*/120, /*c=*/2);
//    Params param(/*m=*/20480, /*p=*/127, /*r=*/1, /*bits=*/120, /*c=*/2);
//    Params param(/*m=*/9344, /*p=*/383, /*r=*/1, /*bits=*/120, /*c=*/2);
//    Params param(/*m=*/9472, /*p=*/127, /*r=*/1, /*bits=*/119, /*c=*/2);
//    Params param(/*m=*/16384, /*p=*/8191, /*r=*/1, /*bits=*/120, /*c=*/2);
//    Params param(/*m=*/16384, /*p=*/6143, /*r=*/1, /*bits=*/180, /*c=*/2);
    // Params param(/*m=*/40960, /*p=*/40961, /*r=*/1, /*bits=*/180, /*c=*/2);
    Params param(/*m=*/16384, /*p=*/8191, /*r=*/1, /*bits=*/120, /*c=*/2);
    
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, po2dim, neededslots) = (" << nrows << ", " << HEmatpar.dim  << ", " << HEmatpar.sqrdim << ", " << HEmatpar.neededslots << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    
    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    
    Meta meta;
    meta(param);
    meta.data->context.printout();
    HEmatrix HEmatrix(HEmatpar, meta);
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
    Mat<long> Amat;
    Mat<long> Bmat;
    
    Amat.SetDims(nrows, ncols);
    Bmat.SetDims(nrows, ncols);
    
    for(long i = 0; i < nrows ; i++){
        for(long j = 0; j < ncols; j++){
            Amat[i][j]= i % 2;
            Bmat[i][j]= j % 2;
        }
    }
    // cout << "Generated matrixA: " << Amat << endl;
    // cout << "Generated matrixB: " << Bmat << endl;
    
    /*---------------------------------------*/
    //  Encryption
    /*---------------------------------------*/
    
    Ctxt Actxt(meta.data->publicKey);
    HEmatrix.encryptZmat(Actxt, Amat);
    Ctxt Bctxt(meta.data->publicKey);
    HEmatrix.encryptZmat(Bctxt, Bmat);

    stringstream ss;
    Actxt.writeTo(ss);
    string ClientTemp = ss.str();
    uint64_t client_totalLength = ClientTemp.size();
    
    cout << "Ciphertext size = : " << ((double) client_totalLength / (double)(1024 * 1024)) << " MB" << endl;
    
    /*---------------------------------------*/
    //  GenPoly
    /*---------------------------------------*/
    
    zzX** Initpoly;
    HELIB_NTIMER_START(GenMulPoly_testMul);
    HEmatrix.genMultPoly(Initpoly);
    HELIB_NTIMER_STOP(GenMulPoly_testMul);
    // printNamedTimer(cout, "GenMulPoly_testMul");
    
    zzX* shiftpoly;
    HELIB_NTIMER_START(GenShiftPoly_testMul);
    HEmatrix.genShiftPoly(shiftpoly);
    HELIB_NTIMER_STOP(GenShiftPoly_testMul);
    // printNamedTimer(cout, "GenShiftPoly_testMul");
    
    /*---------------------------------------*/
    //  Mult
    /*---------------------------------------*/
    
    HELIB_NTIMER_START(Mul_testMul);
    Ctxt Cctxt(meta.data->publicKey);
    HEmatrix.HEmatmul(Cctxt, Actxt, Bctxt, Initpoly, shiftpoly);
    HELIB_NTIMER_STOP(Mul_testMul);
    printNamedTimer(cout, "Mul_testMul");
    
    /*---------------------------------------*/
    //  Decryption
    /*---------------------------------------*/
    Mat<long> HEresmat;
    HEmatrix.decryptZmat(HEresmat, Cctxt);
    
    /*---------------------------------------*/
    //  Result Comparison
    /*---------------------------------------*/
//    Mat<long> resmat;
//    mul(resmat, Amat, Bmat);
    
//    cout << "------------------" << endl;
//    cout << "Plaintext Multiplication:" << endl;
//    cout << resmat << endl;
    cout << "------------------" << endl;
    cout << "Encrypted Multiplication:" << endl;
    cout << HEresmat << endl;
    cout << "------------------" << endl;
}

void HEmatrixTest::testRotation(){
    Params param(/*m=*/18944, /*p=*/127, /*r=*/1, /*bits=*/150, /*c=*/2);
    Meta meta;
    meta(param);
    meta.data->context.printout();
    
    vector<long> cmsg1(meta.data->ea.size());
    vector<long> cmsg2(meta.data->ea.size());
    for (int i = 0; i < 16; i++) {
        cmsg1[i] = i+1;
        cmsg2[i] = i+1;
    }
    
    cout << "cmsg1: " << endl;
    cout << cmsg1 << endl;
    
    // Encryption
    HELIB_NTIMER_START(Enc_testSecKeyEnc);
//    Ctxt ctxt1(meta.data->publicKey);
//    meta.data->ea.encrypt(ctxt1, meta.data->publicKey, cmsg1);
    Ctxt ctxt1(meta.data->secretKey);
    EncodedPtxt eptxt;
    meta.data->ea.encode(eptxt, cmsg1);
    meta.data->secretKey.Encrypt(ctxt1, eptxt);
    HELIB_NTIMER_STOP(Enc_testSecKeyEnc);
    printNamedTimer(cout, "Enc_testSecKeyEnc");
    
    vector<long> res1(meta.data->ea.size());
    
    HELIB_NTIMER_START(Dec_testSecKeyEnc);
    meta.data->ea.decrypt(ctxt1, meta.data->secretKey, res1);
    HELIB_NTIMER_STOP(Dec_testSecKeyEnc);
    printNamedTimer(cout, "Dec_testSecKeyEnc");
    
    cout << "res1: " << endl;
    cout << res1 << endl;
    
    //Rotation
//    start= chrono::steady_clock::now();
    
//    long eadim = meta.data->ea.dimension();
//    cout << "eadim: " << eadim << endl;
//    cout << "rotate1: " << endl;
//    rotate(ctxt1, -5);
//    //shift(ctxt1, -32);
//    cout << "rotate2: " << endl;
//    rotate1D(ctxt2, 0, -5, true);
//    rotate1D(ctxt2, 1, -1, true);
//
//    vector<long> res1(meta.data->ea.size());
//    vector<long> res2(meta.data->ea.size());
//    meta.data->ea.decrypt(ctxt1, meta.data->secretKey, res1);
//    meta.data->ea.decrypt(ctxt2, meta.data->secretKey, res2);
//
//    cout << "res1: " << endl;
//    cout << res1 << endl;
//    cout << "res2: " << endl;
//    cout << res2 << endl;
    
//    end = std::chrono::steady_clock::now();
//    diff = end - start;
//    timeElapsed = chrono::duration <double, milli> (diff).count()/1000.0;
//    cout << "Rotation time= " << timeElapsed << " s" << endl;
//    cout << "------------------" << endl;
}

void HEmatrixTest::testEnc(long nrows){
    long ncols = nrows;
    Params param(/*m=*/12800, /*p=*/127, /*r=*/1, /*bits=*/119, /*c=*/2);
    
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, po2dim, neededslots) = (" << nrows << "," << HEmatpar.dim  << "," << HEmatpar.sqrdim << "," << HEmatpar.neededslots << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    HELIB_NTIMER_START(hematrix_test_setup);
    Meta meta;
    meta(param);
    cout << "Context contents: " << endl;
    meta.data->context.printout();
    HEmatrix HEmatrix(HEmatpar, meta);
    HELIB_NTIMER_STOP(hematrix_test_setup);
    printNamedTimer(cout, "hematrix_test_setup");
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
//    Mat<ZZ> Amat;
    Mat<long> Amat;
    Amat.SetDims(nrows, ncols);
        
    for(long i = 0; i < nrows ; i++){
        for(long j = 0; j < ncols; j++){
            Amat[i][j]= (i + 1);
        }
    }
    cout << "Generated matrix: " << Amat << endl;
    
    /*---------------------------------------*/
    //  Encryption
    /*---------------------------------------*/
    HELIB_NTIMER_START(enc_hematrix);
    Ctxt ctxt(meta.data->publicKey);
    HEmatrix.encryptZmat(ctxt, Amat);
    HELIB_NTIMER_STOP(enc_hematrix);
    printNamedTimer(cout, "enc_hematrix");
    stringstream encss;
    ctxt.writeTo(encss);
    string ctxtr = encss.str();
    cout << "Ctxt size:" << ctxtr.size() << endl;
    
    /*---------------------------------------*/
    //  Decryption
    /*---------------------------------------*/
    Ctxt newCtxt(meta.data->publicKey);
    //newCtxt.read(encss);
    newCtxt.Ctxt::read(encss);
//    Mat<ZZ> resmat;
    Mat<long> resmat;
    HELIB_NTIMER_START(dec_hematrix);
    HEmatrix.decryptZmat(resmat, newCtxt);
    HELIB_NTIMER_STOP(dec_hematrix);
    printNamedTimer(cout, "dec_hematrix");
    
//    Ctxt ctxt(meta.data->publicKey);
//    Mat<ZZ> resmat;
//    for (int i = 0; i < 3; i++) {
//        HEmatrix.encryptZmat(ctxt, Amat);
//        HEmatrix.decryptZmat(resmat, ctxt);
//        cout << "Decrypted matrix: " << resmat << endl;
//    }
    
//    cout << "Decrypt correct: " << ctxt.isCorrect() << endl;
//
    cout << "Decrypted matrix: " << resmat << endl;
}

void HEmatrixTest::testAdd(long nrows) {
    long ncols = nrows;
    Params param(/*m=*/19801, /*p=*/4999, /*r=*/1, /*bits=*/200, /*c=*/2);
    
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, po2dim, nslots) = (" << nrows << "," << HEmatpar.dim  << "," << HEmatpar.sqrdim << "," << HEmatpar.neededslots << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    Meta meta;
    meta(param);
    cout << "Context contents: " << endl;
    meta.data->context.printout();
    HEmatrix HEmatrix(HEmatpar, meta);
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
    Mat<long> Amat;
    Mat<long> Bmat;
    Amat.SetDims(nrows, ncols);
    Bmat.SetDims(nrows, ncols);
        
    for(long i = 0; i < nrows ; i++){
        for(long j = 0; j < ncols; j++){
            Amat[i][j]= i + 1;
            Bmat[i][j]= i * 2 + j;
        }
    }
    cout << "Generated matrixA: " << endl;
    cout << Amat << endl;
    cout << "Generated matrixB: " << endl;
    cout << Bmat << endl;
    
    /*---------------------------------------*/
    //  Encryption
    /*---------------------------------------*/
    HELIB_NTIMER_START(enc_testAdd);
    Ctxt Actxt(meta.data->publicKey);
    HEmatrix.encryptZmat(Actxt, Amat);
    Ctxt Bctxt(meta.data->publicKey);
    HEmatrix.encryptZmat(Bctxt, Bmat);
    HELIB_NTIMER_STOP(enc_testAdd);
    printNamedTimer(cout, "enc_testAdd");
    
    /*---------------------------------------*/
    //  Addition
    /*---------------------------------------*/
    HELIB_NTIMER_START(add_testAdd);
    Actxt += Bctxt;
    HELIB_NTIMER_STOP(add_testAdd);
    printNamedTimer(cout, "add_testAdd");
    
    /*---------------------------------------*/
    //  Decryption
    /*---------------------------------------*/
    Mat<long> resmat;
    HELIB_NTIMER_START(dec_testAdd);
    HEmatrix.decryptZmat(resmat, Actxt);
    HELIB_NTIMER_STOP(dec_testAdd);
    printNamedTimer(cout, "dec_testAdd");
    
//    Mat<long> plain_mat;
//    add(plain_mat, Amat, Bmat);
//    cout << "Plaintext matrix after addition: " << plain_mat << endl;
//    cout << "Decrypted matrix after addition: " << resmat << endl;
}

void HEmatrixTest::testTrans(long nrows) {
    long ncols = nrows;
    //Params param(/*m=*/20480, /*p=*/127, /*r=*/1, /*bits=*/200, /*c=*/2);
    Params param(/*m=*/20191, /*p=*/47, /*r=*/1, /*bits=*/200, /*c=*/2);
    
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, po2dim, neededslots) = (" << nrows << ", " << HEmatpar.dim  << ", " << HEmatpar.sqrdim << ", " << HEmatpar.neededslots << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    
    Meta meta;
    meta(param);
    cout << "Context contents: " << endl;
    meta.data->context.printout();
    HEmatrix HEmatrix(HEmatpar, meta);
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
    Mat<long> Amat;
    Amat.SetDims(nrows, ncols);
        
    for(int i = 0; i < nrows ; i++){
        for(int j = 0; j < ncols; j++){
            Amat[i][j]= i + 1;
        }
    }
    cout << "Generated matrix: " << Amat << endl;
    
    /*---------------------------------------*/
    //  Encryption
    /*---------------------------------------*/
    HELIB_NTIMER_START(enc_testTran);
    Ctxt ctxt(meta.data->publicKey);
    HEmatrix.encryptZmat(ctxt, Amat);
    HELIB_NTIMER_STOP(enc_testTran);
    printNamedTimer(cout, "enc_testTran");
    
    /*---------------------------------------*/
    //  Transposition
    /*---------------------------------------*/
    zzX* transpoly;
    HELIB_NTIMER_START(genTran_testTran);
    HEmatrix.genTransPoly(transpoly);
    HELIB_NTIMER_STOP(genTran_testTran);
    printNamedTimer(cout, "genTran_testTran");
    
    Ctxt Tctxt(meta.data->publicKey);
    
    HELIB_NTIMER_START(tran_testTran);
    HEmatrix.transpose(Tctxt, ctxt, transpoly);
    HELIB_NTIMER_STOP(tran_testTran);
    printNamedTimer(cout, "tran_testTran");
    
    /*---------------------------------------*/
    //  Decryption
    /*---------------------------------------*/
    Mat<long> resmat;
    HELIB_NTIMER_START(dec_testTran);
    HEmatrix.decryptZmat(resmat, Tctxt);
    HELIB_NTIMER_STOP(dec_testTran);
    printNamedTimer(cout, "dec_testTran");
    
    cout << "Decrypt correct: " << Tctxt.isCorrect() << endl;
    
    /*---------------------------------------*/
    //  Result Comparison
    /*---------------------------------------*/
    cout << "Decrypted_matrix after Transpose: " << endl;
    cout << resmat << endl;
    
//    Mat<long> plain_trans;
//    transpose(plain_trans, Amat);
//    cout << "Plain_matrix after Transpose: " << endl;
//    cout << plain_trans << endl;
}


void HEmatrixTest::testShift(long nrows, long k) {

    long ncols = nrows;
    
    //Params param(/*m=*/20480, /*p=*/127, /*r=*/1, /*bits=*/200, /*c=*/2);
    Params param(/*m=*/20191, /*p=*/47, /*r=*/1, /*bits=*/200, /*c=*/2);
    
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, po2dim, neededslots) = (" << nrows << ", " << HEmatpar.dim  << ", " << HEmatpar.sqrdim << ", " << HEmatpar.neededslots << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    
    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    
    Meta meta;
    meta(param);
    HEmatrix HEmatrix(HEmatpar, meta);
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
    Mat<long> Amat;
    Amat.SetDims(nrows, ncols);
        
    for(int i = 0; i < nrows ; i++){
        for(int j = 0; j < ncols; j++){
            Amat[i][j]= j + 1;
        }
    }
    cout << "Generated matrix: " << Amat << endl;
    
    /*---------------------------------------*/
    //  Encryption
    /*---------------------------------------*/
    
    Ctxt Actxt(meta.data->publicKey);
    HEmatrix.encryptZmat(Actxt, Amat);
    
    /*---------------------------------------*/
    //  Transposition
    /*---------------------------------------*/
    
    zzX* shiftpoly;
    HELIB_NTIMER_START(genShift_testShiftPoly);
    HEmatrix.genShiftPoly(shiftpoly);
    HELIB_NTIMER_STOP(genShift_testShiftPoly);
    printNamedTimer(cout, "genShift_testShiftPoly");

    Ctxt Sctxt(meta.data->publicKey);
    
    HELIB_NTIMER_START(shift_testShiftPoly);
    HEmatrix.shiftBycols(Sctxt, Actxt, k, shiftpoly);
    HELIB_NTIMER_STOP(shift_testShiftPoly);
    printNamedTimer(cout, "shift_testShiftPoly");
    
    /*---------------------------------------*/
    //  Decryption
    /*---------------------------------------*/
    Mat<long> resmat;
    HEmatrix.decryptZmat(resmat, Sctxt);
    
    cout << "Decrypt iscorrect: " << Sctxt.isCorrect() << endl;
    
    /*---------------------------------------*/
    //  Result Comparison
    /*---------------------------------------*/
    cout << "Decrypted_matrix after Shift: " << endl;
    cout << resmat << endl;
}



void HEmatrixTest::testRMult(long nrows, long subdim) {
    long ncols = nrows;
    Params param(/*m=*/12800, /*p=*/127, /*r=*/1, /*bits=*/119, /*c=*/2);
    
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols, subdim);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, po2dim, neededslots) = (" << nrows << ", " << HEmatpar.dim  << ", " << HEmatpar.sqrdim << ", " << HEmatpar.neededslots << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    
    /*---------------------------------------*/
    //  Key Generation
    /*---------------------------------------*/
    
    Meta meta;
    meta(param);
    HEmatrix HEmatrix(HEmatpar, meta);
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
    Mat<long> rAmat;
    rAmat.SetDims(subdim, ncols);
    
    Mat<long> Amat;
    Amat.SetDims(nrows, ncols);
    
    Mat<long> Bmat;
    Bmat.SetDims(nrows, ncols);
    
    Mat<long> rBmat;
    rBmat.SetDims(nrows, subdim);
    
    for(long i = 0; i < subdim; i++){
        for(long j = 0; j < ncols; j++){
            rAmat[i][j]= (i % 2);
        }
    }
    
    for(long i = 0; i < nrows ; i++){
        for(long j = 0; j < ncols; j++){
            Amat[i][j]= ((rand() % 2));;
            Bmat[i][j]= ((rand() % 2));
        }
    }
    
    for(long i = 0; i < nrows; i++){
        for(long j = 0; j < subdim; j++){
            rBmat[i][j]= (i % 2);
        }
    }
    
    // replicate
    Mat<long> rAmat1;
    Mat<long> rBmat1;
    rAmat1.SetDims(nrows, ncols);
    rBmat1.SetDims(nrows, ncols);
    for(long i = 0; i < nrows/subdim; ++i){
        for(long j = 0; j < subdim; ++j){
            for(long k = 0; k < ncols; k++){
                rAmat1[i * subdim + j][k] = rAmat[j][k];
            }
        }
    }
    for(long i = 0; i < ncols/subdim; ++i){
        for(long j = 0; j < nrows; ++j){
            for(long k = 0; k < subdim; k++){
                rBmat1[j][i * subdim + k] = rBmat[j][k];
            }
        }
    }
    
    cout << "Generated matrixA: " << Amat << endl;
    cout << "Generated matrixRA: " << rAmat << endl;
    cout << "Generated matrixB: " << Bmat << endl;
    cout << "Generated matrixRB: " << rBmat << endl;
    cout << "Generated matrixRA replicate: " << rAmat1 << endl;
    cout << "Generated matrixRB replicate: " << rBmat1 << endl;
    
    /*---------------------------------------*/
    //  Encryption
    /*---------------------------------------*/
    
    //Note: Mat<long> cannot use transpose function, have to self-implement
//    Mat<ZZ> tempA;
//    transpose(tempA, Amat);
//
//    Mat<ZZ> tempB;
//    transpose(tempB, rBmat1);
    
    HELIB_NTIMER_START(ENc_testRMul);
    Ctxt Actxt(meta.data->publicKey);
    HEmatrix.encryptZmat(Actxt, rAmat1);
    Ctxt Bctxt(meta.data->publicKey);
    HEmatrix.encryptZmat(Bctxt, Bmat);
    HELIB_NTIMER_STOP(ENc_testRMul);
    printNamedTimer(cout, "ENc_testRMul");
    
    /*---------------------------------------*/
    //  GenPoly
    /*---------------------------------------*/
    zzX** Initpoly;
    HEmatrix.genMultPoly(Initpoly);
    
    zzX* shiftpoly;
    HEmatrix.genShiftPoly(shiftpoly);
    
    /*---------------------------------------*/
    //  Mult
    /*---------------------------------------*/
    HELIB_NTIMER_START(RMul_testRMul);
    Ctxt Cctxt(meta.data->publicKey);
    HEmatrix.HErmatmul(Cctxt, Actxt, Bctxt, Initpoly, shiftpoly);
    HELIB_NTIMER_STOP(RMul_testRMul);
    printNamedTimer(cout, "RMul_testRMul");
    
    /*---------------------------------------*/
    //  Decryption
    /*---------------------------------------*/
    Mat<long> HEresmat;
    HEmatrix.decryptZmat(HEresmat, Cctxt);
    
    /*---------------------------------------*/
    //  Result Comparison
    /*---------------------------------------*/
//    Mat<long> resmat;
//    mul(resmat, rAmat, Bmat);
//
//    Mat<long> trans;
//    transpose(trans, HEresmat);
//
//    cout << "------------------" << endl;
//    cout << "Plaintext Multiplication:" << endl;
//    cout << resmat << endl;
    cout << "------------------" << endl;
    cout << "Encrypted Multiplication:" << endl;
    cout << HEresmat << endl;
    cout << "------------------" << endl;
}

void HEmatrixTest::testMult_preprocessing(long nrows) {
    
    long ncols = nrows;
    
//    Params param(/*m=*/12800, /*p=*/127, /*r=*/1, /*bits=*/119, /*c=*/2);
    // Params param(/*m=*/17408, /*p=*/127, /*r=*/1, /*bits=*/180, /*c=*/2);
    Params param(/*m=*/16384, /*p=*/8191, /*r=*/1, /*bits=*/120, /*c=*/2);
    
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, po2dim, neededslots) = (" << nrows << ", " << HEmatpar.dim  << ", " << HEmatpar.sqrdim << ", " << HEmatpar.neededslots << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    
    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    
    Meta meta;
    meta(param);
    HEmatrix HEmatrix(HEmatpar, meta);
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
    Mat<long> Amat;
    Mat<long> Bmat;
    
    Amat.SetDims(nrows, ncols);
    Bmat.SetDims(nrows, ncols);
    
    for(long i = 0; i < nrows ; i++){
        for(long j = 0; j < ncols; j++){
            Amat[i][j]= (i % 2);
            Bmat[i][j]= (i % 3);
        }
    }
    // cout << "Generated matrixA: " << Amat << endl;
    // cout << "Generated matrixB: " << Bmat << endl;
    
    /*---------------------------------------*/
    //  Encryption
    /*---------------------------------------*/
    
    HELIB_NTIMER_START(Enc_testMul_preprocessing);
    vector<Ctxt> Actxts;
    HEmatrix.genInitActxt(Actxts, Amat);
    
    Ctxt Bctxt(meta.data->publicKey);
    HEmatrix.encryptZmat(Bctxt, Bmat);
    HELIB_NTIMER_STOP(Enc_testMul_preprocessing);
    printNamedTimer(cout, "Enc_testMul_preprocessing");

    stringstream ss;
    for (int i = 0; i < Actxts.size(); i++)
    {
        Actxts[i].writeTo(ss);
    }
    string ClientTemp = ss.str();
    uint64_t client_totalLength = ClientTemp.size();
    
    cout << "Ciphertext size = : " << ((double) client_totalLength / (double)(1024 * 1024)) << " MB" << endl;
    
    /*---------------------------------------*/
    //  GenPoly
    /*---------------------------------------*/
    
    zzX* Initpoly;
    HELIB_NTIMER_START(GenMulPoly_testMul_preprocessing);
    HEmatrix.genMultBPoly(Initpoly);
    HELIB_NTIMER_STOP(GenMulPoly_testMul_preprocessing);
    printNamedTimer(cout, "GenMulPoly_testMul_preprocessing");
    
    /*---------------------------------------*/
    //  Mult
    /*---------------------------------------*/
    
    HELIB_NTIMER_START(Mul_testMul_preprocessing);
    Ctxt Cctxt(meta.data->publicKey);
    HEmatrix.HEmatmul_preprocessing(Cctxt, Actxts, Bctxt, Initpoly);
    HELIB_NTIMER_STOP(Mul_testMul_preprocessing);
    printNamedTimer(cout, "Mul_testMul_preprocessing");
    
    /*---------------------------------------*/
    //  Decryption
    /*---------------------------------------*/
    Mat<long> HEresmat;
    HEmatrix.decryptZmat(HEresmat, Cctxt);
    
    /*---------------------------------------*/
    //  Result Comparison
    /*---------------------------------------*/
//    Mat<long> resmat;
//    mul(resmat, Amat, Bmat);
//
//    cout << "------------------" << endl;
//    cout << "Plaintext Multiplication:" << endl;
//    cout << resmat << endl;
    cout << "------------------" << endl;
    cout << "Encrypted Multiplication:" << endl;
    cout << HEresmat << endl;
    cout << "------------------" << endl;
}

void HEmatrixTest::testRMult_preprocessing(long nrows, long subdim) {
    
    long ncols = nrows;
    Params param(/*m=*/20480, /*p=*/127, /*r=*/1, /*bits=*/200, /*c=*/2);
    
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols, subdim);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, po2dim, neededslots) = (" << nrows << ", " << HEmatpar.dim  << ", " << HEmatpar.sqrdim << ", " << HEmatpar.neededslots << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    
    /*---------------------------------------*/
    //  Key Generation
    /*---------------------------------------*/
    
    Meta meta;
    meta(param);
    HEmatrix HEmatrix(HEmatpar, meta);
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
    Mat<long> rAmat;
    rAmat.SetDims(subdim, ncols);
    
    Mat<long> Bmat;
    Bmat.SetDims(nrows, ncols);
    
    for(long i = 0; i < subdim; i++){
        for(long j = 0; j < ncols; j++){
            rAmat[i][j]= (i % 2);
        }
    }
    
    for(long i = 0; i < nrows ; i++){
        for(long j = 0; j < ncols; j++){
            Bmat[i][j]= (i % 2);;
        }
    }
    cout << "Generated matrixRA: " << rAmat << endl;
    cout << "Generated matrixB: " << Bmat << endl;
    
    /*---------------------------------------*/
    //  Encryption
    /*---------------------------------------*/
    
    HELIB_NTIMER_START(Enc_testRMulPreprocess);
    vector<Ctxt> Actxts;
    HEmatrix.genInitRecActxt(Actxts, rAmat);
    
    Ctxt Bctxt(meta.data->publicKey);
    HEmatrix.encryptZmat(Bctxt, Bmat);
    HELIB_NTIMER_STOP(Enc_testRMulPreprocess);
    printNamedTimer(cout, "Enc_testRMulPreprocess");
    
    /*---------------------------------------*/
    //  GenPoly
    /*---------------------------------------*/
    
    zzX* Initpoly;
    HEmatrix.genMultBPoly(Initpoly);
    
    /*---------------------------------------*/
    //  Mult
    /*---------------------------------------*/
    
    HELIB_NTIMER_START(RMul_testRMulPreprocess);
    Ctxt Cctxt(meta.data->publicKey);
    HEmatrix.HErmatmul_preprocessing(Cctxt, Actxts, Bctxt, Initpoly);
    HELIB_NTIMER_STOP(RMul_testRMulPreprocess);
    printNamedTimer(cout, "RMul_testRMulPreprocess");
    
    /*---------------------------------------*/
    //  Decryption
    /*---------------------------------------*/
    Mat<long> HEresmat;
    HEmatrix.decryptZmat(HEresmat, Cctxt);
    
    /*---------------------------------------*/
    //  Result Comparison
    /*---------------------------------------*/
//    Mat<ZZ> resmat;
//    mul(resmat, rAmat, Bmat);
//
//    cout << "------------------" << endl;
//    cout << "Plaintext Multiplication:" << endl;
//    cout << resmat << endl;
    cout << "------------------" << endl;
    cout << "Encrypted Multiplication:" << endl;
    cout << HEresmat << endl;
    cout << "------------------" << endl;
}

//----------------------------------------------------
// SIMD (parallel computation)

void HEmatrixTest::testSIMDAdd(long nrows, long nbatching, const long niter) {
    
    long ncols = nrows;
    Params param(/*m=*/20480, /*p=*/127, /*r=*/1, /*bits=*/200, /*c=*/2);
    
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols, 0, nbatching);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, nslots, nbatching) = (" << nrows << "," << HEmatpar.dim  << "," << HEmatpar.neededslots << "," << HEmatpar.nbatching << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    
    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    
    HELIB_NTIMER_START(SIMDAdd_setup_time);
    Meta meta;
    meta(param);
    cout << "Context contents: " << endl;
    meta.data->context.printout();
    HEmatrix HEmatrix(HEmatpar, meta);
    HELIB_NTIMER_STOP(SIMDAdd_setup_time);
    printNamedTimer(cout, "SIMDAdd_setup_time");
    
    /*---------------------------------------*/
    //  Generate two random matrixs
    /*---------------------------------------*/
    Mat<long>* Amat = new Mat<long>[nbatching];
    Mat<long>* Bmat = new Mat<long>[nbatching];
    
    for(long k = 0; k < nbatching; ++k){
        Amat[k].SetDims(nrows, ncols);
        Bmat[k].SetDims(nrows, ncols);
        
        for(long i = 0; i < nrows ; i++){
            for(long j = 0; j < ncols; j++){
                Amat[k][i][j] = ((i * ncols + j));
                Bmat[k][i][j] = ((i * ncols + j));
            }
        }
    }
    for(long k = 0; k < nbatching; ++k){
        cout << "Generated matrixA " << k << "-th: " << endl;
        cout << Amat[k] << endl;
        cout << "Generated matrixB " << k << "-th: " << endl;
        cout << Bmat[k] << endl;
    }
    
    for(long l = 0; l < niter; ++l){
        /*---------------------------------------*/
        //  Encryption
        /*---------------------------------------*/
        
        HELIB_NTIMER_START(EncM_Mat_time);
        Ctxt ctxtA(meta.data->publicKey);
        HEmatrix.encryptParallelZmat(ctxtA, Amat, nbatching);
        Ctxt ctxtB(meta.data->publicKey);
        HEmatrix.encryptParallelZmat(ctxtB, Bmat, nbatching);
        HELIB_NTIMER_STOP(EncM_Mat_time);
        printNamedTimer(cout, "EncM_Mat_time");
        
        /*---------------------------------------*/
        //  Addition
        /*---------------------------------------*/
        
        HELIB_NTIMER_START(SIMDAdd_time);
        ctxtA += ctxtB;
        HELIB_NTIMER_STOP(SIMDAdd_time);
        printNamedTimer(cout, "SIMDAdd_time");
        
        /*---------------------------------------*/
        //  Decryption
        /*---------------------------------------*/
        Mat<long>* resmat;
        HELIB_NTIMER_START(DecM_Mat_time);
        HEmatrix.decryptParallelZmat(resmat, ctxtA);
        HELIB_NTIMER_STOP(DecM_Mat_time);
        printNamedTimer(cout, "DecM_Mat_time");
        
        for(long k = 0; k < nbatching; ++k){
            cout << "Decrypted matrix_" << k << ": " << resmat[k] << endl;
        }
    }
}

void HEmatrixTest::testSIMDTrans(long nrows, long nbatching, const long niter) {
    
    long ncols = nrows;
    
    Params param(/*m=*/20480, /*p=*/127, /*r=*/1, /*bits=*/200, /*c=*/2);
    
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols, 0, nbatching);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, nslots, nbatching) = (" << nrows << "," << HEmatpar.dim  << "," << HEmatpar.neededslots << "," << HEmatpar.nbatching << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    
    /*---------------------------------------*/
    //  Setup
    /*---------------------------------------*/
    
    Meta meta;
    meta(param);
    HEmatrix HEmatrix(HEmatpar, meta);
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
    Mat<long>* Amat = new Mat<long>[nbatching];
   
    for(long k = 0; k < nbatching; ++k){
        Amat[k].SetDims(nrows, ncols);
        for(long i = 0; i < nrows ; i++){
            for(long j = 0; j < ncols; j++){
                Amat[k][i][j] = ((i + 1));
            }
        }
    }
    for(long k = 0; k < nbatching; ++k){
        cout << "Generated matrixA_" << k << ": " << endl;
        cout << Amat[k] << endl;
    }
    
    for(long l = 0; l < niter; ++l){
        /*---------------------------------------*/
        //  Encryption
        /*---------------------------------------*/
        
        HELIB_NTIMER_START(EncM_testTransParallel);
        Ctxt Actxt(meta.data->publicKey);
        HEmatrix.encryptParallelZmat(Actxt, Amat, nbatching);
        HELIB_NTIMER_STOP(EncM_testTransParallel);
        printNamedTimer(cout, "EncM_testTransParallel");
        
        /*---------------------------------------*/
        //  Transpose
        /*---------------------------------------*/
        zzX* transpoly;
        HEmatrix.genTransPoly_Parallel(transpoly);
        
        Ctxt Tctxt(meta.data->publicKey);
        HELIB_NTIMER_START(tran_testSIMDTran);
        HEmatrix.transpose_Parallel(Tctxt, Actxt, transpoly);
        HELIB_NTIMER_STOP(tran_testSIMDTran);
        printNamedTimer(cout, "tran_testSIMDTran");
        
        /*---------------------------------------*/
        //  Decryption
        /*---------------------------------------*/
        Mat<long>* HEresmat;
        HEmatrix.decryptParallelZmat(HEresmat, Tctxt);
        
        /*---------------------------------------*/
        //  Result Comparison
        /*---------------------------------------*/
//        Mat<ZZ>* resmat = new Mat<ZZ>[nbatching];
//        for(long k = 0; k < nbatching; ++k){
//            transpose(resmat[k], Amat[k]);
//            cout << "Plaintext Matrix_" << k <<  " Transpose:" << endl;
//            cout << resmat[k] << endl;
//            cout << "++++++++++++++++++++++++++++++++++++++++++" << endl;
//            cout << "Decrypted matrix_" << k << ": " << endl;
//            cout << HEresmat[k] << endl;
//        }
    }
}

void HEmatrixTest::testSIMDMult(long nrows, long nbatching, const long niter) {
    
    long ncols = nrows;
    
    Params param(/*m=*/20480, /*p=*/127, /*r=*/1, /*bits=*/200, /*c=*/2);
    
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols, 0, nbatching);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, nslots, nbatching) = (" << nrows << "," << HEmatpar.dim  << "," << HEmatpar.neededslots << "," << HEmatpar.nbatching << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    
    /*---------------------------------------*/
    //  Key Generation
    /*---------------------------------------*/
    
    Meta meta;
    meta(param);
    HEmatrix HEmatrix(HEmatpar, meta);
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
    Mat<long>* Amat = new Mat<long>[nbatching];
    Mat<long>* Bmat = new Mat<long>[nbatching];
    
    for(long k = 0; k < nbatching; ++k){
        Amat[k].SetDims(nrows, ncols);
        Bmat[k].SetDims(nrows, ncols);
        
        for(long i = 0; i < nrows ; i++){
            for(long j = 0; j < ncols; j++){
                Amat[k][i][j] = ((i % 2));
                Bmat[k][i][j] = ((j % 2));
            }
        }
    }
    for(long k = 0; k < nbatching; ++k){
        cout << "Generated matrixA " << k << "-th: " << endl;
        cout << Amat[k] << endl;
        cout << "Generated matrixB " << k << "-th: " << endl;
        cout << Bmat[k] << endl;
    }
    
    for(long l = 0; l < niter; ++l){
        /*---------------------------------------*/
        //  Encryption
        /*---------------------------------------*/
        
        HELIB_NTIMER_START(Enc_testSIMDMul);
        Ctxt ctxtA(meta.data->publicKey);
        HEmatrix.encryptParallelZmat(ctxtA, Amat, nbatching);
        Ctxt ctxtB(meta.data->publicKey);
        HEmatrix.encryptParallelZmat(ctxtB, Bmat, nbatching);
        HELIB_NTIMER_STOP(Enc_testSIMDMul);
        printNamedTimer(cout, "Enc_testSIMDMul");
        
        /*---------------------------------------*/
        //  Genpoly
        /*---------------------------------------*/
        
        zzX** Initpoly;
        HEmatrix.genMultPoly_Parallel(Initpoly);
        
        zzX* shiftpoly;
        HEmatrix.genShiftPoly_Parallel(shiftpoly);
        
        /*---------------------------------------*/
        //  Mult
        /*---------------------------------------*/
        
        Ctxt Cctxt(meta.data->publicKey);
        HELIB_NTIMER_START(Mul_testSIMDMul);
        HEmatrix.HEmatmul_Parallel(Cctxt, ctxtA, ctxtB, Initpoly, shiftpoly);
        HELIB_NTIMER_STOP(Mul_testSIMDMul);
        printNamedTimer(cout, "Mul_testSIMDMul");
        
        /*---------------------------------------*/
        //  Decryption
        /*---------------------------------------*/
        
        Mat<long>* HEresmat;
        HEmatrix.decryptParallelZmat(HEresmat, Cctxt);
        
        /*---------------------------------------*/
        //  Result Comparison
        /*---------------------------------------*/
//        Mat<long>* resmat = new Mat<long>[nbatching];
        for(long k = 0; k < nbatching; ++k){
//            mul(resmat[k], Amat[k], Bmat[k]);
//            cout << "Plaintext Matrix_" << k <<  " Mul:" << endl;
//            cout << resmat[k] << endl;
            cout << "Decrypted matrix_" << k << " Mul: " << endl;
            cout << HEresmat[k] << endl;
        }
    }
}

void HEmatrixTest::testSIMDRMult(long nrows, long subdim, long nbatching, const long niter){
    long ncols = nrows;
    
    Params param(/*m=*/20480, /*p=*/127, /*r=*/1, /*bits=*/200, /*c=*/2);
    
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols, subdim, nbatching);
    
    cout << "-----------------------------------" << endl;
    cout << "(rows, dim, nslots, nbatching) = (" << nrows << "," << HEmatpar.dim  << "," << HEmatpar.neededslots << "," << HEmatpar.nbatching << ")" << endl;
    cout << "-----------------------------------" << endl;
    
    /*---------------------------------------*/
    //  Key Generation
    /*---------------------------------------*/
    
    Meta meta;
    meta(param);
    HEmatrix HEmatrix(HEmatpar, meta);
    
    /*---------------------------------------*/
    //  Generate a random matrix
    /*---------------------------------------*/
    Mat<long> rAmat;
    rAmat.SetDims(subdim, ncols);
    
    Mat<long>* Bmat = new Mat<long>[nbatching];
    
    for(long i = 0; i < subdim; i++){
        for(long j = 0; j < ncols; j++){
            rAmat[i][j]= ((i % 2) + 1);
        }
    }
    
    for(long k = 0; k < nbatching; ++k){
        Bmat[k].SetDims(nrows, ncols);
        
        for(long i = 0; i < nrows ; i++){
            for(long j = 0; j < ncols; j++){
                Bmat[k][i][j] = ((j % 2) + 1);
            }
        }
    }
    
    // replicate
    Mat<long>* rAmat1 = new Mat<long>[nbatching];
    for(long k = 0; k < nbatching; ++k){
        rAmat1[k].SetDims(nrows, ncols);
        for(long i = 0; i < nrows/subdim; ++i){
            for(long j = 0; j < subdim; ++j){
                for(long l = 0; l < ncols; l++){
                    rAmat1[k][i * subdim + j][l] = rAmat[j][l];
                }
            }
        }
    }
    
    cout << "Generated matrixRA: " << rAmat << endl;
    for(long k = 0; k < nbatching; ++k){
        cout << "Generated matrixB " << k << "-th: " << endl;
        cout << Bmat[k] << endl;
        cout << "Generated matrixRAmat1 " << k << "-th: " << endl;
        cout << rAmat1[k] << endl;
    }
    for(long l = 0; l < niter; ++l){
        /*---------------------------------------*/
        //  Encryption
        /*---------------------------------------*/
        Ctxt Actxt(meta.data->publicKey);
        HEmatrix.encryptParallelZmat(Actxt, rAmat1, nbatching);
        Ctxt Bctxt(meta.data->publicKey);
        HEmatrix.encryptParallelZmat(Bctxt, Bmat, nbatching);
        
        /*---------------------------------------*/
        //  GenPoly
        /*---------------------------------------*/
        zzX** Initpoly;
        HEmatrix.genMultPoly_Parallel(Initpoly);
        
        zzX* shiftpoly;
        HEmatrix.genShiftPoly_Parallel(shiftpoly);
        
        /*---------------------------------------*/
        //  Mult
        /*---------------------------------------*/
        
        Ctxt Cctxt(meta.data->publicKey);
        HELIB_NTIMER_START(SIMDRMul_testSIMDRMul);
        HEmatrix.HErmatmul_Parallel(Cctxt, Actxt, Bctxt, Initpoly, shiftpoly);
        HELIB_NTIMER_STOP(SIMDRMul_testSIMDRMul);
        printNamedTimer(cout, "SIMDRMul_testSIMDRMul");
        
        /*---------------------------------------*/
        //  Decryption
        /*---------------------------------------*/
        Mat<long>* HEresmat;
        HEmatrix.decryptParallelZmat(HEresmat, Cctxt);
        
        /*---------------------------------------*/
        //  Result Comparison
        /*---------------------------------------*/
//        Mat<ZZ>* resmat = new Mat<ZZ>[nbatching];
//        for(long k = 0; k < nbatching; ++k){
//            mul(resmat[k], rAmat, Bmat[k]);
//            cout << "Plaintext Matrix_" << k <<  " Mul:" << endl;
//            cout << resmat[k] << endl;
//            cout << "Decrypted matrix_" << k << " Mul: " << endl;
//            cout << HEresmat[k] << endl;
//        }
    }
}

//Fixme: Mat<ZZ> can be replaced with self-defined matrix.h inside the enclave
bool HEmatrixTest::LoadData(Mat<long> &rawData, long &rawdim, long &padzeros, const string &filename){
    ifstream fin;
    fin.open(filename);
    if (!fin) {
        cout << "Unable to read data file." << endl;
        return false;
    }
    
    long n;
    int data; // 0 or 1
    int numAttr;
    fin >> numAttr >> rawdim >> n;
    
    //int padzeros = 0;
    long dimpow = 0;
    for (int i = 1; i < 7; i++) {
        dimpow = 1 << i;
        if (rawdim < dimpow) {
            padzeros = dimpow - rawdim;
            //rawdim = dimpow;
            break;
        }
    }
    rawData.SetDims(n, dimpow);
    
    for (int i = 0; i < n; i++) {
        for (unsigned j = 0; j < rawdim; j++) {
            fin >> data;
            rawData[i][j] = (data);
        }
    }
    // pad zeros
    for (long i = 0; i < n; i++) {
        for (long j = rawdim; j < dimpow; j++) {
            rawData[i][j] = 0;
        }
    }
    return true;
}

void HEmatrixTest::Mult_Rawdata(){
    Mat<long> rawData;
    long rawdim = 0;
    long padzeros = 0;
    const string datafile = "/Users/qifanwang/code/ePPDSC/ePPDSC/scripts/onehot_kohkiloyeh.dat";
    if (!HEmatrixTest::LoadData(rawData, rawdim, padzeros, datafile)) {
        return;
    }
    cout << "rawData: " << endl;
    cout << rawData << endl;
    
    long nrows = rawdim + padzeros;
    long ncols = nrows;
    long subdim = 4;
    cout << "nrows: " << nrows << endl;
    Params param(/*m=*/20480, /*p=*/127, /*r=*/1, /*bits=*/200, /*c=*/2);
    HEMatpar HEmatpar;
    readHEMatpar(HEmatpar, nrows, ncols, subdim);
    Meta meta;
    meta(param);
    HEmatrix HEmatrix(HEmatpar, meta);
    
    long numMat = 0;
    Mat<long>* Amat;
    if (rawData.NumRows()%nrows != 0) {
        //last matrix is padded with zero, in ePPDSC, we can set the buffer size to the power of two to avoid pad zero in the last matrix
        numMat = rawData.NumRows()/nrows + 1;
        cout << "numMat: " << numMat << endl;
        Amat = new Mat<long>[numMat];
        NTL_EXEC_RANGE(numMat - 1, first, last);
        for(long k = first; k < last; ++k){
            Amat[k].SetDims(nrows, ncols);
            for(long i = 0; i < nrows ; i++){
                for(long j = 0; j < ncols; j++){
                    Amat[k][i][j] = rawData[k * nrows + i][j];
                }
            }
        }
        NTL_EXEC_RANGE_END;
        Amat[numMat - 1].SetDims(nrows, ncols);
        for(long i = 0; i < nrows ; i++){
            for(long j = 0; j < ncols; j++){
                if (((numMat - 1) * nrows + i) < rawData.NumRows()){
                    Amat[numMat - 1][i][j] = rawData[(numMat - 1) * nrows + i][j];
                }else{
                    Amat[numMat - 1][i][j] = 0;
                }
            }
        }
    }
    else{
        //dont need to pad with zero
        numMat = rawData.NumRows()/nrows;
        Amat = new Mat<long>[numMat];
        NTL_EXEC_RANGE(numMat, first, last);
        for(long k = first; k < last; ++k){
            Amat[k].SetDims(nrows, ncols);
            for(long i = 0; i < nrows ; i++){
                for(long j = 0; j < ncols; j++){
                    Amat[k][i][j] = rawData[k * nrows + i][j];
                }
            }
        }
        NTL_EXEC_RANGE_END;
    }
    
    Mat<long> rBmat;
    rBmat.SetDims(nrows, subdim);
    for(long i = 0; i < nrows ; i++){
        for(long j = 0; j < subdim; j++){
            rBmat[i][j]= (i % 2);
        }
    }
    // replicate
    Mat<long> Bmat;
    Bmat.SetDims(nrows, ncols);
    for(long i = 0; i < ncols/subdim; ++i){
        for(long j = 0; j < nrows; ++j){
            for(long k = 0; k < subdim; k++){
                Bmat[j][i * subdim + k] = rBmat[j][k];
            }
        }
    }
    cout << "Generated matrixRB: " << rBmat << endl;
    cout << "Generated matrixB: " << Bmat << endl;
    
    Mat<long> Bmat_trans;
//    transpose(Bmat_trans, Bmat);
    Ctxt Bctxt(meta.data->publicKey);
    HEmatrix.encryptZmat(Bctxt, Bmat);
    
    cout << "generate poly:" << endl;
    zzX** Initpoly;
    HEmatrix.genMultPoly(Initpoly);
    
    zzX* shiftpoly;
    HEmatrix.genShiftPoly(shiftpoly);
    
    //encrypt the matrix
    cout << "encryption, mul, decryption: " << endl;
    vector<Ctxt> ctxtA(numMat, Ctxt(meta.data->publicKey));
    vector<Ctxt> Cctxt(numMat, Ctxt(meta.data->publicKey));
    Mat<long>* HEresmat = new Mat<long>[numMat];
    HELIB_NTIMER_START(Mul_testMulRawData);
    NTL_EXEC_RANGE(numMat, first, last);
    for(long i = first; i < last; ++i){
        HEmatrix.encryptZmat(ctxtA[i], Amat[i]);
        HEmatrix.HErmatmul(Cctxt[i], ctxtA[i], Bctxt, Initpoly, shiftpoly);
        HEmatrix.decryptZmat(HEresmat[i], Cctxt[i]);
    }
    NTL_EXEC_RANGE_END;
    HELIB_NTIMER_STOP(Mul_testMulRawData);
    printNamedTimer(cout, "Mul_testMulRawData");
    
//    Mat<ZZ>* resmat = new Mat<ZZ>[numMat];
//    for (int i = 0; i < numMat; i++) {
//        mul(resmat[i], Amat[i], rBmat);
//        cout << i << "-th Plaintext Multiplication:" << endl;
//        cout << resmat[i] << endl;
//        cout << "------------------" << endl;
//        cout << i << "-th Encrypted Multiplication:" << endl;
//        cout << HEresmat[i] << endl;
//    }
}
