//
//  HEmatrix_test.h
//  cloneHElib
//
//  Created by Qifan Wang on 25/04/21.
//

#ifndef HEmatrix_test_h
#define HEmatrix_test_h

class HEmatrixTest {
//private:
//    testS* tes;
    
public:
    static void FindCorrectM();
    static void testRotation();
    //static void testBGV();
    static void testEnc(long dim);
    static void testAdd(long dim);
    static void testTrans(long nrows);
    static void testShift(long nrows, long k);
    static void testMult(long nrows);
    static void testRMult(long nrows, long subdim);
    static void testMult_preprocessing(long nrows);
    static void testRMult_preprocessing(long nrows, long subdim);
    static void testSIMDAdd(long dim, long nbatching, const long nitr);
    static void testSIMDTrans(long nrows, long nbatching, const long niter);
    static void testSIMDMult(long nrows, long nbatching, const long niter);
    static void testSIMDRMult(long nrows, long subdim, long nbatching, const long niter);
    static bool LoadData(Mat<long> &rawData, long &rawdim, long &padzeros, const string &filename);
    static void Mult_Rawdata();
};

#endif /* HEmatrix_test_h */
