//
//  CKKSmatrix_test.h
//  ePPDSC
//
//  Created by Qifan Wang on 2/01/22.
//

#ifndef CKKSmatrix_test_h
#define CKKSmatrix_test_h

class CKKSmatrixTest{
public:
    static void testMult(long dim);
    static void testRMult(long nrows, long subdim);
    static void testMult_preprocessing(long dim);
    static void testRMult_preprocessing(long nrows, long subdim);
};

#endif /* CKKSmatrix_test_h */
