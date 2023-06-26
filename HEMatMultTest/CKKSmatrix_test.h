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
