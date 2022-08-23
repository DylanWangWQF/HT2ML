//
//  main.cpp
//  cloneHElib
//
//  Created by Qifan Wang on 19/04/21.
//

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <sys/time.h>

// BGV Matrix tests
#include "hematrix.h"
#include "HEmatrix_test.h"

// CKKS Matrix tests
#include "CKKSmatrix.h"
#include "CKKSmatrix_test.h"

// Tests
#include "../Tests/BGV_test.h"
#include "../Tests/Unit_Test.h"
#include "../Tests/Noise_Test.h"
#include "../Tests/CKKS_Test.h"
#include "../Tests/Bootstrapping_params.h"

// Linear Regression
#include "../LR/Test_Regression.h"

using namespace std;

static test tes;
int main(int argc, const char * argv[]) {
//    tes.testStruct();
//    tes.testTrans();
//    tes.close();
//    return 0;
    SetNumThreads(8);

    while (true) {
        int selection = 0;
        bool invalid = true;
        do
        {
            cout << endl << "> Run example (1~40) or exit (0): ";
            if (!(cin >> selection))
            {
                invalid = false;
            }
            else if (selection < 0 || selection > 40)
            {
                invalid = false;
            }
            else
            {
                invalid = true;
            }
            if (!invalid)
            {
                cout << "  Invalid option: type 0 ~ 40" << endl;
                cin.clear();
            }
        } while (!invalid);

        switch (selection)
        {
            case 1:
                cout << "test BGV_binary_arithmetric starts: " << endl;
                BGVTest::BGV_binary_arithmetric();
                break;

            case 2:
                cout << "test BGV_packed_arithmetric starts: " << endl;
                BGVTest::BGV_packed_arithmetric();
                break;
            case 3:
                cout << "test enc matrix: " << endl;
                HEmatrixTest::testEnc(16);
                break;
            case 4:
                cout << "test encryptedarray functions starts: " << endl;
                BGVTest::encryptedarray_test();
                break;
            case 5:
                cout << "test SIMDAddWithMultipleMatrixs starts: " << endl;
                HEmatrixTest::testSIMDAdd(4, 2, 1);
                break;
            case 6:
                cout << "test FindM: " << endl;
                HEmatrixTest::FindCorrectM();
                break;
            case 7:
                cout << "test Single Matrix Add: " << endl;
                HEmatrixTest::testAdd(16);
                break;
            case 8:
                cout << "test Single Matrix Transpose: " << endl;
                HEmatrixTest::testTrans(16);
                break;
            case 9:
                cout << "test Rotation: " << endl;
                HEmatrixTest::testRotation();
                break;
            case 10:
                cout << "test ShiftByCol: " << endl;
                HEmatrixTest::testShift(16, 3);
                break;
            case 11:
                cout << "test Mul: " << endl;
                HEmatrixTest::testMult(128);
                break;
            case 12:
                cout << "test RMul: " << endl;
                HEmatrixTest::testRMult(16, 4);
                break;
            case 13:
                cout << "test Mult_preprocessing: " << endl;
                HEmatrixTest::testMult_preprocessing(32);
                break;
            case 14:
                cout << "test RMult_preprocessing: " << endl;
                HEmatrixTest::testRMult_preprocessing(16, 4);
                break;
            case 15:
                cout << "test SIMD transpose: " << endl;
                HEmatrixTest::testSIMDTrans(8, 4, 1);
                break;
            case 16:
                cout << "test SIMD Mul: " << endl;
                HEmatrixTest::testSIMDMult(8, 4, 1);
                break;
            case 17:
                cout << "test SIMD RMul: " << endl;
                HEmatrixTest::testSIMDRMult(8, 4, 4, 1);
                break;
            case 18:
                cout << "test Mul with rawdata: " << endl;
                HEmatrixTest::Mult_Rawdata();
                break;
            case 19:
                cout << "test noise estimation: " << endl;
                NoiseTest::Test_Noise();
                break;
            case 20:
                cout << "test CKKS basics: " << endl;
                CKKSTest::CKKS_basics();
                break;
            case 21:
                cout << "test CKKS depth: " << endl;
                CKKSTest::CKKS_depth();
                break;
            case 22:
                cout << "test CKKS data move: " << endl;
                CKKSTest::CKKS_data_move();
                break;
            case 23:
                cout << "test CKKS matmul: " << endl;
                CKKSTest::CKKS_matmul();
                break;
            case 24:
                cout << "test bootstrapping with mul in BGV: " << endl;
                BGVTest::Bootstrapping_Test();
                break;
            case 25:
                cout << "Find ParamsForBootstrapping: " << endl;
                BootsParams::Find_ParamsForBootstrapping();
                break;
            case 26:
                cout << "Check Context information: " << endl;
                BGVTest::Check_Context();
                break;
            case 27:
                cout << "Linear Regression: " << endl;
                LRTest::RunRegressionTest();
                break;
            case 28:
                cout << "Test CKKS Matrix Enc: " << endl;
                CKKSmatrixTest::testEnc(64);
                break;
            case 29:
                cout << "Test CKKS basics: " << endl;
                CKKSmatrixTest::testCKKS();
                break;
            case 30:
                cout << "Test CKKS Matrix Add: " << endl;
                CKKSmatrixTest::testAdd(16);
                break;
            case 31:
                cout << "Test CKKS Matrix Transpose: " << endl;
                CKKSmatrixTest::testTrans(16);
                break;
            case 32:
                cout << "CKKSmatrixTest ShiftByCol: " << endl;
                CKKSmatrixTest::testShift(16, 3);
                break;
            case 33:
                cout << "CKKSmatrixTest testMult: " << endl;
                CKKSmatrixTest::testMult(64);
                break;
            case 34:
                cout << "CKKSmatrixTest testRMult: " << endl;
                CKKSmatrixTest::testRMult(16, 4);
                break;
            case 35:
                cout << "CKKSmatrixTest testMult_preprocessing: " << endl;
                CKKSmatrixTest::testMult_preprocessing(16);
                break;
            case 36:
                cout << "CKKSmatrixTest testRMult_preprocessing: " << endl;
                CKKSmatrixTest::testRMult_preprocessing(16, 4);
                break;
            case 37:
                cout << "CKKSmatrixTest testSIMDAdd: " << endl;
                CKKSmatrixTest::testSIMDAdd(4, 2, 1);
                break;

            case 0:
                return 0;
        }
    }
    return 0;
}
