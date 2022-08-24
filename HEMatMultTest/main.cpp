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
#include "BGVmatrix.h"
#include "BGVmatrix_test.h"

// CKKS Matrix tests
#include "CKKSmatrix.h"
#include "CKKSmatrix_test.h"

using namespace std;

int main(int argc, const char * argv[]) {
    SetNumThreads(16);

    while (true) {
        int selection = 0;
        bool invalid = true;
        do
        {
            cout << endl << "> Run example (1~8) or exit (0): ";
            if (!(cin >> selection))
            {
                invalid = false;
            }
            else if (selection < 0 || selection > 8)
            {
                invalid = false;
            }
            else
            {
                invalid = true;
            }
            if (!invalid)
            {
                cout << "  Invalid option: type 0 ~ 8" << endl;
                cin.clear();
            }
        } while (!invalid);

        switch (selection)
        {
            case 1:
                cout << "test Mul: " << endl;
                HEmatrixTest::testMult(64);
                break;
            case 2:
                cout << "test RMul: " << endl;
                HEmatrixTest::testRMult(16, 4);
                break;
            case 3:
                cout << "test Mult_preprocessing: " << endl;
                HEmatrixTest::testMult_preprocessing(64);
                break;
            case 4:
                cout << "test RMult_preprocessing: " << endl;
                HEmatrixTest::testRMult_preprocessing(16, 4);
                break;
            case 5:
                cout << "CKKSmatrixTest testMult: " << endl;
                CKKSmatrixTest::testMult(64);
                break;
            case 6:
                cout << "CKKSmatrixTest testRMult: " << endl;
                CKKSmatrixTest::testRMult(16, 4);
                break;
            case 7:
                cout << "CKKSmatrixTest testMult_preprocessing: " << endl;
                CKKSmatrixTest::testMult_preprocessing(64);
                break;
            case 8:
                cout << "CKKSmatrixTest testRMult_preprocessing: " << endl;
                CKKSmatrixTest::testRMult_preprocessing(16, 4);
                break;

            case 0:
                return 0;
        }
    }
    return 0;
}
