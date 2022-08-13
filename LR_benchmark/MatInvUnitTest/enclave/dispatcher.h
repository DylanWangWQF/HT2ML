#pragma once

#include <cstring>
#include <sstream>
#include <math.h>
#include <iostream>
// #include <openenclave/enclave.h>

// #define numAttr 2
// #define numAttr 4
// #define numAttr 8
// #define numAttr 6
#define numAttr 32

using namespace std;

class ecall_dispatcher
{
  public:
    ecall_dispatcher();
    int MatInv();
    void close();
  private:
    void getCofactor(int A[numAttr][numAttr], int temp[numAttr][numAttr], int p, int q, int n);
    int determinant(int A[numAttr][numAttr], int n);
    void adjoint(int A[numAttr][numAttr],int adj[numAttr][numAttr]);
    bool inverse(int A[numAttr][numAttr], int inverse[numAttr][numAttr]);
};
