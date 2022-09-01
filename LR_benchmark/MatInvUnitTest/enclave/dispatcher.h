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
#define numAttr 64

using namespace std;

class ecall_dispatcher
{
  public:
    ecall_dispatcher();
    int MatInv();
    void close();
  private:
    // Adjoint method
    void getCofactor(int A[numAttr][numAttr], int temp[numAttr][numAttr], int p, int q, int n);
    int determinant(int A[numAttr][numAttr], int n);
    void adjoint(int A[numAttr][numAttr],int adj[numAttr][numAttr]);
    bool inverse(int A[numAttr][numAttr], int inverse[numAttr][numAttr]);

    // LU decomposition
    // void LUP_Descomposition(int A[numAttr][numAttr], int L[numAttr][numAttr], int U[numAttr][numAttr], int P[numAttr]);
    // void LUP_Solve(int L[numAttr][numAttr], int U[numAttr][numAttr], int P[numAttr], double b[numAttr], int x[numAttr]);
    // double * LUP_solve_inverse(double A[numAttr*numAttr]);
    void LUP_Descomposition(double A[numAttr*numAttr],double L[numAttr*numAttr],double U[numAttr*numAttr],int P[numAttr]);
    double* LUP_Solve(double L[numAttr*numAttr],double U[numAttr*numAttr],int P[numAttr],double b[numAttr]);
    int getNext(int i, int m, int n);
    int getPre(int i, int m, int n);
    void movedata(double *mtx, int i, int m, int n);
    void transpose(double *mtx, int m, int n);
    double* LUP_solve_inverse(double A[numAttr*numAttr]);
    double* mul(double A[numAttr*numAttr],double B[numAttr*numAttr]);
};
