#ifndef _PTMATRIX_H
#define _PTMATRIX_H

#include <vector>
#include <iostream>
#include <functional>

using namespace std;

template<class T>
class PtMatrix;

template <typename T>
ostream& operator<<(ostream &, const PtMatrix<T> &);

template<class T>
class PtMatrix {
  T dummy;
  public:
    PtMatrix() : dummy(T()) {
      transpose = false;
    }
    PtMatrix(const T &dummy) : dummy(dummy) {
      transpose = false;
    }
    
    PtMatrix(unsigned nRows, unsigned nCols, const T &dummy) : dummy(dummy) {
      Resize(nRows, nCols);
    }
    
    PtMatrix(unsigned nRows, unsigned nCols) : dummy(T()) {
      Resize(nRows, nCols);
    }
    
    PtMatrix operator+(const PtMatrix &other) const;
    PtMatrix &operator+=(const PtMatrix &other);
    
    PtMatrix operator-(const PtMatrix &other) const;
    PtMatrix &operator-=(const PtMatrix &other);
    
    PtMatrix operator*(PtMatrix &other) const;
    PtMatrix operator*(vector<T> &other) const;
    PtMatrix &operator*=(PtMatrix &other);
    PtMatrix &operator*=(vector<T> &other);
    PtMatrix &operator*=(T &other);
    
    PtMatrix &operator=(const PtMatrix &other);
    
    T &operator()(unsigned row, unsigned col);
    T const &operator()(unsigned row, unsigned col) const;
    
    vector<T> &operator[](unsigned row);
    vector<T> const &operator[](unsigned row) const;
    
    void MultByTranspose();
    void Transpose();
    void Invert(T& det, std::function<void(T&)> reduce = NULL);
    
    void Determinant(T& det, std::function<void(T&)> reduce = NULL) const;
    
    unsigned NumRows() const;
    unsigned NumCols() const;
    
    void Resize(unsigned nRows, unsigned nCols);
    
    void AddRow(vector<T> &row);
    void Concatenate(PtMatrix<T> &other);
    void Clear();
    
    void MapAll(std::function<void(T&)> func);
    
    friend ostream &operator<< <>(ostream &os, const PtMatrix &m);
  private:
    vector<vector<T>> mat;
    bool transpose;
    
    T &ElemAt(unsigned row, unsigned col);
    T const &ElemAt(unsigned row, unsigned col) const;
    
    void Determinant(T& det, vector<bool> &usedRows,
                     vector<bool> &usedCols, unsigned dim,
                     std::function<void(T&)> reduce = NULL) const;
};

#endif
