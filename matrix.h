#pragma once
// #include "ap_fixed.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>
#include <iostream>

#include "help.h"

template <typename T = float, size_t X = 1, size_t Y = 1> class Matrix {
public:
  using type = T;
  Matrix() {}
  type &operator()(size_t i, size_t j) { return data_[i][j]; }
  const type &operator()(size_t i, size_t j) const { return data_[i][j]; }
  type &operator[](size_t i) { return data_[i][0]; }
  const type &operator[](size_t i) const { return data_[i][0]; }
  type data_[X][Y] = {};

public:
  Matrix<type, Y, X> transpose() {
    Matrix<type, Y, X> res;
    MyMatrixTranspose<type, type, X, Y>(this->data_, res.data_);
    return res;
  }
  template <typename U, size_t M, size_t N>
  Matrix<T, X, N> operator*(const Matrix<U, M, N> &matrix) {
    Matrix<T, X, N> res;
    static_assert(M == Y);
    MyMatrixMultiple<T, U, T, X, Y, N>(this->data_, matrix.data_, res.data_);
    return res;
  }
  static Matrix<T, X, Y> Identity() {
    static_assert(X == Y);
    Matrix<T, X, Y> res;
    for (int i = 0; i < X; i++) {
      res(i, i) = 1;
    }
  }
  template <typename U> Matrix<T, X, Y> operator*(U value) {
    Matrix<T, X, Y> res;
    for (int i = 0; i < X; i++) {
      for (int j = 0; j < Y; j++) {
        res(i, j) = data_[i][j] * value;
      }
    }
    return res;
  }
  // 一定注意相加和相减都是原地操作
  void operator+(const Matrix<T, X, Y> &matrix) {
    MyMatrixAdd<T, T, X, Y, X, Y, 0, 0>(this->data_, matrix.data_);
  }
  void operator-(const Matrix<T, X, Y> &matrix) {
    MyMatrixSub<T, T, X, Y, X, Y, 0, 0>(this->data_, matrix.data_);
  }
  template <typename U> void SetValue(U value) {
    for (int i = 0; i < X; i++) {
      for (int j = 0; j < Y; j++) {
        data_[i][j] = value;
      }
    }
  }
  void normalized() {
    for (size_t j = 0; j < Y; j++) {
      // 小心这里norm相加会溢出
      T norm = 0;
      for (size_t i = 0; i < X; i++) {
        norm += data_[i][j] * data_[i][j];
      }
      norm = sqrt(norm);
      for (size_t i = 0; i < X; i++) {
        data_[i][j] /= norm;
      }
    }
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const Matrix<T, X, Y> &matrix) {
    MATRIXDEBUG(matrix.data_);
    return os;
  }
};

template <typename T, size_t N> using Vector = Matrix<T, N, 1>;




void MatrixMultiple(unsigned int A[100], unsigned int B[100],
                    unsigned int C[100], unsigned int A_a[100],
                    unsigned int B_b[100]);
