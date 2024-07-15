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
  template <typename U> Matrix<T, X, Y> operator/(U value) {
    Matrix<T, X, Y> res;
    for (int i = 0; i < X; i++) {
      for (int j = 0; j < Y; j++) {
        res(i, j) = data_[i][j] / value;
      }
    }
    return res;
  }
  // // 一定注意相加和相减都是原地操作
  // void operator+(const Matrix<T, X, Y> &matrix) {
  //   MyMatrixAdd<T, T, X, Y, X, Y, 0, 0>(this->data_, matrix.data_);
  // }
  Matrix<T, X, Y> operator+(const Matrix<T, X, Y> &matrix) {
    Matrix<T, X, Y> res;
    for (int i = 0; i < X; i++) {
      for (int j = 0; j < Y; j++) {
        res(i, j) = (*this)(i, j) + matrix(i, j);
      }
    }
    return res;
  }
  // void operator-(const Matrix<T, X, Y> &matrix) {
  //   MyMatrixSub<T, T, X, Y, X, Y, 0, 0>(this->data_, matrix.data_);
  // }
  Matrix<T, X, Y> operator-(const Matrix<T, X, Y> &matrix) {
    Matrix<T, X, Y> res;
    for (int i = 0; i < X; i++) {
      for (int j = 0; j < Y; j++) {
        res(i, j) = (*this)(i, j) - matrix(i, j);
      }
    }
    return res;
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
  type &x() {
    static_assert(X == 3 && Y == 1, "only X==3&&Y==1 canuse it");
    return (*this)[0];
  }
  type &y() {
    static_assert(X == 3 && Y == 1, "only X==3&&Y==1 canuse it");
    return (*this)[1];
  }
  type &z() {
    static_assert(X == 3 && Y == 1, "only X==3&&Y==1 canuse it");
    return (*this)[2];
  }
  template <size_t N> Matrix<type, N, 1> &head() {
    static_assert(Y == 1 && N <= X && N > 0, "must Y == 1 && N<=X && N>0");
    Matrix<type, N, 1> res;
    for (int i = 0; i < N; i++) {
      res[i] = (*this)[i];
    }
    return res;
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const Matrix<T, X, Y> &matrix) {
    MATRIXDEBUG(matrix.data_);
    return os;
  }
};

template <typename T_q> class Quaternion : public Matrix<T_q, 4, 1> {
public:
  template <typename T_t>
  Matrix<T_t, 3, 1> operator*(const Matrix<T_t, 3, 1> &t) {
    Matrix<T_t, 3, 1> res;
    QuaternionMultiply<T_q, T_t>(this->data_, t.data_, res.data_);
    return res;
  }
  Quaternion<T_q> inverse() {
    T_q norm = x() * x() + y() * y() + z() * z() + w() * w();
    if (norm < std::numeric_limits<T_q>::epsilon()) {
      // 四元数范数接近 0，无法计算逆
      return *this;
    }

    T_q inv_norm = static_cast<T_q>(1.0) / norm;
    Quaternion<T_q> inv;
    inv.x() = -x() * inv_norm;
    inv.y() = -y() * inv_norm;
    inv.z() = -z() * inv_norm;
    inv.w() = w() * inv_norm;
    return inv;
  }
  T_q &x() { return (*this)[1]; }
  T_q &y() { return (*this)[2]; }
  T_q &z() { return (*this)[3]; }
  T_q &w() { return (*this)[0]; }
  const T_q &x() const { return (*this)[1]; }
  const T_q &y() const { return (*this)[2]; }
  const T_q &z() const { return (*this)[3]; }
  const T_q &w() const { return (*this)[0]; }
};

template <typename T, size_t N> using Vector = Matrix<T, N, 1>;

void MatrixMultiple(unsigned int A[100], unsigned int B[100],
                    unsigned int C[100], unsigned int A_a[100],
                    unsigned int B_b[100]);
