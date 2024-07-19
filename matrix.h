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
  type &operator()(size_t i) { return data_[i][0]; }
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
  template <typename U> void AddDiagonal(U value) {
    static_assert(X == Y);
    for (int i = 0; i < X; i++) {
      (*this)(i, i) += value;
    }
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
  template <size_t N> void SetleftCols(const Matrix<T, X, N> &matrix) {
    for (int i = 0; i < X; i++) {
      for (int j = 0; j < N; j++) {
        (*this)(i, j) = matrix(i, j);
      }
    }
  }
  template <size_t N> void SetrightCols(const Matrix<T, X, N> &matrix) {
    static_assert(N >= 0 && Y >= N);
    for (int i = 0; i < X; i++) {
      for (int j = Y - N; j < Y; j++) {
        (*this)(i, j) = matrix(i, j - (Y - N));
      }
    }
  }
  T MaxDiagonal() {
    static_assert(X == Y);
    auto res = MyMatrixAMax<T, T, X>(this->data_);
    return res;
  }

  Matrix<T, 3, 3> hat() {
    static_assert(X == 3 && Y == 1);
    Matrix<T, 3, 3> Omega;
    Omega(0, 0) = 0.0;
    Omega(0, 1) = -(*this)(2);
    Omega(0, 2) = (*this)(1);
    Omega(1, 0) = (*this)(2);
    Omega(1, 1) = 0.0;
    Omega(1, 2) = -(*this)(0);
    Omega(2, 0) = -(*this)(1);
    Omega(2, 1) = (*this)(0);
    Omega(2, 2) = 0.0;
    // Omega << 0.0, -omega(2), omega(1),
    //          omega(2), 0.0, -omega(0),
    //         -omega(1), omega(0), 0.0;
    return Omega;
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
      // NOTE: 小心这里norm相加会溢出
      float norm = 0;
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
  template <size_t N> Matrix<type, N, 1> head() {
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
    // Matrix<T_t, 3, 1> res;
    // QuaternionMultiply<T_q, T_t>(this->data_, t.data_, res.data_);
    // return res;
    return toRotationMatrix() * t;
  }
  Matrix<T_q, 3, 3> toRotationMatrix() {
    T_q q0 = w();
    T_q q1 = x();
    T_q q2 = y();
    T_q q3 = z();

    Matrix<T_q, 3, 3> R;
    R(0, 0) = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
    R(0, 1) = 2 * (q1 * q2 - q0 * q3);
    R(0, 2) = 2 * (q1 * q3 + q0 * q2);

    R(1, 0) = 2 * (q1 * q2 + q0 * q3);
    R(1, 1) = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
    R(1, 2) = 2 * (q2 * q3 - q0 * q1);

    R(2, 0) = 2 * (q1 * q3 - q0 * q2);
    R(2, 1) = 2 * (q2 * q3 + q0 * q1);
    R(2, 2) = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;

    return R;
  }
  Quaternion<T_q> inverse() {
    //    T_q norm = x() * x() + y() * y() + z() * z() + w() * w();
    //    if (norm < std::numeric_limits<T_q>::epsilon()) {
    //      // 四元数范数接近 0，无法计算逆
    //      return *this;
    //    }

    //    T_q inv_norm = static_cast<T_q>(1.0) / norm;
    //    Quaternion<T_q> inv;
    //    inv.x() = -x() * inv_norm;
    //    inv.y() = -y() * inv_norm;
    //    inv.z() = -z() * inv_norm;
    //    inv.w() = w() * inv_norm;
    Quaternion<T_q> inv;
    inv.x() = -x();
    inv.y() = -y();
    inv.z() = -z();
    inv.w() = w();
    return inv;
  }
  // eigen3中也是按照wxyz存储的
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
