#pragma once

#include "matrix.h"

// 专门设置顶点，预先定义好
class VertexCurveABC {
public:
  // float parameters[3][1] = {{0}, {0}, {0}}; // abc
  Matrix<float, 3, 1> parameters;
};


/**
 * 以逆深度形式存储的顶点
 */
class VertexInverseDepth {
public:
  VertexInverseDepth() {}
  Matrix<float, 1, 1> parameters;
  void Plus(const Matrix<float, 1, 1> &delta) { parameters[0] += delta[0]; }
};

class VertexPose {
public:
  VertexPose() {}
  Matrix<float, 7, 1> parameters;  // xyz xyzw

  void Plus(const Matrix<float, 6, 1> &delta) {
    parameters(0, 0) += delta[0];
    parameters(1, 0) += delta[1];
    parameters(2, 0) += delta[2];

    float q[4][1] = {};
    q[0][0] = parameters[3];
    q[1][0] = parameters[4];
    q[2][0] = parameters[5];
    q[3][0] = parameters[6];
    // Matrix<double,4,1> dq;
    float dq[4][1] = {};
    dq[0][0] = delta[3] / 2;
    dq[1][0] = delta[4] / 2;
    dq[2][0] = delta[5] / 2;
    dq[3][0] = 1;
    Matrix<double, 4, 1> m;
    // std::cout << "q:" << q[0][0] << " " << q[1][0] << " " << q[2][0] << " "
    //           << q[3][0] << " " << std::endl;
    // std::cout << "dq:" << dq[0] << " " << dq[1] << " " << dq[2] << " "
    //           << dq[3] << " " << std::endl;
    MyQuaternionMultiple(q, dq, m.data_);
    m.normalized();
    parameters[3] = m[0];
    parameters[4] = m[1];
    parameters[5] = m[2];
    parameters[6] = m[3];
  }
};
