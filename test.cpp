#include "problem.h"
#include <exception>
#include <iostream>

float my_exp(float x) {
  float result = 1.0;
  float term = 1.0;
  int i = 1;

  while (fabs(term) > 1e-10) {
    term *= x / i;
    result += term;
    i++;
  }

  return result;
}

int main() {
  {
    VertexPose tmp;
    tmp.parameters[0] = 1;
    tmp.parameters[1] = 2;
    tmp.parameters[2] = 3;
    tmp.parameters[3] = 0.1;
    tmp.parameters[4] = 0.2;
    tmp.parameters[5] = 0.3;
    tmp.parameters[6] = sqrt(1 - (0.01 + 0.04 + 0.09));
    Matrix<float, 6, 1> delta;
    delta[0] = 1;
    delta[1] = 2;
    delta[2] = 3;
    delta[3] = 0.1;
    delta[4] = 0.2;
    delta[5] = 0.1;
    // delta[6] = sqrt(1 - (0.01 + 0.04 + 0.01));
    tmp.Plus(delta);
    std::cout << tmp.parameters << std::endl;
    // std::terminate();
  }

  {
    float A[100] = {};
    float B[100] = {};
    float C[100] = {};
    float A_a[100] = {};
    float B_b[100] = {};
    float a = 1, b = 2, c = 1;
    for (int i = 0; i < 100; i++) {
      A[i] = i / 100.0;
      float x = A[i];
      B[i] = my_exp(a * x * x + b * x + c);
    }
    unsigned int A_[100] = {};
    unsigned int B_[100] = {};
    unsigned int C_[100] = {};
    unsigned int A_a_[100] = {};
    unsigned int B_b_[100] = {};
    memcpy(A_, A, 100 * sizeof(float));
    memcpy(B_, B, 100 * sizeof(float));
    memcpy(C_, C, 100 * sizeof(float));
    memcpy(A_a_, A_a, 100 * sizeof(float));
    memcpy(B_b_, B_b, 100 * sizeof(float));
    MatrixMultiple(A_, B_, C_, A_a_, B_b_);
    memcpy(C, C_, 100 * sizeof(float));
    memcpy(A_a, A_a_, 100 * sizeof(float));
    memcpy(B_b, B_b_, 100 * sizeof(float));
    for (int i = 0; i < 3; i++) {
      std::cout << C[i] << std::endl;
    }
    for (int i = 0; i < 3; i++) {
      std::cout << A[i] << std::endl;
      std::cout << A_a[i] << std::endl;
      std::cout << B[i] << std::endl;
      std::cout << B_b[i] << std::endl;
    }
  }
}
