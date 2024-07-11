#pragma once


#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>

#define COUNT_ARGS_IMPL(_null, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N,     \
                        ...)                                                   \
  N
#define COUNT_ARGS(...)                                                        \
  COUNT_ARGS_IMPL(0 __VA_OPT__(, ) __VA_ARGS__, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, \
                  0)

#define VARIMP(x) #x << "=" << x

#define G_VAR0() ""
#define G_VAR1(_1) VARIMP(_1)
#define G_VAR2(_1, _2) VARIMP(_1) << "," << G_VAR1(_2)
#define G_VAR3(_1, _2, _3) VARIMP(_1) << "," << G_VAR2(_2, _3)
#define G_VAR4(_1, _2, _3, _4) VARIMP(_1) << "," << G_VAR3(_2, _3, _4)
#define G_VAR5(_1, _2, _3, _4, _5) VARIMP(_1) << "," << G_VAR4(_2, _3, _4, _5)
#define G_VAR6(_1, _2, _3, _4, _5, _6)                                         \
  VARIMP(_1) << "," << G_VAR5(_2, _3, _4, _5, _6)
#define G_VAR7(_1, _2, _3, _4, _5, _6, _7)                                     \
  VARIMP(_1) << "," << G_VAR6(_2, _3, _4, _5, _6, _7)
#define G_VAR8(_1, _2, _3, _4, _5, _6, _7, _8)                                 \
  VARIMP(_1) << "," << G_VAR7(_2, _3, _4, _5, _6, _7, _8)
#define G_VAR9(_1, _2, _3, _4, _5, _6, _7, _8, _9)                             \
  VARIMP(_1) << "," << G_VAR8(_2, _3, _4, _5, _6, _7, _8, _9)
#define G_VAR10(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10)                       \
  VARIMP(_1) << "," << G_VAR9(_2, _3, _4, _5, _6, _7, _8, _9, _10)

#define G_VARHELPIMP(N, ...) G_VAR##N(__VA_ARGS__)
#define G_VARHELP(N, ...) G_VARHELPIMP(N __VA_OPT__(, ) __VA_ARGS__)

// Usage: gDebug() << VAR(a,b) // stdout: a = ${a} , b = ${b}
// #define VAR(...) G_VARHELP(COUNT_ARGS(__VA_ARGS__) __VA_OPT__(, )
// __VA_ARGS__)
#define VAR(...) ""


template <typename T_a, typename T_b, typename T_c, int x, int y, int z>
void MyMatrixMultiple(const T_a A[x][y],const T_b B[y][z], T_c C[x][z]) {
#pragma HLS ARRAY_RESHAPE variable = B complete dim = 1
#pragma HLS ARRAY_RESHAPE variable = A complete dim = 2
  for (int i = 0; i < x; i++)
    for (int j = 0; j < z; j++) {
#pragma HLS PIPELINE II = 1
      C[i][j] = 0;
      // 循环乘四次，并进行相加
      for (int k = 0; k < y; k++) {
        C[i][j] = C[i][j] + A[i][k] * B[k][j];
      }
    }
}
template <typename T_a, typename T_b, int x, int y>
void MyMatrixMultipleNumber(T_a A[x][y], T_b val, T_a C[x][y]) {
  for (int i = 0; i < x; i++)
    for (int j = 0; j < y; j++) {
      C[i][j] = A[i][j] * val;
    }
}
template <typename T_a, typename T_b, int x, int y>
void MyMatrixTranspose(T_a A[x][y], T_a B[y][x]) {
// #pragma HLS ARRAY_RESHAPE variable = A complete dim = 2
// #pragma HLS ARRAY_RESHAPE variable = B complete dim = 2
#pragma HLS PIPELINE II = 1
  for (int i = 0; i < x; i++) {
#pragma HLS PIPELINE II = 1
    for (int j = 0; j < y; j++) {
      B[j][i] = A[i][j];
    }
  }
}
template <typename T_a, typename T_b, int x, int y, int z, int w, int index_i,
          int index_j>
void MyMatrixAdd(T_a A[x][y],const T_b B[z][w]) {
  for (int i = 0; i < z; i++) {
    for (int j = 0; j < w; j++) {
      A[i + index_i][j + index_j] += B[i][j];
    }
  }
}
template <typename T_a, typename T_b, int x, int y, int z, int w, int index_i,
          int index_j>
void MyMatrixSub(T_a A[x][y],const T_b B[z][w]) {
  for (int i = 0; i < z; i++) {
    for (int j = 0; j < w; j++) {
      A[i + index_i][j + index_j] -= B[i][j];
    }
  }
}
template <typename T_a, typename T_b, int x, int y, int index_s, int cnts>
void MyMatrixAddNumber(T_a A[x][y], T_b val) {
  for (int i = index_s; i < index_s + cnts; i++) {
    A[i][i] += val;
  }
}

template <typename T_a, int x> double MySquaredNorm(T_a A[x][1]) {
  double norm = 0.0;
#pragma HLS PIPELINE II = 1
  for (int i = 0; i < x; i++) {
    norm += A[i][0] * A[i][0];
  }
  return norm;
}
template <typename T_a, typename T_res, int x>
double MyMatrixAMax(T_a A[x][x]) {
  T_a res = std::abs(A[0][0]);
  for (int i = 1; i < x; i++) {
    res = std::max(res, std::abs(A[i][i]));
  }
  return res;
}
template <typename T_a, typename T_b, int x, int y>
void MyMatrixSet(T_a A[x][y], T_b val) {
  for (int i = 0; i < x; i++) {
    for (int j = 0; j < y; j++) {
      A[i][j] = val;
    }
  }
}

template <typename T_a, int x, int y>
void MyMatrixDebug(T_a (&A)[x][y], std::string str = "") {
  std::string str_out = "============" + str + "============";
  std::cout << str_out << std::endl;
  for (int i = 0; i < x; i++) {
    for (int j = 0; j < y; j++) {
      std::cout << A[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::string(str_out.size(), '=') << std::endl;
}
#define MATRIXDEBUG(m) MyMatrixDebug(m, #m)
// #define MATRIXDEBUG(m)

template <typename T_a, int N>
void ldl(T_a A[N][N], float L[N][N], float D[N][N]) {
  int n = N;
  int i, j, k;
  double sum;

  for (i = 0; i < n; i++) {
    for (j = 0; j <= i; j++) {
      sum = A[i][j];
      for (k = 0; k < j; k++) {
        sum -= L[i][k] * D[k][k] * L[j][k];
      }
      if (i == j) {
        D[i][i] = sum;
        L[i][i] = 1;
      } else {
        L[i][j] = sum / D[j][j];
      }
    }
    for (j = i + 1; j < n; j++) {
      sum = A[j][i];
      for (k = 0; k < i; k++) {
        sum -= L[j][k] * D[k][k] * L[i][k];
      }
      L[j][i] = sum / D[i][i];
    }
  }
}

template <int N>
void solve(float L[N][N], float D[N][N], float b[N][1], float x[N][1]) {
  int n = N;
  int i, j;
  double y[N][1] = {}, sum = 0;

  for (i = 0; i < n; i++) {
    sum = b[i][0];
    for (j = 0; j < i; j++) {
      sum -= L[i][j] * y[j][0];
    }
    y[i][0] = sum / L[i][i];
  }

  // MATRIXDEBUG(y);

  double z[N][1] = {};
  for (i = 0; i < n; i++) {
    z[i][0] = y[i][0] / D[i][i];
  }

  // MATRIXDEBUG(z);

  for (i = n - 1; i >= 0; i--) {
    sum = z[i][0];
    for (j = i + 1; j < n; j++) {
      sum -= L[j][i] * x[j][0];
    }
    x[i][0] = sum;
  }
}

template <typename T_a, typename T_x, typename T_b, int N>
void LdltSolve(T_a A[N][N], T_x x[N][1], T_b b[N][1], bool log = false) {
  T_a L[N][N] = {}, D[N][N] = {};
  int i, j;

  ldl(A, L, D);
 if (log) {
   std::cout << __PRETTY_FUNCTION__ << "log" << std::endl;
   T_a L_transpose[N][N] = {};
   MyMatrixTranspose<float, float, N, N>(L, L_transpose);
   MATRIXDEBUG(L);
   MATRIXDEBUG(L_transpose);
   MATRIXDEBUG(D);
   T_a L_D[N][N] = {};
   MyMatrixMultiple<T_a, T_a, T_a, N, N, N>(L, D, L_D);
   T_a L_D_LT[N][N] = {};
   MyMatrixMultiple<T_a, T_a, T_a, N, N, N>(L_D, L_transpose, L_D_LT);
   std::cout << "L_D_LT的结果应该就是hessian矩阵" << std::endl;
   MATRIXDEBUG(L_D_LT);
 }
  solve(L, D, b, x);
}

