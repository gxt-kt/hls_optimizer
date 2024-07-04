#pragma once
// #include "ap_fixed.h"

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
#define VAR(...) G_VARHELP(COUNT_ARGS(__VA_ARGS__) __VA_OPT__(, ) __VA_ARGS__)

// constexpr int A_D1 =  4;
// constexpr int A_D2 =  4;
// constexpr int B_D1 =  A_D2;
// constexpr int B_D2 =  4;
// constexpr int C_D1 =  A_D1;
// constexpr int C_D2 =  B_D2;

////using A_T = ap_int<8>;
// using A_T = ap_int<8>;
// using B_T = ap_int<8>;
// using C_T = ap_int<8>;

// void MatrixMultiple(A_T A[A_D1][A_D2], B_T B[B_D1][B_D2], C_T C[C_D1][C_D2]);
// void MatrixMultiple(ap_int<8> A[4][4], ap_int<8> B[4][4], ap_int<16>
// C[4][4]);

template <typename T_a, typename T_b, typename T_c, int x, int y, int z>
void MyMatrixMultiple(T_a A[x][y], T_a B[y][z], T_a C[x][z]) {
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
void MyMatrixAdd(T_a A[x][y], T_b B[z][w]) {
  for (int i = 0; i < z; i++) {
    for (int j = 0; j < w; j++) {
      A[i + index_i][j + index_j] += B[i][j];
    }
  }
}
template <typename T_a, typename T_b, int x, int y, int z, int w, int index_i,
          int index_j>
void MyMatrixSub(T_a A[x][y], T_b B[z][w]) {
  for (int i = 0; i < z; i++) {
    for (int j = 0; j < w; j++) {
      A[i + index_i][j + index_j] -= B[i][j];
    }
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
// template <typename T_a, int x, int y>
// void MyMatrixDebug(T_a A[x][y]) {
//   std::cout << "============MyMatrixDebug============" << std::endl;
//   for (int i = 0; i < x; i++) {
//     for (int j = 0; j < y; j++) {
//       std::cout << A[i][j] << " ";
//     }
//     std::cout << std::endl;
//   }
// }
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

template <typename T_a, int N>
void ldl(T_a A[N][N], double L[N][N], double D[N][N]) {
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
void solve(double L[N][N], double D[N][N], double b[N][1], double x[N][1]) {
  int n = N;
  int i, j;
  double y[N], sum;

  for (i = 0; i < n; i++) {
    sum = b[i][0];
    for (j = 0; j < i; j++) {
      sum -= L[i][j] * y[j];
    }
    y[i] = sum / L[i][i];
    // y[i] = (sum - L[i][i] * y[i]) / L[i][i]; // 修改此行
  }

  for (i = 0; i < n; i++) {
    x[i][0] = y[i] / sqrt(D[i][i]);
  }

  for (i = n - 1; i >= 0; i--) {
    sum = x[i][0];
    for (j = i + 1; j < n; j++) {
      sum -= L[j][i] * x[j][0];
    }
    x[i][0] = sum / L[i][i];
  }
}

template <typename T_a, typename T_x, typename T_b, int N>
void LdltSolve(T_a A[N][N], T_x x[N][1], T_b b[N][1], bool log = false) {
  // T_a A[N][N] = {{4, -2, 2}, {-2, 2, -4}, {2, -4, 11}};
  //   double b[N] = {6, -10, 27};
  // double L[N][N], D[N], x[N];
  T_a L[N][N] = {}, D[N][N] = {};
  int i, j;

  ldl(A, L, D);
  if (log) {
    std::cout << __PRETTY_FUNCTION__ << "log" << std::endl;
    T_a L_transpose[N][N] = {};
    MyMatrixTranspose<double, double, N, N>(L, L_transpose);
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

// template <typename T_a,typename T_b,typename T_c,int k,int n>
// void MyGevm(T_a A[k], T_a B[k][n], T_a C[1][n]) {
//	#pragma HLS ARRAY_RESHAPE variable=B complete dim=1
//	#pragma HLS ARRAY_RESHAPE variable=A complete dim=2
//		for(int j=0;j<n;j++){
//			#pragma HLS PIPELINE II=1
//			C[0][j]=0;
//			// 循环乘四次，并进行相加
//			for(int k_i=0;_i<k;k_i++){
//				C[0][j]=C[0][j]+A[0][k_i]*B[k_i][j];
//			}
//		}
//}
//
// 专门设置顶点，预先定义好
class VertexCurveABC {
public:
  bool enable = false;
  double parameters[3][1] = {{0}, {0}, {0}}; // abc
};

template <int residual_dimension, int num_verticies> class Edge {
public:
  explicit Edge(){};

  ~Edge(){};

  /// 计算残差，由子类实现
  // virtual void ComputeResidual() = 0;

  /// 计算雅可比，由子类实现
  /// 本后端不支持自动求导，需要实现每个子类的雅可比计算方法
  // virtual void ComputeJacobians() = 0;

  //    ///计算该edge对Hession矩阵的影响，由子类实现
  //    virtual void ComputeHessionFactor() = 0;

  /// 计算平方误差，会乘以信息矩阵
  double Chi2() {
    //	  return 1;
    double result[1][residual_dimension];
    MyMatrixMultiple<double, double, double, 1, residual_dimension,
                     residual_dimension>(residual_, information_, result);
    double res[1][1];
    MyMatrixMultiple<double, double, double, 1, residual_dimension, 1>(
        result, residual_, res);
    // std::cout << VAR(residual_[0][0]) << std::endl;
    // std::cout << VAR(information_[0][0]) << std::endl;
    // std::cout << VAR(res[0][0]) << std::endl;
    return res[0][0];
  }

  void
  SetInformation(double information[residual_dimension][residual_dimension]) {
#pragma HLS PIPELINE
#pragma HLS UNROLL
    for (int i = 0; i < residual_dimension; i++) {
      for (int j = 0; j < residual_dimension; j++) {
        information_[i][j] = information[i][j];
      }
    }
  }

public:
  //  unsigned long id_;                          // edge id
  int ordering_id_; // edge id in problem
  //  std::vector<std::shared_ptr<Vertex>> verticies_;  // 该边对应的顶点
  //  double residual_;
  double residual_[1][residual_dimension] = {0}; // 残差
                                                 //  std::vector<MatXX>
  //      jacobians_;  // 雅可比，每个雅可比维度是 residual x vertex[i]
  //  MatXX information_;      // 信息矩阵
  double information_[residual_dimension][residual_dimension] = {{0}};
};

// template <int residual_dimension, int num_verticies>
class CurveFittingEdge : public Edge<1, 1> {
public:
  static const int residual_dimension = 1;
  static const int num_verticies = 1;
  double x_ = 0, y_ = 0;
  void ComputeResidual() {
    // GetResidual();
    // residual_[0][0] = 10;
    // std::cout <<
    // VAR(verticies_[0].parameters[0][0],verticies_[0].parameters[1][0],verticies_[0].parameters[2][0])
    // << std::endl;
    residual_[0][0] = std::exp(verticies_[0].parameters[0][0] * x_ * x_ +
                               verticies_[0].parameters[1][0] * x_ +
                               verticies_[0].parameters[2][0]) -
                      y_;
    // std::cout << VAR(residual_[0][0]) << std::endl;
  }
  void ComputeJacobians() {
    double exp_x = std::exp(verticies_[0].parameters[0][0] * x_ * x_ +
                            verticies_[0].parameters[1][0] * x_ +
                            verticies_[0].parameters[2][0]);
    // std::cout << VAR(exp_y) << std::endl;
    jacobians_0[0][0] = x_ * x_ * exp_x;
    jacobians_0[0][1] = x_ * exp_x;
    jacobians_0[0][2] = exp_x;
    MATRIXDEBUG(jacobians_0);
  }

  double jacobians_0[1][3];
  // const int residual_dimension=1;
  // const int num_verticies=1;
  VertexCurveABC verticies_[num_verticies];
  const int num_verticies_ = num_verticies;
};

class Problem {
public:
  Problem() {}
  // 最大支持200条边
  static const int max_edge_curves_size = 200;
  int edge_curves_size = 0;
  CurveFittingEdge edge_curves[max_edge_curves_size];

  void CalculateResidual() {
    for (int i = 0; i < max_edge_curves_size; i++) {
      if (i >= edge_curves_size)
        continue;
      edge_curves[i].ComputeResidual();
    }
  }
  double Chi2() {
    double res;
    for (int i = 0; i < max_edge_curves_size; i++) {
      // std::cout << VAR(res) << std::endl;
      if (i >= edge_curves_size)
        continue;
      res += edge_curves[i].Chi2();
      // res += edge_curves[i].Chi2();
      // std::cout << VAR(res) << std::endl;
    }
    return res;
  }
  void MakeHessian() {
    for (int i = 0; i < max_edge_curves_size; i++) {
      if (i >= edge_curves_size)
        continue;
      edge_curves[i].ComputeResidual();
      edge_curves[i].ComputeJacobians();

      // 遍历所有顶点
      for (int j = 0; j < edge_curves[i].num_verticies_; j++) {
      }

      double JtW[3][1];
      double transpose[3][1];
      MyMatrixTranspose<double, double, 1, 3>(edge_curves[i].jacobians_0,
                                              transpose);
      // MATRIXDEBUG(edge_curves[i].jacobians_0);
      MyMatrixMultiple<double, double, double, 3, 1, 1>(
          transpose, edge_curves[i].information_, JtW);

      // double jacobian_j[1][3] = jacobians[j];
      double hessian[3][3];
      MyMatrixMultiple<double, double, double, 3, 1, 3>(
          JtW, edge_curves[i].jacobians_0, hessian);
      MyMatrixAdd<double, double, 3, 3, 3, 3, 0, 0>(hessian_, hessian);

      double bb[3][1];
      MyMatrixMultiple<double, double, double, 3, 1, 1>(
          JtW, edge_curves[i].residual_, bb);
      MyMatrixSub<double, double, 3, 1, 3, 1, 0, 0>(b_, bb);
      // MATRIXDEBUG(hessian_);
      // MATRIXDEBUG(edge_curves[i].residual_);
      // = jacobian_i.transpose() * edge.second->Information(); //3*1 * 1*1
    }
    MATRIXDEBUG(hessian_);
    MATRIXDEBUG(b_);
  }
  void SolveLinearSystem() {
    // TODO:
    // 在这里我们假设已经求解完成了
    // 输入是hessian_ 和 b_
    // 输出是 delta_x_
    // delta_x_ = Hessian_.inverse() * b_;
    //
    LdltSolve<double, double, double, 3>(hessian_, delta_x_, b_, true);
    MATRIXDEBUG(delta_x_);
    double h_d[3][1] = {};
    MyMatrixMultiple<double, double, double, 3, 3, 1>(hessian_, delta_x_, h_d);
    MATRIXDEBUG(delta_x_);
    MATRIXDEBUG(h_d);
  }
  void UpdateStates() {}
  void Solve(int it_cnts = 10) {
    MakeHessian();
    SolveLinearSystem();
    std::terminate();
    int cnt = 0;
    int false_cnt = 0;
    bool stop = false;
    for (int i = 0; cnt < it_cnts; i++) {
      // 优化退出条件1： delta_x_ 很小则退出
      if (MySquaredNorm<double, 3>(delta_x_) <= 1e-6) {
        std::cout << "stop!: SquaredNorm(delta_x)<=number" << std::endl;
        stop = true;
        break;
      }
      if (false_cnt > 10) {
        std::cout << "stop!: false_cnt > 10" << std::endl;
        stop = true;
        break;
      }
      UpdateStates();
    }
  }
  double hessian_[3][3] = {};
  double b_[3][1];
  double chi2_;
  double delta_x_[3][1];
};
