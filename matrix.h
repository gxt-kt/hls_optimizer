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
void MyMatrixMultiple(T_a A[x][y], T_b B[y][z], T_c C[x][z]) {
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
    // MATRIXDEBUG(jacobians_0);
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

  void CalculateResidual() {
    for (int i = 0; i < max_edge_curves_size; i++) {
      if (i >= edge_curves_size)
        continue;
      edge_curves[i].ComputeResidual();
    }
  }
  double Chi2() {
    CalculateResidual();
    double res = 0;
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
    LdltSolve<double, double, double, 3>(hessian_, delta_x_, b_, false);
    double h_dx[3][1] = {};
    MyMatrixMultiple<double, double, double, 3, 3, 1>(hessian_, delta_x_, h_dx);
    // MATRIXDEBUG(delta_x_);
    // MATRIXDEBUG(h_dx);
  }
  void UpdateStates() {
    // std::cout << __PRETTY_FUNCTION__ << "begin" << std::endl;
    // MATRIXDEBUG(edge_curves[0].verticies_[0].parameters);
    MyMatrixAdd<double, double, 3, 1, 3, 1, 0, 0>(
        edge_curves[0].verticies_[0].parameters, delta_x_);
    // MATRIXDEBUG(edge_curves[0].verticies_[0].parameters);
    // std::cout << __PRETTY_FUNCTION__ << "end" << std::endl;
  }
  void RollbackStates() {
    // std::cout << __PRETTY_FUNCTION__ << "begin" << std::endl;
    // MATRIXDEBUG(edge_curves[0].verticies_[0].parameters);
    // MATRIXDEBUG(delta_x_);
    MyMatrixSub<double, double, 3, 1, 3, 1, 0, 0>(
        edge_curves[0].verticies_[0].parameters, delta_x_);
    // MATRIXDEBUG(edge_curves[0].verticies_[0].parameters);
    // MATRIXDEBUG(delta_x_);
    // std::cout << __PRETTY_FUNCTION__ << "end" << std::endl;
  }

  /// Hessian 对角线加上或者减去  Lambda
  void AddLambdatoHessianLM() {
    MyMatrixAddNumber<double, double, 3, 3, 0, 3>(hessian_, currentLambda_);
    // unsigned int size = Hessian_.cols();
    // assert(Hessian_.rows() == Hessian_.cols() && "Hessian is not square");
    // for (unsigned long i = 0; i < size; ++i) {
    //   Hessian_(i, i) += my_type{currentLambda_};
    // }
  }
  void RemoveLambdatoHessianLM() {
    MyMatrixAddNumber<double, double, 3, 3, 0, 3>(hessian_, -currentLambda_);
    // unsigned int size = Hessian_.cols();
    // assert(Hessian_.rows() == Hessian_.cols() && "Hessian is not square");
    // for (unsigned long i = 0; i < size; ++i) {
    //   Hessian_(i, i) += my_type{currentLambda_};
    // }
  }
  void Solve(int it_cnts = 10) {
    MakeHessian();
    ComputeLambdaInitLM();
    bool stop = false;
    int iter = 0;
    while (!stop && iter < it_cnts) {
      std::cout << "iter: " << iter << " , chi= " << currentChi_
                << " , Lambda= " << currentLambda_ << std::endl;
      bool one_step_success=false;
      int false_cnt = 0;
      while (!one_step_success) { // 不断尝试 Lambda, 直到成功迭代一步
        // 更新Hx=b为(H+uI)x=b也就是H变为H+uI
        AddLambdatoHessianLM();
        // 解线性方程 Hx=b
        SolveLinearSystem();
        // 把H+uI恢复到原来的H
        RemoveLambdatoHessianLM();
        // 优化退出条件1： delta_x_ 很小则退出
        if (MySquaredNorm<double, 3>(delta_x_) <= 1e-6) {
          std::cout << "stop!: SquaredNorm(delta_x)<=number" << std::endl;
          stop = true;
          break;
        }
        // 优化退出条件2： 尝试很多次都不行则退出
        if (false_cnt > 10) {
          std::cout << "stop!: false_cnt > 10" << std::endl;
          stop = true;
          break;
        }
        // 更新状态量 X = X+ delta_x
        UpdateStates();
        // 判断当前步是否可行以及 LM 的 lambda 怎么更新
        bool isgood = IsGoodStepInLM();
        std::cout << "good??:" << isgood << std::endl;
        if (isgood) {
          // 为下一次优化构建 hessian
          MakeHessian();
          false_cnt = 0;
        } else {
          false_cnt++;
          // 误差没下降，回滚
          RollbackStates();
        }
      }
      iter++;
      // 优化退出条件3： 残差下降满足阈值
      if (std::sqrt(currentChi_) <= stopThresholdLM_) {
        stop = true;
      }
    }
  }
  /// LM 算法中用于判断 Lambda 在上次迭代中是否可以，以及Lambda怎么缩放
  bool IsGoodStepInLM() {
    double scale = 0;
    // scale = delta_x_.transpose() * (my_type{currentLambda_} * delta_x_ + b_);
    double dot[1][1] = {};
    double delta_x_transpose[1][3] = {};
    MyMatrixTranspose<double, double, 3, 1>(delta_x_, delta_x_transpose);
    double delta_x_lambda_[3][1] = {};
    MyMatrixMultipleNumber<double, double, 3, 1>(delta_x_, currentLambda_,
                                                 delta_x_lambda_);
    MyMatrixAdd<double, double, 3, 1, 3, 1, 0, 0>(delta_x_lambda_, b_);
    MyMatrixMultiple<double, double, double, 1, 3, 1>(delta_x_transpose,
                                                      delta_x_lambda_, dot);
    scale = dot[0][0];
    // scale = static_cast<double>((delta_x_.transpose() *
    //                          (my_type{currentLambda_} * delta_x_ +
    //                          b_))(0, 0));
    // my_type scale_tmp = delta_x_.transpose() *
    // (my_type{currentLambda_} * delta_x_ + b_); gDebugCol3(delta_x_);
    // gDebugCol3(currentLambda_);
    // gDebugCol3(b_);
    // gDebugCol3(scale);
    // gDebugCol3(scale_tmp);
    // gDebugCol4() << G_SPLIT_LINE;
    // gDebugCol4(my_type{currentLambda_} * delta_x_ + b_);
    // gDebugCol4(delta_x_.transpose());
    // gDebugCol4(delta_x_.transpose() *
    //            (my_type{currentLambda_} * delta_x_ + b_));
    scale += 1e-3; // make sure it's non-zero :)

    std::cout << VAR(scale) << std::endl;
    // recompute residuals after update state
    // 统计所有的残差
    double tempChi = Chi2();
    // std::cout << "tempChi:" << tempChi << " "
    //           << "currentChi_:" << currentChi_ << std::endl;
    std::cout << VAR(currentChi_, tempChi) << std::endl;
    // for (auto edge : edges_) {
    //   edge.second->ComputeResidual();
    //   tempChi += edge.second->Chi2();
    // }

    // gDebugCol5(tempChi);
    // gDebugCol5(currentChi_);

    double rho = (currentChi_ - tempChi) / scale;
    std::cout << VAR(rho) << std::endl;
    // gDebugCol5(rho);

    // std::terminate();

    if (rho > 0 && std::isfinite(tempChi)) { // last step was good, 误差在下降
      double alpha = 1. - std::pow((2 * rho - 1), 3);
      alpha = std::min(alpha, 2.0 / 3.0);
      double scaleFactor = std::max(1.0 / 3.0, alpha);
      currentLambda_ *= scaleFactor;
      ni_ = 2;
      currentChi_ = tempChi;
      return true;
    } else {
      currentLambda_ *= ni_;
      ni_ *= 2;
      return false;
    }
  }

  double ni_ = 0;
  double currentLambda_ = 0;
  double currentChi_ = 0;
  /// Levenberg
  /// 计算LM算法的初始Lambda
  void ComputeLambdaInitLM() {
    ni_ = 2.;
    currentLambda_ = -1.;
    currentChi_ = 0.0;

    // // 计算出当前的总残差
    currentChi_ = Chi2();
    // for (const auto& edge : edges_) {
    //   currentChi_ += edge.second->Chi2();
    // }

    // 1. 第一步计算停止迭代条件stopThresholdLM_
    stopThresholdLM_ = 1e-6 * currentChi_; // 迭代条件为 误差下降 1e-6 倍

    // 取出H矩阵对角线的最大值
    double maxDiagonal = 0.;
    maxDiagonal = MyMatrixAMax<double, double, 3>(hessian_);
    double tau = 1e-5;
    // 2. 根据对角线最大值计算出currentLambda_
    currentLambda_ = tau * maxDiagonal; // 给到u0的初值
  }
  double hessian_[3][3] = {};
  double b_[3][1];
  double chi2_;
  double delta_x_[3][1];
  // 最大支持200条边
  static const int max_edge_curves_size = 200;
  int edge_curves_size = 0;
  CurveFittingEdge edge_curves[max_edge_curves_size];
  double stopThresholdLM_ = 0;
  double current_chi2_;
};
