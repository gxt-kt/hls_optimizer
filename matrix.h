#pragma once
// #include "ap_fixed.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>

#include "help.h"

template <typename T = double, size_t X = 1, size_t Y = 1> class MyMatrix {
public:
  using type = T;
  // type **data() { return data_; }
  type &operator()(size_t i, size_t j) { return data_[i][j]; }
  const type &operator()(size_t i, size_t j) const { return data_[i][j]; }
  type data_[X][Y] = {};

public:
  operator type **() { return data_; }
  MyMatrix<type, Y, X> Transpose() {
    MyMatrix<type, Y, X> res;
    MyMatrixTranspose<type, type, X, Y>(this->data_, res.data_);
    return res;
  }
  template <typename U, size_t M, size_t N>
  MyMatrix<T, X, N> operator*(const MyMatrix<U, M, N> &matrix) {
    MyMatrix<T, X, N> res;
    static_assert(M == Y);
    MyMatrixMultiple<T, U, T, X, Y, N>(this->data_, matrix.data_, res.data_);
    return res;
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const MyMatrix<T, X, Y> &matrix) {
    MATRIXDEBUG(matrix.data_);
    return os;
  }
};

template <typename T, size_t N> using Vector = MyMatrix<T, N, 1>;

// 专门设置顶点，预先定义好
class VertexCurveABC {
public:
  double parameters[3][1] = {{0}, {0}, {0}}; // abc
};

static VertexCurveABC verticies_[10];

/**
 * 以逆深度形式存储的顶点
 */
class VertexInverseDepth {
public:
  VertexInverseDepth() {}
  double parameters[7][1] = {}; // xyz xyzw

  void Plus(double delta[7][1]) {
    parameters[0][0] += delta[0][0];
    parameters[1][0] += delta[1][0];
    parameters[2][0] += delta[2][0];

    double q[4] = {};
    double dq[4] = {};
    dq[0] = delta[3][0] / 2;
    dq[1] = delta[4][0] / 2;
    dq[2] = delta[5][0] / 2;
    dq[3] = 1;
  }
};

template <int residual_dimension, int num_verticies> class Edge {
public:
  explicit Edge(){};

  ~Edge(){};

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

class EdgeReprojection : public Edge<2, 3> {
  static const int residual_dimension = 2;
  static const int num_verticies = 3;

  // 需要知道当前边对应的顶点
  // 第一个顶点是逆深度
  size_t inver_depth_idx = 0;
  // 第二三个顶点分别是i和j的位姿
  size_t pose_i_idx = 0;
  size_t pose_j_idx = 0;

  double pts_i_[3][1] = {};
  double pts_j_[3][1] = {};

  void ComputeResidual() {}
  void ComputeJacobians() {}

  // 总共三个雅可比
  double jacobians_0[1][3];
  double jacobians_1[1][3];
  double jacobians_2[1][3];
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
  // VertexCurveABC *verticies_[num_verticies];
  const int num_verticies_ = num_verticies;
};

class Problem {
public:
  Problem() {
    verticies_[0].parameters[0][0] = 0;
    verticies_[0].parameters[1][0] = 0;
    verticies_[0].parameters[2][0] = 0;
  }

public:
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
    // 需要先清空原有的hessian矩阵
    MyMatrixSet<double, double, 3, 3>(hessian_, 0);
    MyMatrixSet<double, double, 3, 1>(b_, 0);
    // 清空delta_x
    MyMatrixSet<double, double, 3, 1>(delta_x_, 0);

    for (int i = 0; i < max_edge_curves_size; i++) {
      if (i >= edge_curves_size)
        continue;
      edge_curves[i].ComputeResidual();
      edge_curves[i].ComputeJacobians();

      // 遍历所有顶点
      // for (int j = 0; j < edge_curves[i].num_verticies_; j++) {
      // }

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
    MyMatrixAdd<double, double, 3, 1, 3, 1, 0, 0>(verticies_[0].parameters,
                                                  delta_x_);
    // MATRIXDEBUG(edge_curves[0].verticies_[0].parameters);
    // std::cout << __PRETTY_FUNCTION__ << "end" << std::endl;
  }
  void RollbackStates() {
    // std::cout << __PRETTY_FUNCTION__ << "begin" << std::endl;
    // MATRIXDEBUG(edge_curves[0].verticies_[0].parameters);
    // MATRIXDEBUG(delta_x_);
    MyMatrixSub<double, double, 3, 1, 3, 1, 0, 0>(verticies_[0].parameters,
                                                  delta_x_);
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
      bool one_step_success = false;
      int false_cnt = 0;
      while (!one_step_success) { // 不断尝试 Lambda, 直到成功迭代一步
        // 更新Hx=b为(H+uI)x=b也就是H变为H+uI
        AddLambdatoHessianLM();
        // 解线性方程 Hx=b
        SolveLinearSystem();
        // 把H+uI恢复到原来的H
        RemoveLambdatoHessianLM();
        // 优化退出条件1： delta_x_ 很小则退出
        if (MySquaredNorm<double, 3>(delta_x_) <= 1e-8) {
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
        one_step_success = IsGoodStepInLM();
        std::cout << VAR(one_step_success) << std::endl;
        if (one_step_success) {
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
        std::cout << "std::sqrt(currentChi_) <= stopThresholdLM_" << std::endl;
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

void MatrixMultiple(unsigned int A[100], unsigned int B[100],
                    unsigned int C[100], unsigned int A_a[100],
                    unsigned int B_b[100]);
