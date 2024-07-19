#pragma once

#include "edge.h"

class Problem {
public:
  Problem() { verticies_[0].parameters.SetValue(0); }

public:
  void CalculateResidual() {
    for (int i = 0; i < max_edge_curves_size; i++) {
      if (i >= edge_curves_size)
        continue;
      edge_curves[i].ComputeResidual();
    }
  }
  float Chi2() {
    CalculateResidual();
    float res = 0;
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
    // hessian_=hessian_*0;
    // 需要先清空原有的hessian矩阵
    hessian_.SetValue(0);
    b_.SetValue(0);
    // 清空delta_x
    delta_x_.SetValue(0);

    for (int i = 0; i < max_edge_curves_size; i++) {
      if (i >= edge_curves_size)
        continue;
      edge_curves[i].ComputeResidual();
      edge_curves[i].ComputeJacobians();

      Matrix<float, 3, 1> JtW;

      JtW =
          edge_curves[i].jacobians_0.transpose() * edge_curves[i].information_;

      Matrix<float, 3, 3> hessian = JtW * edge_curves[i].jacobians_0;

      hessian_ = hessian_ + hessian;

      Matrix<float, 3, 1> bb = JtW * edge_curves[i].residual_;

      b_ = b_ - bb;
    }
    // MATRIXDEBUG(hessian_);
    // MATRIXDEBUG(b_);
  }
  void SolveLinearSystem() {
    // TODO:
    // 在这里我们假设已经求解完成了
    // 输入是hessian_ 和 b_
    // 输出是 delta_x_
    // delta_x_ = Hessian_.inverse() * b_;
    //
    LdltSolve<float, float, float, 3>(hessian_.data_, delta_x_.data_, b_.data_,
                                      false);
    // float h_dx[3][1] = {};
    Matrix<float, 3, 1> h_dx = hessian_ * delta_x_;
    // MyMatrixMultiple<float, float, float, 3, 3, 1>(hessian_.data_,
    //                                                   delta_x_.data_, h_dx);
    // MATRIXDEBUG(delta_x_);
    // MATRIXDEBUG(h_dx);
  }
  void UpdateStates() {
    verticies_[0].parameters = verticies_[0].parameters + delta_x_;
  }
  void RollbackStates() {
    verticies_[0].parameters = verticies_[0].parameters - delta_x_;
  }

  /// Hessian 对角线加上或者减去  Lambda
  void AddLambdatoHessianLM() {
    MyMatrixAddNumber<float, float, 3, 3, 0, 3>(hessian_.data_, currentLambda_);
  }
  void RemoveLambdatoHessianLM() {
    MyMatrixAddNumber<float, float, 3, 3, 0, 3>(hessian_.data_,
                                                -currentLambda_);
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
        if (MySquaredNorm<float, 3>(delta_x_.data_) <= 1e-8) {
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
    float scale = 0;

    Matrix<float, 1, 1> dot;
    Matrix<float, 1, 3> delta_x_transpose;

    MyMatrixTranspose<float, float, 3, 1>(delta_x_.data_,
                                          delta_x_transpose.data_);

    Matrix<float, 3, 1> delta_x_lambda_ = delta_x_ * currentLambda_;
    delta_x_lambda_ + b_;
    dot = delta_x_.transpose() * delta_x_lambda_;

    scale = dot(0, 0);

    scale += 1e-3; // make sure it's non-zero :)

    std::cout << VAR(scale) << std::endl;

    // 统计所有的残差
    float tempChi = Chi2();

    std::cout << VAR(currentChi_, tempChi) << std::endl;

    float rho = (currentChi_ - tempChi) / scale;
    std::cout << VAR(rho) << std::endl;

    if (rho > 0 && std::isfinite(tempChi)) { // last step was good, 误差在下降
      float alpha = 1. - std::pow((2 * rho - 1), 3);
      alpha = std::min(alpha, float(2.0 / 3.0));
      float scaleFactor = std::max(float(1.0 / 3.0), alpha);
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

  float ni_ = 0;
  float currentLambda_ = 0;
  float currentChi_ = 0;

  /// 计算LM算法的初始Lambda
  void ComputeLambdaInitLM() {
    ni_ = 2.;
    currentLambda_ = -1.;
    currentChi_ = 0.0;

    // // 计算出当前的总残差
    currentChi_ = Chi2();

    // 1. 第一步计算停止迭代条件stopThresholdLM_
    stopThresholdLM_ = 1e-6 * currentChi_; // 迭代条件为 误差下降 1e-6 倍

    // 取出H矩阵对角线的最大值
    float maxDiagonal = 0.;
    maxDiagonal = MyMatrixAMax<float, float, 3>(hessian_.data_);
    float tau = 1e-5;
    // 2. 根据对角线最大值计算出currentLambda_
    currentLambda_ = tau * maxDiagonal; // 给到u0的初值
  }
  // float hessian_[3][3] = {};
  Matrix<float, 3, 3> hessian_;
  // float b_[3][1];
  Matrix<float, 3, 1> b_;
  float chi2_;
  // float delta_x_[3][1];
  Matrix<float, 3, 1> delta_x_;
  // 最大支持200条边
  static const int max_edge_curves_size = 200;
  int edge_curves_size = 0;
  CurveFittingEdge edge_curves[max_edge_curves_size];
  float stopThresholdLM_ = 0;
  float current_chi2_;
};


const int dimen=3*6+1*20;
class BAProblem {
public:
  // BAProblem() { }

public:
  void CalculateResidual() {
    for (int i = 0; i < max_edge_reproject_size; i++) {
      if (i >= edge_reproject_size)
        continue;
      // edge_curves[i].ComputeResidual();
      e_reproject[i].ComputeResidual();
    }
  }
  float Chi2() {
    CalculateResidual();
    float res = 0;
    for (int i = 0; i < max_edge_reproject_size; i++) {
      // std::cout << VAR(res) << std::endl;
      if (i >= edge_reproject_size)
        continue;
      res += e_reproject[i].Chi2();
      // res += edge_curves[i].Chi2();
      // std::cout << VAR(res) << std::endl;
    }
    return res;
  }
  void MakeHessian() {
    // hessian_=hessian_*0;
    // 需要先清空原有的hessian矩阵
    hessian_.SetValue(0);
    b_.SetValue(0);
    // 清空delta_x
    delta_x_.SetValue(0);

    for (int i = 0; i < max_edge_reproject_size; i++) {
      if (i >= edge_reproject_size)
        continue;
      e_reproject[i].ComputeResidual();
      // std::cout << "residual i=" << i << " "
      //           << e_reproject[i].residual_.transpose() << std::endl;
      e_reproject[i].ComputeJacobians();
      // std::cout << "jacobians0 i=" << i << " " << e_reproject[i].jacobians_0
      //           << std::endl;
      // std::cout << "jacobians1 i=" << i << " " << e_reproject[i].jacobians_1
      //           << std::endl;
      // std::cout << "jacobians2 i=" << i << " " << e_reproject[i].jacobians_2
      //           << std::endl;

      // Matrix<float, 3, 1> JtW;
      //
      // JtW =
      //     e_reproject[i].jacobians_0.transpose() *
      //     e_reproject[i].information_;
      //
      // Matrix<float, 3, 3> hessian = JtW * e_reproject[i].jacobians_0;
      //
      // hessian_ = hessian_ + hessian;
      //
      // Matrix<float, 3, 1> bb = JtW * edge_curves[i].residual_;
      //
      // b_ = b_ - bb;
    }
    // MATRIXDEBUG(hessian_);
    // MATRIXDEBUG(b_);
  }
  void SolveLinearSystem() {
    // TODO:
    // 在这里我们假设已经求解完成了
    // 输入是hessian_ 和 b_
    // 输出是 delta_x_
    // delta_x_ = Hessian_.inverse() * b_;
    //
    // LdltSolve<float, float, float, 3>(hessian_.data_, delta_x_.data_, b_.data_,
    //                                   false);
    // float h_dx[3][1] = {};
    auto h_dx = hessian_ * delta_x_;
    // MyMatrixMultiple<float, float, float, 3, 3, 1>(hessian_.data_,
    //                                                   delta_x_.data_, h_dx);
    // MATRIXDEBUG(delta_x_);
    // MATRIXDEBUG(h_dx);
  }
  void UpdateStates() {
    // verticies_[0].parameters = verticies_[0].parameters + delta_x_;
  }
  void RollbackStates() {
    // verticies_[0].parameters = verticies_[0].parameters - delta_x_;
  }

  /// Hessian 对角线加上或者减去  Lambda
  void AddLambdatoHessianLM() {
    // MyMatrixAddNumber<float, float, 3, 3, 0, 3>(hessian_.data_, currentLambda_);
    hessian_.AddDiagonal(currentLambda_);
  }
  void RemoveLambdatoHessianLM() {
    // MyMatrixAddNumber<float, float, 3, 3, 0, 3>(hessian_.data_,
    //                                             -currentLambda_);
    hessian_.AddDiagonal(-currentLambda_);
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
        if (MySquaredNorm<float, 3>(delta_x_.data_) <= 1e-8) {
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
    float scale = 0;

    Matrix<float, 1, 1> dot;

    auto delta_x_lambda_ = delta_x_ * currentLambda_;
    delta_x_lambda_ + b_;
    dot = delta_x_.transpose() * delta_x_lambda_;

    scale = dot(0, 0);

    scale += 1e-3; // make sure it's non-zero :)

    std::cout << VAR(scale) << std::endl;

    // 统计所有的残差
    float tempChi = Chi2();

    std::cout << VAR(currentChi_, tempChi) << std::endl;

    float rho = (currentChi_ - tempChi) / scale;
    std::cout << VAR(rho) << std::endl;

    if (rho > 0 && std::isfinite(tempChi)) { // last step was good, 误差在下降
      float alpha = 1. - std::pow((2 * rho - 1), 3);
      alpha = std::min(alpha, float(2.0 / 3.0));
      float scaleFactor = std::max(float(1.0 / 3.0), alpha);
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

  float ni_ = 0;
  float currentLambda_ = 0;
  float currentChi_ = 0;

  /// 计算LM算法的初始Lambda
  void ComputeLambdaInitLM() {
    ni_ = 2.;
    currentLambda_ = -1.;
    currentChi_ = 0.0;

    // // 计算出当前的总残差
    currentChi_ = Chi2();

    // 1. 第一步计算停止迭代条件stopThresholdLM_
    stopThresholdLM_ = 1e-6 * currentChi_; // 迭代条件为 误差下降 1e-6 倍

    // 取出H矩阵对角线的最大值
    float maxDiagonal = 0.;
    maxDiagonal = hessian_.MaxDiagonal();
    float tau = 1e-5;
    // 2. 根据对角线最大值计算出currentLambda_
    currentLambda_ = tau * maxDiagonal; // 给到u0的初值
  }
  // 维度和顶点个数有关 3*6+1*20
  Matrix<float, dimen, dimen> hessian_;
  Matrix<float, dimen, 1> b_;
  float chi2_;
  Matrix<float, dimen, 1> delta_x_;
  // 最大支持200条边
  static const int max_edge_reproject_size = 100;
  int edge_reproject_size = 0;
  float stopThresholdLM_ = 0;
  float current_chi2_;
};
