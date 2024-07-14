#pragma once

#include "matrix.h"
#include "vertex.h"

extern VertexCurveABC verticies_[];

extern VertexInverseDepth v_points[];

extern VertexPose v_poses[];

template <int residual_dimension, int num_verticies> class Edge {
public:
  explicit Edge(){};

  ~Edge(){};

  /// 计算平方误差，会乘以信息矩阵
  float Chi2() {
    Matrix<float, 1, 1> res;
    res = residual_ * information_ * residual_;
    return res(0, 0);
  }

  void
  SetInformation(float information[residual_dimension][residual_dimension]) {
#pragma HLS PIPELINE
    for (int i = 0; i < residual_dimension; i++) {
      for (int j = 0; j < residual_dimension; j++) {
        information_.data_[i][j] = information[i][j];
      }
    }
  }

public:
  //  unsigned long id_;                          // edge id
  int ordering_id_; // edge id in problem
  //  std::vector<std::shared_ptr<Vertex>> verticies_;  // 该边对应的顶点
  //  float residual_;
  // float residual_[1][residual_dimension] = {0}; // 残差
  Matrix<float, residual_dimension, 1> residual_; // 残差
                                                  //  std::vector<MatXX>
  //      jacobians_;  // 雅可比，每个雅可比维度是 residual x vertex[i]
  //  MatXX information_;      // 信息矩阵
  // float information_[residual_dimension][residual_dimension] = {{0}};
  Matrix<float, residual_dimension, residual_dimension> information_; // 残差
};

// template <int residual_dimension, int num_verticies>
class CurveFittingEdge : public Edge<1, 1> {
public:
  static const int residual_dimension = 1;
  static const int num_verticies = 1;
  float x_ = 0, y_ = 0;
  void ComputeResidual() {
    // GetResidual();
    // residual_[0][0] = 10;
    // std::cout <<
    // VAR(verticies_[0].parameters[0][0],verticies_[0].parameters[1][0],verticies_[0].parameters[2][0])
    // << std::endl;
    residual_(0, 0) = std::exp(verticies_[0].parameters(0, 0) * x_ * x_ +
                               verticies_[0].parameters(1, 0) * x_ +
                               verticies_[0].parameters(2, 0)) -
                      y_;
    // std::cout << VAR(residual_[0][0]) << std::endl;
  }
  void ComputeJacobians() {
    float exp_x = std::exp(verticies_[0].parameters(0, 0) * x_ * x_ +
                           verticies_[0].parameters(1, 0) * x_ +
                           verticies_[0].parameters(2, 0));
    // std::cout << VAR(exp_y) << std::endl;
    jacobians_0(0, 0) = x_ * x_ * exp_x;
    jacobians_0(0, 1) = x_ * exp_x;
    jacobians_0(0, 2) = exp_x;
    // MATRIXDEBUG(jacobians_0);
  }

  // float jacobians_0[1][3];
  Matrix<float, 1, 3> jacobians_0;
  // const int residual_dimension=1;
  // const int num_verticies=1;
  // VertexCurveABC *verticies_[num_verticies];
  const int num_verticies_ = num_verticies;
};

class EdgeReprojection : public Edge<2, 3> {
  static const int residual_dimension = 2;
  static const int num_verticies = 3;

  // 需要知道当前边对应的顶点
  // 第一个顶点是逆深度
  size_t inver_depth_idx = 0;
  // 第0个顶点是逆深度
  size_t pose_idx0 = 0;
  // 第二三个顶点分别是i和j的位姿
  size_t pose_idx1 = 0;
  size_t pose_idx2 = 0;

  void ComputeResidual() {
    float inv_dep_i = v_points[pose_idx0].parameters[0];

    auto param_i = v_poses[pose_idx1].parameters;
    Quaternion<float> Qi;
    Qi[0] = param_i[6];
    Qi[1] = param_i[3];
    Qi[2] = param_i[4];
    Qi[3] = param_i[5];
    Matrix<float, 3, 1> Pi;
    Pi[0] = param_i[0];
    Pi[1] = param_i[1];
    Pi[2] = param_i[2];

    auto param_j = v_poses[pose_idx2].parameters;
    Quaternion<float> Qj;
    Qj[0] = param_i[6];
    Qj[1] = param_i[3];
    Qj[2] = param_i[4];
    Qj[3] = param_i[5];
    Matrix<float, 3, 1> Pj;
    Pi[0] = param_j[0];
    Pi[1] = param_j[1];
    Pi[2] = param_j[2];

    auto pts_camera_i = pts_i_ / inv_dep_i;
    auto pts_imu_i = qic * pts_camera_i + tic;
    auto pts_w = Qi * pts_imu_i + Pi;
    auto pts_imu_j = Qj.inverse() * (pts_w - Pj);
    auto pts_camera_j = qic.inverse() * (pts_imu_j - tic);

    float dep_j = pts_camera_j.z();
    residual_ = (pts_camera_j / dep_j).head<2>() -
                pts_j_.head<2>(); /// J^t * J * delta_x = - J^t * r
  }
  void SetTranslationImuFromCamera(const Quaternion<float> &qic_,
                                   const Matrix<float, 3, 1> &tic_) {
    qic = qic_;
    tic = tic_;
  }
  void ComputeJacobians() {}

  Matrix<float, 3, 1> pts_i_;
  Matrix<float, 3, 1> pts_j_;

  Quaternion<float> qic;
  Matrix<float, 3, 1> tic;

  // 总共三个雅可比
  Matrix<float, 1, 3> jacobians_0;
  Matrix<float, 1, 3> jacobians_1;
  Matrix<float, 1, 3> jacobians_2;
};
