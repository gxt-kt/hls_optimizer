#pragma once

#include "matrix.h"
#include "vertex.h"

extern VertexCurveABC verticies_[];

extern VertexInverseDepth v_points[];

extern VertexPose v_poses[];

class EdgeReprojection;
extern EdgeReprojection e_reproject[];

template <int residual_dimension, int num_verticies> class Edge {
public:
  explicit Edge(){};

  ~Edge(){};

  /// 计算平方误差，会乘以信息矩阵
  float Chi2() {
    Matrix<float, 1, 1> res;
    res = residual_.transpose() * information_ * residual_;
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
public:
  static const int residual_dimension = 2;
  // 三个顶点分别是：
  // 1. 特征点的逆深度
  // 2. 特阵点在第0帧相机的投影
  // 3. 特阵点在第k帧相机的投影
  static const int num_verticies = 3;

  // 需要知道当前边对应的顶点
  // 第一个顶点是逆深度
  size_t inver_depth_idx = 0;
  // 第0个顶点是逆深度
  size_t v_idx0 = 0;
  // 第二三个顶点分别是i和j的位姿
  size_t v_idx1 = 0;
  size_t v_idx2 = 0;

  EdgeReprojection() {}

  void ComputeResidual() {
    float inv_dep_i = v_points[v_idx0].parameters[0];

    auto param_i = v_poses[v_idx1].parameters;
    Quaternion<float> Qi;
    Qi[0] = param_i[6];
    Qi[1] = param_i[3];
    Qi[2] = param_i[4];
    Qi[3] = param_i[5];
    Matrix<float, 3, 1> Pi = param_i.head<3>();

    auto param_j = v_poses[v_idx2].parameters;
    Quaternion<float> Qj;
    Qj[0] = param_j[6];
    Qj[1] = param_j[3];
    Qj[2] = param_j[4];
    Qj[3] = param_j[5];
    Matrix<float, 3, 1> Pj = param_j.head<3>();

    auto pts_camera_i = pts_i_ / inv_dep_i;
    auto pts_imu_i = qic * pts_camera_i + tic;
    auto pts_w = Qi * pts_imu_i + Pi;
    auto pts_imu_j = Qj.inverse() * (pts_w - Pj);
    // std::cout << VAR(Qj) << std::endl;
    // std::cout << VAR(Qj.inverse()) << std::endl;
    // std::cout << VAR(pts_w-Pj) << std::endl;
    // std::cout << VAR(pts_imu_j) << std::endl;
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
  void ComputeJacobians() {
    float inv_dep_i = v_points[v_idx0].parameters[0];

    auto param_i = v_poses[v_idx1].parameters;
    Quaternion<float> Qi;
    Qi[0] = param_i[6];
    Qi[1] = param_i[3];
    Qi[2] = param_i[4];
    Qi[3] = param_i[5];
    Matrix<float, 3, 1> Pi = param_i.head<3>();

    auto param_j = v_poses[v_idx2].parameters;
    Quaternion<float> Qj;
    Qj[0] = param_j[6];
    Qj[1] = param_j[3];
    Qj[2] = param_j[4];
    Qj[3] = param_j[5];
    Matrix<float, 3, 1> Pj = param_j.head<3>();

    auto pts_camera_i = pts_i_ / inv_dep_i;
    auto pts_imu_i = qic * pts_camera_i + tic;
    // gxt:
    auto pts_imu_i_double = pts_imu_i;
    auto pts_w = Qi * pts_imu_i + Pi;
    auto pts_imu_j = Qj.inverse() * (pts_w - Pj);
    // gxt:
    auto pts_imu_j_double = pts_imu_j;
    auto pts_camera_j = qic.inverse() * (pts_imu_j - tic);
    //
    auto dep_j = pts_camera_j.z();
    //
    auto Ri = Qi.toRotationMatrix();
    auto Rj = Qj.toRotationMatrix();
    auto ric = qic.toRotationMatrix();
    Matrix<float, 2, 3> reduce;
    // // gDebugWarn() << G_FILE_LINE;
    reduce(0, 0) = 1. / dep_j;
    reduce(0, 1) = 0;
    reduce(0, 2) = -pts_camera_j(0) / (dep_j * dep_j);
    reduce(1, 0) = 0;
    reduce(1, 1) = 1. / dep_j;
    reduce(1, 2) = -pts_camera_j(1) / (dep_j * dep_j);
    // // gDebugWarn() << G_FILE_LINE;
    // reduce = information_ * reduce;
    //
    Matrix<float, 2, 6> jacobian_pose_i;
    Matrix<float, 3, 6> jaco_i;
    jaco_i.SetleftCols<3>(ric.transpose() * Rj.transpose());
    jaco_i.SetrightCols<3>(ric.transpose() * Rj.transpose() * Ri * (-1) *
                           pts_imu_i_double.hat());
    // std::cout << VAR(pts_imu_i_double) << std::endl;
    // std::cout << VAR(pts_imu_i_double.hat()) << std::endl;
    // std::cout << VAR(ric.transpose() * Rj.transpose() * Ri * (-1) * pts_imu_i_double.hat()) << std::endl;
    // std::cout << VAR(jaco_i) << std::endl;
    jacobian_pose_i.SetleftCols<6>(reduce * jaco_i);
    //
    Matrix<float, 2, 6> jacobian_pose_j;
    Matrix<float, 3, 6> jaco_j;
    jaco_j.SetleftCols<3>(ric.transpose() * (-1) * Rj.transpose());
    jaco_j.SetrightCols<3>(ric.transpose() * pts_imu_j_double.hat());
    jacobian_pose_j.SetleftCols<6>(reduce * jaco_j);

    Matrix<float, 2, 1> jacobian_feature;
    jacobian_feature = reduce * ric.transpose() * Rj.transpose() * Ri * ric *
                       pts_i_ * -1.0 / (inv_dep_i * inv_dep_i);
    //
    jacobians_0 = jacobian_feature;
    jacobians_1 = jacobian_pose_i;
    jacobians_2 = jacobian_pose_j;
  }

  Matrix<float, 3, 1> pts_i_;
  Matrix<float, 3, 1> pts_j_;

  Quaternion<float> qic;
  Matrix<float, 3, 1> tic;

  // 总共三个雅可比
  Matrix<float, 2, 1> jacobians_0;
  Matrix<float, 2, 6> jacobians_1;
  Matrix<float, 2, 6> jacobians_2;
};
