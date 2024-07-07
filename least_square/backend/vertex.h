#pragma once

#include "eigen_types.hpp"

class Vertex {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /**
   * 构造函数
   * @param num_dimension 顶点自身维度
   * @param local_dimension 本地参数化维度，为-1时认为与本身维度一样
   */
  explicit Vertex(int num_dimension, int local_dimension = -1);

  virtual ~Vertex(){};

  /// 返回变量维度
  int Dimension() const { return parameters_.rows(); };

  /// 返回变量本地维度
  int LocalDimension() const { return local_dimension_; };

  /// 该顶点的id
  unsigned long Id() const { return id_; }

  /// 返回参数值
  VecX Parameters() const { return parameters_; }

  /// 返回参数值的引用
  VecX& Parameters() { return parameters_; }

  /// 设置参数值
  void SetParameters(const VecX& params) { parameters_ = params; }

  /// 加法，可重定义
  /// 默认是向量加
  virtual void Plus(const VecX& delta) { parameters_ += delta; };

  /// 返回顶点的名称，在子类中实现
  virtual std::string TypeInfo() const = 0;  // 纯虚函数来了

  /// 获取Id
  int OrderingId() const { return ordering_id_; }

  /// 设置在problem排序的位置Id,id代表在矩阵H中的第几列
  void SetOrderingId(unsigned int id) { ordering_id_ = id; }

  /// 固定该点的估计值
  void SetFixed(bool fixed = true) { fixed_ = fixed; }

  /// 测试该点是否被固定
  bool IsFixed() const { return fixed_; }

 protected:
  VecX parameters_;      // 实际存储的变量值
  int local_dimension_;  // 局部参数化维度
  unsigned long id_;     // 顶点的id，自动生成

  /// ordering id是在problem中排序后的id，用于寻找雅可比对应块
  /// ordering id带有维度信息，例如ordering_id=6则对应Hessian中的第6列
  /// 从零开始
  unsigned long ordering_id_ = 0;
  bool fixed_ = false;  // 是否固定
};
