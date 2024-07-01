#pragma once
#include "ap_fixed.h"

#include <iostream>

//constexpr int A_D1 =  4;
//constexpr int A_D2 =  4;
//constexpr int B_D1 =  A_D2;
//constexpr int B_D2 =  4;
//constexpr int C_D1 =  A_D1;
//constexpr int C_D2 =  B_D2;

////using A_T = ap_int<8>;
//using A_T = ap_int<8>;
//using B_T = ap_int<8>;
//using C_T = ap_int<8>;


//void MatrixMultiple(A_T A[A_D1][A_D2], B_T B[B_D1][B_D2], C_T C[C_D1][C_D2]);
void MatrixMultiple(ap_int<8> A[4][4], ap_int<8> B[4][4], ap_int<16> C[4][4]);

template <typename T_a,typename T_b,typename T_c,int x,int y,int z>
void MyMatrixMultiple(T_a A[x][y], T_a B[y][z], T_a C[x][z]) {
	#pragma HLS ARRAY_RESHAPE variable=B complete dim=1
	#pragma HLS ARRAY_RESHAPE variable=A complete dim=2
	for(int i=0;i<x;i++)
		for(int j=0;j<z;j++){
			#pragma HLS PIPELINE II=1
			C[i][j]=0;
			// 循环乘四次，并进行相加
			for(int k=0;k<y;k++){
				C[i][j]=C[i][j]+A[i][k]*B[k][j];
			}
		}
}

//template <typename T_a,typename T_b,typename T_c,int k,int n>
//void MyGevm(T_a A[k], T_a B[k][n], T_a C[1][n]) {
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

template <int residual_dimension,int num_verticies>
class Edge {
 public:
  explicit Edge(){};

  virtual ~Edge(){};

  /// 计算残差，由子类实现
  virtual void ComputeResidual() = 0;

  /// 计算雅可比，由子类实现
  /// 本后端不支持自动求导，需要实现每个子类的雅可比计算方法
  virtual void ComputeJacobians() = 0;

  //    ///计算该edge对Hession矩阵的影响，由子类实现
  //    virtual void ComputeHessionFactor() = 0;

  /// 计算平方误差，会乘以信息矩阵
  double Chi2() {
//	  return 1;
	  double result[1][residual_dimension];
	  MyMatrixMultiple<double,double,double,1,residual_dimension,residual_dimension>(residual_,information_,result);
	  double res[1][1];
	  MyMatrixMultiple<double,double,double,1,residual_dimension,1>(result,residual_,res);
	  return res[0][0];
  }

  void SetInformation(double information[residual_dimension][residual_dimension]) {
#pragma HLS PIPELINE
#pragma HLS UNROLL
	  for(int i=0;i<residual_dimension;i++) {
		  for(int j=0;j<residual_dimension;j++) {
			  information_[i][j]=information[i][j];
		  }
	  }
  }

 public:
//  unsigned long id_;                          // edge id
  int ordering_id_;                           // edge id in problem
//  std::vector<std::shared_ptr<Vertex>> verticies_;  // 该边对应的顶点
//  double residual_;
  double residual_[1][residual_dimension];                 // 残差
//  std::vector<MatXX>
//      jacobians_;  // 雅可比，每个雅可比维度是 residual x vertex[i]
//  MatXX information_;      // 信息矩阵
  double information_[residual_dimension][residual_dimension];
};


class CurveFittingEdge : public Edge<1,1>{
public:
	  double x_, y_;
	  void ComputeResidual() override {
		  residual_[0][0]=10;
	  }
	  void ComputeJacobians() override {
	  }
	  const int residual_dimension=1;
	  const int num_verticies=1;
};
