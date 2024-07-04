#include "matrix.h"
#include <exception>
#include <iostream>
// #include "ap_fixed.h"

// #include "hls_stream.h"
// #include "xf_blas.hpp"

// using namespace xf::blas;

// void TestStream() {
// 	using t_DataType = ap_int<8>;
// 	t_DataType in[20]={10,53,2,34,-6,-123,-10,3};
// 	t_DataType out[20]={10,53,2,34,-6,-123,-10,3};
// 	const int t_LogParEntries =2;
// 	hls::stream<typename WideType<t_DataType, (1 <<
// t_LogParEntries)>::t_TypeInt> l_str; 	readVec2Stream<t_DataType, 1 <<
// t_LogParEntries>(in, 16, l_str); 	writeStream2Vec<t_DataType, 1 <<
// t_LogParEntries>(l_str,16,out); 	for(int i=0;i<20;i++) {
// std::cout << "i=" << out[i] << std::endl;
// 	}
//
// }

// void Test() {
//     unsigned int p_n = 8;
// 	using t_DataType = ap_int<8>;
//
// 	t_DataType p_x[20]={10,53,2,34,-6,-123,-10,3};
// 	// 定义数据类型和并行处理元素个数
// 	constexpr unsigned int t_LogParEntries = 2;
// 	std::cout << "1 << t_LogParEntries=" << (1 << t_LogParEntries) <<
// std::endl; 	using t_IndexType = unsigned int;
//
// 	// 创建输入向量并写入 HLS 流
// 	hls::stream<typename WideType<t_DataType, (1 <<
// t_LogParEntries)>::t_TypeInt> l_str; 	readVec2Stream<t_DataType, 1 <<
// t_LogParEntries>(p_x, p_n, l_str);
//
// 	// 调用 amin 函数并检查结果
// 	t_IndexType l_res;
// 	amin<t_DataType, t_LogParEntries, t_IndexType>(p_n, l_str, l_res);
// 	std::cout << "Minimum element index: " << l_res << std::endl;
// }

// template <typename BLAS_dataType,int BLAS_matrixSizeA=100,int
// BLAS_matrixSizeB=100,int BLAS_matrixSizeC=100> void fun(uint32_t p_m,
//             uint32_t p_n,
//             uint32_t p_k,
//             BLAS_dataType p_alpha,
//             BLAS_dataType p_beta,
//             BLAS_dataType p_A[BLAS_matrixSizeA],
//             BLAS_dataType p_B[BLAS_matrixSizeB],
//             BLAS_dataType p_C[BLAS_matrixSizeC],
//             BLAS_dataType p_R[BLAS_matrixSizeC]) {
// 	   const int BLAS_parEntries =1;
//
//    hls::stream<typename WideType<BLAS_dataType, BLAS_parEntries>::t_TypeInt>
//    l_strA; hls::stream<typename WideType<BLAS_dataType,
//    BLAS_parEntries>::t_TypeInt> l_strB; hls::stream<typename
//    WideType<BLAS_dataType, BLAS_parEntries>::t_TypeInt> l_strC;
//    hls::stream<typename WideType<BLAS_dataType, BLAS_parEntries>::t_TypeInt>
//    l_strSum;
// //#pragma HLS DATAFLOW
//    gemmMatAMover<BLAS_dataType, BLAS_parEntries>(p_A, p_m, p_n, p_k, l_strA);
//    gemmMatBMover<BLAS_dataType, BLAS_parEntries>(p_B, p_m, p_n, p_k, l_strB);
//    readVec2Stream<BLAS_dataType, BLAS_parEntries>(p_C, p_m * p_n, l_strC);
//    const int BLAS_k =1;
//    gemm<BLAS_dataType, BLAS_k, BLAS_parEntries, BLAS_matrixSizeC>(p_m, p_n,
//    p_k, p_alpha, l_strA, l_strB, p_beta,
//                                                                   l_strC,
//                                                                   l_strSum);
//    writeStream2Vec<BLAS_dataType, BLAS_parEntries>(l_strSum, p_m * p_n, p_R);
// };
//
// void TestGemm() {
// 	using BLAS_dataType = ap_int<32>;
// //	using BLAS_dataType = float;
// 	uint32_t p_m = 3;
// 	uint32_t p_n = 3;
// 	const uint32_t p_k = 1;
// 	const unsigned int BLAS_parEntries = 1;
//
// 	BLAS_dataType p_A[9]={1,2,3,4,5,6,7,8,9}; //3*2
// 	BLAS_dataType p_B[9]={1,2,3,4,5,6,7,8,9};//2*4
// 	BLAS_dataType p_C[9]={1,2,3,4,5,6,7,8,9};//3*4
// 	BLAS_dataType p_R[9]={1,2,3,4,5,6,7,8,9};//3*4
//
// #if 1
//     hls::stream<typename WideType<BLAS_dataType, BLAS_parEntries>::t_TypeInt>
//     l_strA; hls::stream<typename WideType<BLAS_dataType,
//     BLAS_parEntries>::t_TypeInt> l_strB; hls::stream<typename
//     WideType<BLAS_dataType, BLAS_parEntries>::t_TypeInt> l_strC;
//     hls::stream<typename WideType<BLAS_dataType, BLAS_parEntries>::t_TypeInt>
//     l_strSum; gemmMatAMover<BLAS_dataType, 1>(p_A, p_m, p_n, p_k, l_strA);
//     gemmMatBMover<BLAS_dataType, 1>(p_B, p_m, p_n, p_k, l_strB);
// //    readVec2Stream<BLAS_dataType, BLAS_parEntries>(p_A, p_m * p_n, l_strA);
// //    readVec2Stream<BLAS_dataType, BLAS_parEntries>(p_B, p_n * p_k, l_strB);
// //    gemm<BLAS_dataType, p_k, BLAS_parEntries>(p_m, p_n,
// p_k,l_strA,l_strB,l_strC);
// //    gemm<BLAS_dataType, 1, BLAS_parEntries, 3>(p_m, p_n, p_k, 1, l_strA,
// l_strB, 1,
// //                                                                   l_strC,
// l_strSum);
//
// //    gemm<BLAS_dataType, 2, BLAS_parEntries, 100>(p_m, p_n,p_k, 1, l_strA,
// l_strB, 1,
// //                                                                   l_strC,
// l_strSum);
//     writeStream2Vec<BLAS_dataType, BLAS_parEntries>(l_strA, p_m * p_n, p_A);
//     writeStream2Vec<BLAS_dataType, BLAS_parEntries>(l_strB, p_m * p_n, p_B);
//     writeStream2Vec<BLAS_dataType, BLAS_parEntries>(l_strC, p_m * p_n, p_C);
//     std::cout <<"=======p_A========" << std::endl;
//     for(int i=0;i<9;i++) {
//     	std::cout << "i=" << p_A[i] << std::endl;
//     }
//     std::cout <<"=======p_B========" << std::endl;
//         for(int i=0;i<9;i++) {
//         	std::cout << "i=" << p_B[i] << std::endl;
//         }
// 	std::cout <<"=======p_C========" << std::endl;
// 		   for(int i=0;i<9;i++) {
// 			std::cout << "i=" << p_C[i] << std::endl;
// 		   }
// #endif
// //    fun<BLAS_dataType>(p_m,p_n,p_k,1,1,p_A,p_B,p_C,p_R);
// //    for(int i=0;i<20;i++) {
// //        	std::cout << "i=" << p_C[i] << std::endl;
// //        }
// //    for(int i=0;i<20;i++) {
// //            	std::cout << "i=" << p_R[i] << std::endl;
// //            }
// }

// constexpr int A_D1 =  4;
// constexpr int A_D2 =  4;
// constexpr int B_D1 =  A_D2;
// constexpr int B_D2 =  4;
// constexpr int C_D1 =  A_D1;
// constexpr int C_D2 =  B_D2;

// using A_T = ap_int<8>;
// using A_T = ap_int<8>;
// using B_T = ap_int<8>;
// using C_T = ap_int<16>;

void TestMyMatrixMultiple() {
  using type_a = float;
  using type_b = type_a;
  using type_c = type_a;

  type_a p_A[3][2] = {{1, 2}, {3, 4}, {5, 6}}; // 3*2
  type_b p_B[2][3] = {{1, 2, 3}, {4, 5, 6}};   // 2*4
  type_c p_C[3][3];                            // 3*4
  MyMatrixMultiple<type_a, type_b, type_c, 3, 2, 3>(p_A, p_B, p_C);
  for (int i = 0; i < 3; i++) {
    std::cout << "0: ";
    for (int j = 0; j < 3; j++) {
      std::cout << p_C[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

// void TestMyGevm() {
//	using type_a=ap_int<8>;
//	using type_b=type_a;
//	using type_c=type_a;
//
//	type_a p_A[1]={1,2}; //3*2
//	type_b p_B[2][3]={{1,2,3},{4,5,6}}; //2*4
//	type_c p_C[3]; //3*4
//	MyMatrixMultiple<type_a,type_b,type_c,3,2,3>(p_A,p_B,p_C);
//	for(int i=0;i<3;i++) {
//		std::cout << "0: ";
//		for(int j=0;j<3;j++) {
//			std::cout << p_C[i][j] << " ";
//		}
//		std::cout << std::endl;
//	}
//}
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int main() {
  // ap_int<8> a[4][4]={{0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15}};
  // ap_int<8> b[4][4]={{0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15}};
  // ap_int<16> c[4][4];
  // MatrixMultiple(a,b,c);
  // for(int i=0;i<4;i++) {
  // 	for(int j=0;j<4;j++) {
  // 		std::cout << "C["<<i<<","<<j<<"]="<<c[i][j] << "  ";
  // 	}
  // 	std::cout << std::endl;
  // }
  double a = 1, b = 2, c = 1;
  VertexCurveABC abc;
  abc.parameters[0][0] = 0;
  abc.parameters[1][0] = 0;
  abc.parameters[2][0] = 0;

  // 构建问题
  Problem problem;

  // 构建100条边
  for (int i = 0; i < 100; i++) {
    problem.edge_curves[i].verticies_[0].parameters[0][0] = 0;
    problem.edge_curves[i].verticies_[0].parameters[1][0] = 0;
    problem.edge_curves[i].verticies_[0].parameters[2][0] = 0;
    double info[1][1];
    info[0][0] = 1;
    problem.edge_curves[i].SetInformation(info);
    double x = i / 100.0;
    problem.edge_curves[i].x_ = x;
    problem.edge_curves[i].y_ = std::exp(a * x * x + b * x + c);
    problem.edge_curves_size++;
  }

  problem.CalculateResidual();
  // std::terminate();
  std::cout << VAR(problem.Chi2()) << std::endl;
  // problem.MakeHessian();
  problem.Solve();
  // CurveFittingEdge edge;
  // for (int i = 0; i < 100; i++) {
  //   edge[i].ComputeResidual();
  //   edge[i].ComputeJacobians();
  //   std::cout << VAR(edge[i].Chi2()) << std::endl;
  // }
  // std::cout << "edge.Chi2()=" << edge.Chi2() << std::endl;;

  // Test();
  // TestStream();
  // TestMyMatrixMultiple();
  //	TestGemm();
  //	edge.Fun();
  //	edge.Print();
}
