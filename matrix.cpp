#include "matrix.h"

// 参考 https://blog.csdn.net/qq_45364953/article/details/127641923
//
//void MatrixMultiple(A_T A[A_D1][A_D2], B_T B[B_D1][B_D2], C_T C[C_D1][C_D2]) {
//#pragma HLS INTERFACE s_axilite port=B
//#pragma HLS INTERFACE s_axilite port=C
//#pragma HLS INTERFACE s_axilite port=A
//#pragma HLS INTERFACE s_axilite port=return
//    // partiton 优化语句
//	#pragma HLS ARRAY_PARTITION variable=B complete dim=1
//	#pragma HLS ARRAY_PARTITION variable=A complete dim=2
////#pragma HLS INTERFACE s_axilite port=return
////#pragma HLS INTERFACE s_axilite port=B
////#pragma HLS INTERFACE s_axilite port=C
////#pragma HLS INTERFACE s_axilite port=A
////#pragma HLS ARRAY_RESHAPE variable=B complete dim=1
////#pragma HLS ARRAY_RESHAPE variable=A complete dim=2
//	for(int i = 0; i < A_D1; ++ i){
//		for(int j = 0; j < B_D2; ++ j){
//			#pragma HLS PIPELINE II=1
//			C[i][j] = 0;
//			for(int k = 0; k < B_D1; ++ k){
//				C[i][j] += (C_T)A[i][k]*B[k][j];
//			}
//		}
//	}
//}

//void MatrixMultiple(ap_int<8> A[4][4], ap_int<8> B[4][4], ap_int<16> C[4][4]){
//	#pragma HLS INTERFACE s_axilite port=return
//	#pragma HLS INTERFACE s_axilite port=B
//	#pragma HLS INTERFACE s_axilite port=C
//	#pragma HLS INTERFACE s_axilite port=A
//	#pragma HLS ARRAY_RESHAPE variable=B complete dim=1
//	#pragma HLS ARRAY_RESHAPE variable=A complete dim=2
//	for(int i=0;i<4;i++)
//		for(int j=0;j<4;j++){
//			#pragma HLS PIPELINE II=1
//			C[i][j]=0;
//			// 循环乘四次，并进行相加
//			for(int k=0;k<4;k++){
//				C[i][j]=C[i][j]+A[i][k]*B[k][j];
//			}
//		}
//}


void MatrixMultiple(double A[100], double B[100],double C[100]){
	#pragma HLS INTERFACE s_axilite port=return
	#pragma HLS INTERFACE s_axilite port=B
	#pragma HLS INTERFACE s_axilite port=C
	#pragma HLS INTERFACE s_axilite port=A
//	#pragma HLS ARRAY_RESHAPE variable=B complete dim=1
//	#pragma HLS ARRAY_RESHAPE variable=A complete dim=2
//	double a = 1, b = 2, c = 1;
//		  VertexCurveABC abc;
//		  abc.parameters[0][0] = 0;
//		  abc.parameters[1][0] = 0;
//		  abc.parameters[2][0] = 0;

	  // 构建问题
		  Problem problem;

		  // 构建100条边
		  for (int i = 0; i < 100; i++) {
//			verticies_[0]=&abc;
			// problem.edge_curves[i].verticies_[0].parameters[0][0] = 0;
			// problem.edge_curves[i].verticies_[0].parameters[1][0] = 0;
			// problem.edge_curves[i].verticies_[0].parameters[2][0] = 0;
			double info[1][1];
			info[0][0] = 1;
			problem.edge_curves[i].SetInformation(info);
//			double x = i / 100.0;
//			problem.edge_curves[i].x_ = x;
//			problem.edge_curves[i].y_ = std::exp(a * x * x + b * x + c);
			problem.edge_curves[i].x_ = A[i];
			problem.edge_curves[i].y_ = B[i];
			problem.edge_curves_size++;
		  }

		  	  // problem.CalculateResidual();
		  	  // std::terminate();
		  	  // std::cout << VAR(problem.Chi2()) << std::endl;
		  	  // problem.MakeHessian();
		  	  problem.Solve();
		  	  C[0]=verticies_[0].parameters[0][0];
		  	  C[1]=verticies_[0].parameters[1][0];
		  	  C[2]=verticies_[0].parameters[2][0];
//		  	  MATRIXDEBUG(abc.parameters);
}


