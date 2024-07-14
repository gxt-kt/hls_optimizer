#include "matrix.h"
#include <iostream>


float my_exp(float x) {
  float result = 1.0;
  float term = 1.0;
  int i = 1;

  while (fabs(term) > 1e-10) {
    term *= x / i;
    result += term;
    i++;
  }

  return result;
}

int main() {
  //	ap_int<8> a[4][4]={{0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15}};
  //	ap_int<8> b[4][4]={{0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15}};
  //	ap_int<16> c[4][4];
  //	MatrixMultiple(a,b,c);
  //	for(int i=0;i<4;i++) {
  //		for(int j=0;j<4;j++) {
  //			std::cout << "C["<<i<<","<<j<<"]="<<c[i][j] << "  ";
  //		}
  //		std::cout << std::endl;
  //	}

  //	CurveFittingEdge edge;
  //	double info[1][1];
  //	info[0][0]=1;
  //	edge.ComputeResidual();
  //	edge.ComputeJacobians();
  //	edge.SetInformation(info);
  //	std::cout << "edge.Chi2()=" << edge.Chi2() << std::endl;;

  //	Test();
  //	TestStream();
  //	TestMyMatrixMultiple();
  //	TestGemm();
  //	edge.Fun();
  //	edge.Print();

  //	  double a = 1, b = 2, c = 1;
  //	  VertexCurveABC abc;
  //	  abc.parameters[0][0] = 0;
  //	  abc.parameters[1][0] = 0;
  //	  abc.parameters[2][0] = 0;
  //
  //	  // 构建问题
  //	  Problem problem;
  //
  //	  // 构建100条边
  //	  for (int i = 0; i < 100; i++) {
  //	    problem.edge_curves[i].verticies_[0]=&abc;
  //	    // problem.edge_curves[i].verticies_[0].parameters[0][0] = 0;
  //	    // problem.edge_curves[i].verticies_[0].parameters[1][0] = 0;
  //	    // problem.edge_curves[i].verticies_[0].parameters[2][0] = 0;
  //	    double info[1][1];
  //	    info[0][0] = 1;
  //	    problem.edge_curves[i].SetInformation(info);
  //	    double x = i / 100.0;
  //	    problem.edge_curves[i].x_ = x;
  //	    problem.edge_curves[i].y_ = std::exp(a * x * x + b * x + c);
  //	    problem.edge_curves_size++;
  //	  }
  //
  //	  // problem.CalculateResidual();
  //	  // std::terminate();
  //	  // std::cout << VAR(problem.Chi2()) << std::endl;
  //	  // problem.MakeHessian();
  //	  problem.Solve();
  //	  MATRIXDEBUG(abc.parameters);
  //	double A[100]={};
  //	double B[100]={};
  //	double C[100]={};
  //	double a = 1, b = 2, c = 1;
  //	for(int i=0;i <100;i++) {
  //		A[i]=i/100.0;
  //		double x=A[i];
  //		B[i]=my_exp(a * x * x + b * x + c);
  //	}
  //	MatrixMultiple(A,B,C);
  //	for(int i=0;i<3;i++) {
  //		std::cout << C[i]  << std::endl;
  //	}
  {
    float A[100] = {};
    float B[100] = {};
    float C[100] = {};
    float A_a[100] = {};
    float B_b[100] = {};
    float a = 1, b = 2, c = 1;
    for (int i = 0; i < 100; i++) {
      A[i] = i / 100.0;
      float x = A[i];
      B[i] = my_exp(a * x * x + b * x + c);
    }
    unsigned int A_[100] = {};
    unsigned int B_[100] = {};
    unsigned int C_[100] = {};
    unsigned int A_a_[100] = {};
    unsigned int B_b_[100] = {};
    memcpy(A_, A, 100 * sizeof(float));
    memcpy(B_, B, 100 * sizeof(float));
    memcpy(C_, C, 100 * sizeof(float));
    memcpy(A_a_, A_a, 100 * sizeof(float));
    memcpy(B_b_, B_b, 100 * sizeof(float));
    MatrixMultiple(A_, B_, C_, A_a_, B_b_);
    memcpy(C, C_, 100 * sizeof(float));
    memcpy(A_a, A_a_, 100 * sizeof(float));
    memcpy(B_b, B_b_, 100 * sizeof(float));
    for (int i = 0; i < 3; i++) {
      std::cout << C[i] << std::endl;
    }
    for (int i = 0; i < 3; i++) {
      std::cout << A[i] << std::endl;
      std::cout << A_a[i] << std::endl;
      std::cout << B[i] << std::endl;
      std::cout << B_b[i] << std::endl;
    }
  }
}
