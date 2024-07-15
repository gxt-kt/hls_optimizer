#include "problem.h"

// 参考 https://blog.csdn.net/qq_45364953/article/details/127641923
//

VertexCurveABC verticies_[10];

VertexInverseDepth v_points[30];

VertexPose v_poses[5];

void MatrixMultiple(unsigned int A_[100], unsigned int B_[100],
                    unsigned int C_[100], unsigned int A_a_[100],
                    unsigned int B_b_[100]) {
#pragma HLS INTERFACE s_axilite port = return
#pragma HLS INTERFACE s_axilite port = A_
#pragma HLS INTERFACE s_axilite port = B_
#pragma HLS INTERFACE s_axilite port = C_
#pragma HLS INTERFACE s_axilite port = A_a_
#pragma HLS INTERFACE s_axilite port = B_b_

  float A[100] = {};
  float B[100] = {};
  float C[100] = {};
  float A_a[100] = {};
  float B_b[100] = {};
  memcpy(A, A_, 100 * sizeof(float));
  memcpy(B, B_, 100 * sizeof(float));

  // 构建问题
  Problem problem;

  // 构建100条边
  for (int i = 0; i < 100; i++) {
    float info[1][1];
    info[0][0] = 1;
    problem.edge_curves[i].SetInformation(info);
    problem.edge_curves[i].x_ = A[i];
    problem.edge_curves[i].y_ = B[i];
    problem.edge_curves_size++;
  }

  problem.Solve();

  C[0] = verticies_[0].parameters(0, 0);
  C[1] = verticies_[0].parameters(1, 0);
  C[2] = verticies_[0].parameters(2, 0);
  for (int i = 0; i < 100; i++) {
    A_a[i] = A[i];
  }
  for (int i = 0; i < 100; i++) {
    B_b[i] = B[i] + 2;
  }
  memcpy(C_, C, 100 * sizeof(float));
  memcpy(A_a_, A_a, 100 * sizeof(float));
  memcpy(B_b_, B_b, 100 * sizeof(float));
}
