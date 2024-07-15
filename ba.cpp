#include <exception>
#include <iostream>
#include <unordered_map>

#include "edge.h"
#include "matrix.h"
#include "problem.h"


VertexInverseDepth v_points[30];

VertexPose v_poses[5];


/*
 * Frame : 保存每帧的姿态和观测
 */
// struct Frame {
//   Frame(const Matrix<float,3,3>& R, const Matrix<float,3,1>& t) : Rwc(R),
//   qwc(R), twc(t){}; Matrix<float,3,3> Rwc; Quaternion<float> qwc;
//   Matrix<float,3,1> twc;
//
//   std::unordered_map<int, Matrix<float,3,1>> featurePerId;  //
//   该帧观测到的特征以及特征id
// };

int main() {
  {
    VertexPose tmp;
    tmp.parameters[0] = 1;
    tmp.parameters[1] = 2;
    tmp.parameters[2] = 3;
    tmp.parameters[3] = 0.1;
    tmp.parameters[4] = 0.2;
    tmp.parameters[5] = 0.3;
    tmp.parameters[6] = sqrt(1 - (0.01 + 0.04 + 0.09));
    Matrix<float, 7, 1> delta;
    delta[0] = 1;
    delta[1] = 2;
    delta[2] = 3;
    delta[3] = 0.1;
    delta[4] = 0.2;
    delta[5] = 0.1;
    delta[6] = sqrt(1 - (0.01 + 0.04 + 0.01));
    tmp.Plus(delta);
    std::cout << tmp.parameters << std::endl;
    // std::terminate();
  }

  v_poses[0].parameters[0] = 0;
  v_poses[0].parameters[1] = 0;
  v_poses[0].parameters[2] = 0;
  v_poses[0].parameters[3] = 0;
  v_poses[0].parameters[4] = 0;
  v_poses[0].parameters[5] = 0;
  v_poses[0].parameters[6] = 1;

  v_poses[1].parameters[0] = -1.0718;
  v_poses[1].parameters[1] = 4;
  v_poses[1].parameters[2] = 0.866025;
  v_poses[1].parameters[3] = 0;
  v_poses[1].parameters[4] = 0;
  v_poses[1].parameters[5] = 0.258819;
  v_poses[1].parameters[6] = 0.965926;

  v_poses[2].parameters[0] = -4;
  v_poses[2].parameters[1] = 6.9282;
  v_poses[2].parameters[2] = 0.866025;
  v_poses[2].parameters[3] = 0;
  v_poses[2].parameters[4] = 0;
  v_poses[2].parameters[5] = 0.5;
  v_poses[2].parameters[6] = 0.866025;

  v_points[0].parameters[0] = -3.31974;
  v_points[0].parameters[1] = 3.13289;
  v_points[0].parameters[2] = 4.75876;
  v_points[1].parameters[0] = 3.73289;
  v_points[1].parameters[1] = -0.591593;
  v_points[1].parameters[2] = 6.61199;
  v_points[2].parameters[0] = -0.541924;
  v_points[2].parameters[1] = 2.75942;
  v_points[2].parameters[2] = 7.09139;
  v_points[3].parameters[0] = 3.93844;
  v_points[3].parameters[1] = -2.68862;
  v_points[3].parameters[2] = 5.33939;
  v_points[4].parameters[0] = 3.92643;
  v_points[4].parameters[1] = 2.96511;
  v_points[4].parameters[2] = 6.69061;
  v_points[5].parameters[0] = 1.06045;
  v_points[5].parameters[1] = -2.78982;
  v_points[5].parameters[2] = 7.79545;
  v_points[6].parameters[0] = -2.61078;
  v_points[6].parameters[1] = 2.01196;
  v_points[6].parameters[2] = 7.29967;
  v_points[7].parameters[0] = 0.479352;
  v_points[7].parameters[1] = 0.0333801;
  v_points[7].parameters[2] = 5.24006;
  v_points[8].parameters[0] = -1.36016;
  v_points[8].parameters[1] = -1.4589;
  v_points[8].parameters[2] = 7.91143;
  v_points[9].parameters[0] = -3.70806;
  v_points[9].parameters[1] = 2.64461;
  v_points[9].parameters[2] = 7.61457;
  v_points[10].parameters[0] = -1.94;
  v_points[10].parameters[1] = 0.884393;
  v_points[10].parameters[2] = 4.78602;
  v_points[11].parameters[0] = 3.95182;
  v_points[11].parameters[1] = 1.59009;
  v_points[11].parameters[2] = 7.3894;
  v_points[12].parameters[0] = 1.58711;
  v_points[12].parameters[1] = -2.34515;
  v_points[12].parameters[2] = 4.72905;
  v_points[13].parameters[0] = 0.770056;
  v_points[13].parameters[1] = 3.23086;
  v_points[13].parameters[2] = 5.32613;
  v_points[14].parameters[0] = -3.62124;
  v_points[14].parameters[1] = 1.40829;
  v_points[14].parameters[2] = 4.07279;
  v_points[15].parameters[0] = 3.14945;
  v_points[15].parameters[1] = 1.81695;
  v_points[15].parameters[2] = 5.28606;
  v_points[16].parameters[0] = 2.32188;
  v_points[16].parameters[1] = -2.93249;
  v_points[16].parameters[2] = 5.52182;
  v_points[17].parameters[0] = -2.32672;
  v_points[17].parameters[1] = 1.07221;
  v_points[17].parameters[2] = 6.55827;
  v_points[18].parameters[0] = 1.00828;
  v_points[18].parameters[1] = 0.716532;
  v_points[18].parameters[2] = 7.86602;
  v_points[19].parameters[0] = -1.39402;
  v_points[19].parameters[1] = 0.114874;
  v_points[19].parameters[2] = 4.5346;
}
