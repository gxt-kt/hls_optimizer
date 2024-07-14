#include <exception>
#include <iostream>
#include <unordered_map>

#include "edge.h"
#include "matrix.h"
#include "problem.h"


/*
 * Frame : 保存每帧的姿态和观测
 */
struct Frame {
  Frame(const Matrix<float,3,3>& R, const Matrix<float,3,1>& t) : Rwc(R), qwc(R), twc(t){};
  Matrix<float,3,3> Rwc;
  Quaternion<float> qwc;
  Matrix<float,3,1> twc;

  std::unordered_map<int, Matrix<float,3,1>> featurePerId;  // 该帧观测到的特征以及特征id
};

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
  
  v_poses[0].parameters[0]=0;
  v_poses[0].parameters[1]=0;
  v_poses[0].parameters[2]=0;
  v_poses[0].parameters[3]=0;
  v_poses[0].parameters[4]=0;
  v_poses[0].parameters[5]=0;
  v_poses[0].parameters[6]=1;

  v_poses[1].parameters[0]=-1.0718;
  v_poses[1].parameters[1]=4;
  v_poses[1].parameters[2]=0.866025;
  v_poses[1].parameters[3]=0;
  v_poses[1].parameters[4]=0;
  v_poses[1].parameters[5]=0.258819;
  v_poses[1].parameters[6]=0.965926;

  v_poses[2].parameters[0]=-4;
  v_poses[2].parameters[1]=6.9282;
  v_poses[2].parameters[2]=0.866025;
  v_poses[2].parameters[3]=0;
  v_poses[2].parameters[4]=0;
  v_poses[2].parameters[5]=0.5;
  v_poses[2].parameters[6]=0.866025;

v_points[0].parameters[0]=-3.31974;
v_points[0].parameters[1]= 3.13289;
v_points[0].parameters[2]= 4.75876;
v_points[0].parameters[0]= 3.73289;
v_points[0].parameters[1]= -0.591593;
v_points[0].parameters[2]= 6.61199;
i:2  -0.541924
  2.75942
  7.09139
i:3   3.93844
-2.68862
 5.33939
i:4  3.92643
2.96511
6.69061
i:5   1.06045
-2.78982
 7.79545
i:6  -2.61078
 2.01196
 7.29967
i:7   0.479352
0.0333801
  5.24006
i:8  -1.36016
 -1.4589
 7.91143
i:9  -3.70806
 2.64461
 7.61457
i:10     -1.94
0.884393
 4.78602
i:11  3.95182
1.59009
 7.3894
i:12   1.58711
-2.34515
 4.72905
i:13  0.770056
 3.23086
 5.32613
i:14  -3.62124
 1.40829
 4.07279
i:15  3.14945
1.81695
5.28606
i:16   2.32188
-2.93249
 5.52182
i:17  -2.32672
 1.07221
 6.55827
i:18   1.00828
0.716532
 7.86602
i:19  -1.39402
0.114874
  4.5346

}
