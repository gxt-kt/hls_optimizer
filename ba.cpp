#include <exception>
#include <iostream>
#include <unordered_map>

#include "edge.h"
#include "matrix.h"
#include "problem.h"

VertexInverseDepth v_points[30];

VertexPose v_poses[5];

EdgeReprojection e_reproject[50];

const int e_reproject_size = 40;

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
  // 测试用的
  {
    VertexPose tmp;
    tmp.parameters[0] = 1;
    tmp.parameters[1] = 2;
    tmp.parameters[2] = 3;
    tmp.parameters[3] = 0.1;
    tmp.parameters[4] = 0.2;
    tmp.parameters[5] = 0.3;
    tmp.parameters[6] = sqrt(1 - (0.01 + 0.04 + 0.09));
    Matrix<float, 6, 1> delta;
    delta[0] = 1;
    delta[1] = 2;
    delta[2] = 3;
    delta[3] = 0.1;
    delta[4] = 0.2;
    delta[5] = 0.1;
    // delta[6] = sqrt(1 - (0.01 + 0.04 + 0.01));
    tmp.Plus(delta);
    std::cout << tmp.parameters << std::endl;
    // std::terminate();
  }

  // 位姿总共三个顶点，添加完成
  {
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
  }

  // 特征点
  {
    v_points[0].parameters[0] = 0.28417;
    e_reproject[0].pts_i_[0] = -0.698224;
    e_reproject[0].pts_i_[1] = 0.659816;
    e_reproject[0].pts_i_[2] = 1.000000;
    e_reproject[0].pts_j_[0] = -0.611221;
    e_reproject[0].pts_j_[1] = 0.097160;
    e_reproject[0].pts_j_[2] = 1.000000;
    e_reproject[0].v_idx0 = 0;
    e_reproject[0].v_idx1 = 0;
    e_reproject[0].v_idx2 = 1;
    e_reproject[1].pts_i_[0] = -0.698224;
    e_reproject[1].pts_i_[1] = 0.659816;
    e_reproject[1].pts_i_[2] = 1.000000;
    e_reproject[1].pts_j_[0] = -0.756972;
    e_reproject[1].pts_j_[1] = -0.638630;
    e_reproject[1].pts_j_[2] = 1.000000;
    e_reproject[1].v_idx0 = 0;
    e_reproject[1].v_idx1 = 0;
    e_reproject[1].v_idx2 = 2;
    v_points[1].parameters[0] = 0.161173;
    e_reproject[2].pts_i_[0] = 0.562903;
    e_reproject[2].pts_i_[1] = -0.090161;
    e_reproject[2].pts_i_[2] = 1.000000;
    e_reproject[2].pts_j_[0] = 0.324683;
    e_reproject[2].pts_j_[1] = -1.108658;
    e_reproject[2].pts_j_[2] = 1.000000;
    e_reproject[2].v_idx0 = 1;
    e_reproject[2].v_idx1 = 0;
    e_reproject[2].v_idx2 = 1;
    e_reproject[3].pts_i_[0] = 0.562903;
    e_reproject[3].pts_i_[1] = -0.090161;
    e_reproject[3].pts_i_[2] = 1.000000;
    e_reproject[3].pts_j_[0] = -0.460380;
    e_reproject[3].pts_j_[1] = -1.820738;
    e_reproject[3].pts_j_[2] = 1.000000;
    e_reproject[3].v_idx0 = 1;
    e_reproject[3].v_idx1 = 0;
    e_reproject[3].v_idx2 = 2;
    v_points[2].parameters[0] = 0.114898;
    e_reproject[4].pts_i_[0] = -0.077052;
    e_reproject[4].pts_i_[1] = 0.389698;
    e_reproject[4].pts_i_[2] = 1.000000;
    e_reproject[4].pts_j_[0] = -0.027186;
    e_reproject[4].pts_j_[1] = -0.214825;
    e_reproject[4].pts_j_[2] = 1.000000;
    e_reproject[4].v_idx0 = 2;
    e_reproject[4].v_idx1 = 0;
    e_reproject[4].v_idx2 = 1;
    e_reproject[5].pts_i_[0] = -0.077052;
    e_reproject[5].pts_i_[1] = 0.389698;
    e_reproject[5].pts_i_[2] = 1.000000;
    e_reproject[5].pts_j_[0] = -0.303893;
    e_reproject[5].pts_j_[1] = -0.816691;
    e_reproject[5].pts_j_[2] = 1.000000;
    e_reproject[5].v_idx0 = 2;
    e_reproject[5].v_idx1 = 0;
    e_reproject[5].v_idx2 = 2;
    v_points[3].parameters[0] = 0.174241;
    e_reproject[6].pts_i_[0] = 0.737074;
    e_reproject[6].pts_i_[1] = -0.503583;
    e_reproject[6].pts_i_[2] = 1.000000;
    e_reproject[6].pts_j_[0] = 0.222552;
    e_reproject[6].pts_j_[1] = -1.855382;
    e_reproject[6].pts_j_[2] = 1.000000;
    e_reproject[6].v_idx0 = 3;
    e_reproject[6].v_idx1 = 0;
    e_reproject[6].v_idx2 = 1;
    e_reproject[7].pts_i_[0] = 0.737074;
    e_reproject[7].pts_i_[1] = -0.503583;
    e_reproject[7].pts_i_[2] = 1.000000;
    e_reproject[7].pts_j_[0] = -0.974847;
    e_reproject[7].pts_j_[1] = -2.612872;
    e_reproject[7].pts_j_[2] = 1.000000;
    e_reproject[7].v_idx0 = 3;
    e_reproject[7].v_idx1 = 0;
    e_reproject[7].v_idx2 = 2;
    v_points[4].parameters[0] = 0.123863;
    e_reproject[8].pts_i_[0] = 0.587402;
    e_reproject[8].pts_i_[1] = 0.441654;
    e_reproject[8].pts_i_[2] = 1.000000;
    e_reproject[8].pts_j_[0] = 0.654956;
    e_reproject[8].pts_j_[1] = -0.584492;
    e_reproject[8].pts_j_[2] = 1.000000;
    e_reproject[8].v_idx0 = 4;
    e_reproject[8].v_idx1 = 0;
    e_reproject[8].v_idx2 = 1;
    e_reproject[9].pts_i_[0] = 0.587402;
    e_reproject[9].pts_i_[1] = 0.441654;
    e_reproject[9].pts_i_[2] = 1.000000;
    e_reproject[9].pts_j_[0] = 0.090912;
    e_reproject[9].pts_j_[1] = -1.519212;
    e_reproject[9].pts_j_[2] = 1.000000;
    e_reproject[9].v_idx0 = 4;
    e_reproject[9].v_idx1 = 0;
    e_reproject[9].v_idx2 = 2;
    v_points[5].parameters[0] = 0.127425;
    e_reproject[10].pts_i_[0] = 0.135846;
    e_reproject[10].pts_i_[1] = -0.359391;
    e_reproject[10].pts_i_[2] = 1.000000;
    e_reproject[10].pts_j_[0] = -0.222120;
    e_reproject[10].pts_j_[1] = -1.000498;
    e_reproject[10].pts_j_[2] = 1.000000;
    e_reproject[10].v_idx0 = 5;
    e_reproject[10].v_idx1 = 0;
    e_reproject[10].v_idx2 = 1;
    e_reproject[11].pts_i_[0] = 0.135846;
    e_reproject[11].pts_i_[1] = -0.359391;
    e_reproject[11].pts_i_[2] = 1.000000;
    e_reproject[11].pts_j_[0] = -0.848144;
    e_reproject[11].pts_j_[1] = -1.332617;
    e_reproject[11].pts_j_[2] = 1.000000;
    e_reproject[11].v_idx0 = 5;
    e_reproject[11].v_idx1 = 0;
    e_reproject[11].v_idx2 = 2;
    v_points[6].parameters[0] = 0.156359;
    e_reproject[12].pts_i_[0] = -0.356537;
    e_reproject[12].pts_i_[1] = 0.274822;
    e_reproject[12].pts_i_[2] = 1.000000;
    e_reproject[12].pts_j_[0] = -0.361557;
    e_reproject[12].pts_j_[1] = -0.149174;
    e_reproject[12].pts_j_[2] = 1.000000;
    e_reproject[12].v_idx0 = 6;
    e_reproject[12].v_idx1 = 0;
    e_reproject[12].v_idx2 = 1;
    e_reproject[13].pts_i_[0] = -0.356537;
    e_reproject[13].pts_i_[1] = 0.274822;
    e_reproject[13].pts_i_[2] = 1.000000;
    e_reproject[13].pts_j_[0] = -0.553713;
    e_reproject[13].pts_j_[1] = -0.570962;
    e_reproject[13].pts_j_[2] = 1.000000;
    e_reproject[13].v_idx0 = 6;
    e_reproject[13].v_idx1 = 0;
    e_reproject[13].v_idx2 = 2;
    v_points[7].parameters[0] = 0.140633;
    e_reproject[14].pts_i_[0] = 0.090997;
    e_reproject[14].pts_i_[1] = 0.004534;
    e_reproject[14].pts_i_[2] = 1.000000;
    e_reproject[14].pts_j_[0] = -0.146570;
    e_reproject[14].pts_j_[1] = -0.963064;
    e_reproject[14].pts_j_[2] = 1.000000;
    e_reproject[14].v_idx0 = 7;
    e_reproject[14].v_idx1 = 0;
    e_reproject[14].v_idx2 = 1;
    e_reproject[15].pts_i_[0] = 0.090997;
    e_reproject[15].pts_i_[1] = 0.004534;
    e_reproject[15].pts_i_[2] = 1.000000;
    e_reproject[15].pts_j_[0] = -0.850387;
    e_reproject[15].pts_j_[1] = -1.673371;
    e_reproject[15].pts_j_[2] = 1.000000;
    e_reproject[15].v_idx0 = 7;
    e_reproject[15].v_idx1 = 0;
    e_reproject[15].v_idx2 = 2;
    v_points[8].parameters[0] = 0.159974;
    e_reproject[16].pts_i_[0] = -0.172521;
    e_reproject[16].pts_i_[1] = -0.184307;
    e_reproject[16].pts_i_[2] = 1.000000;
    e_reproject[16].pts_j_[0] = -0.424206;
    e_reproject[16].pts_j_[1] = -0.648827;
    e_reproject[16].pts_j_[2] = 1.000000;
    e_reproject[16].v_idx0 = 8;
    e_reproject[16].v_idx1 = 0;
    e_reproject[16].v_idx2 = 1;
    e_reproject[17].pts_i_[0] = -0.172521;
    e_reproject[17].pts_i_[1] = -0.184307;
    e_reproject[17].pts_i_[2] = 1.000000;
    e_reproject[17].pts_j_[0] = -0.843437;
    e_reproject[17].pts_j_[1] = -0.919438;
    e_reproject[17].pts_j_[2] = 1.000000;
    e_reproject[17].v_idx0 = 8;
    e_reproject[17].v_idx1 = 0;
    e_reproject[17].v_idx2 = 2;
    v_points[9].parameters[0] = 0.144373;
    e_reproject[18].pts_i_[0] = -0.488228;
    e_reproject[18].pts_i_[1] = 0.347568;
    e_reproject[18].pts_i_[2] = 1.000000;
    e_reproject[18].pts_j_[0] = -0.440352;
    e_reproject[18].pts_j_[1] = 0.020390;
    e_reproject[18].pts_j_[2] = 1.000000;
    e_reproject[18].v_idx0 = 9;
    e_reproject[18].v_idx1 = 0;
    e_reproject[18].v_idx2 = 1;
    e_reproject[19].pts_i_[0] = -0.488228;
    e_reproject[19].pts_i_[1] = 0.347568;
    e_reproject[19].pts_i_[2] = 1.000000;
    e_reproject[19].pts_j_[0] = -0.528224;
    e_reproject[19].pts_j_[1] = -0.356169;
    e_reproject[19].pts_j_[2] = 1.000000;
    e_reproject[19].v_idx0 = 9;
    e_reproject[19].v_idx1 = 0;
    e_reproject[19].v_idx2 = 2;
    v_points[10].parameters[0] = 0.205688;
    e_reproject[20].pts_i_[0] = -0.406218;
    e_reproject[20].pts_i_[1] = 0.185198;
    e_reproject[20].pts_i_[2] = 1.000000;
    e_reproject[20].pts_j_[0] = -0.588817;
    e_reproject[20].pts_j_[1] = -0.578023;
    e_reproject[20].pts_j_[2] = 1.000000;
    e_reproject[20].v_idx0 = 10;
    e_reproject[20].v_idx1 = 0;
    e_reproject[20].v_idx2 = 1;
    e_reproject[21].pts_i_[0] = -0.406218;
    e_reproject[21].pts_i_[1] = 0.185198;
    e_reproject[21].pts_i_[2] = 1.000000;
    e_reproject[21].pts_j_[0] = -1.072671;
    e_reproject[21].pts_j_[1] = -1.227675;
    e_reproject[21].pts_j_[2] = 1.000000;
    e_reproject[21].v_idx0 = 10;
    e_reproject[21].v_idx1 = 0;
    e_reproject[21].v_idx2 = 2;
    v_points[11].parameters[0] = 0.112831;
    e_reproject[22].pts_i_[0] = 0.534608;
    e_reproject[22].pts_i_[1] = 0.215162;
    e_reproject[22].pts_i_[2] = 1.000000;
    e_reproject[22].pts_j_[0] = 0.481905;
    e_reproject[22].pts_j_[1] = -0.705979;
    e_reproject[22].pts_j_[2] = 1.000000;
    e_reproject[22].v_idx0 = 11;
    e_reproject[22].v_idx1 = 0;
    e_reproject[22].v_idx2 = 1;
    e_reproject[23].pts_i_[0] = 0.534608;
    e_reproject[23].pts_i_[1] = 0.215162;
    e_reproject[23].pts_i_[2] = 1.000000;
    e_reproject[23].pts_j_[0] = -0.100814;
    e_reproject[23].pts_j_[1] = -1.465065;
    e_reproject[23].pts_j_[2] = 1.000000;
    e_reproject[23].v_idx0 = 11;
    e_reproject[23].v_idx1 = 0;
    e_reproject[23].v_idx2 = 2;
    v_points[12].parameters[0] = 0.207199;
    e_reproject[24].pts_i_[0] = 0.335303;
    e_reproject[24].pts_i_[1] = -0.496942;
    e_reproject[24].pts_i_[2] = 1.000000;
    e_reproject[24].pts_j_[0] = -0.223372;
    e_reproject[24].pts_j_[1] = -1.766851;
    e_reproject[24].pts_j_[2] = 1.000000;
    e_reproject[24].v_idx0 = 12;
    e_reproject[24].v_idx1 = 0;
    e_reproject[24].v_idx2 = 1;
    e_reproject[25].pts_i_[0] = 0.335303;
    e_reproject[25].pts_i_[1] = -0.496942;
    e_reproject[25].pts_i_[2] = 1.000000;
    e_reproject[25].pts_j_[0] = -1.355604;
    e_reproject[25].pts_j_[1] = -2.452312;
    e_reproject[25].pts_j_[2] = 1.000000;
    e_reproject[25].v_idx0 = 12;
    e_reproject[25].v_idx1 = 0;
    e_reproject[25].v_idx2 = 2;
    v_points[13].parameters[0] = 0.225542;
    e_reproject[26].pts_i_[0] = 0.143386;
    e_reproject[26].pts_i_[1] = 0.607676;
    e_reproject[26].pts_i_[2] = 1.000000;
    e_reproject[26].pts_j_[0] = 0.270200;
    e_reproject[26].pts_j_[1] = -0.358021;
    e_reproject[26].pts_j_[2] = 1.000000;
    e_reproject[26].v_idx0 = 13;
    e_reproject[26].v_idx1 = 0;
    e_reproject[26].v_idx2 = 1;
    e_reproject[27].pts_i_[0] = 0.143386;
    e_reproject[27].pts_i_[1] = 0.607676;
    e_reproject[27].pts_i_[2] = 1.000000;
    e_reproject[27].pts_j_[0] = -0.182491;
    e_reproject[27].pts_j_[1] = -1.339294;
    e_reproject[27].pts_j_[2] = 1.000000;
    e_reproject[27].v_idx0 = 13;
    e_reproject[27].v_idx1 = 0;
    e_reproject[27].v_idx2 = 2;
    v_points[14].parameters[0] = 0.260145;
    e_reproject[28].pts_i_[0] = -0.889074;
    e_reproject[28].pts_i_[1] = 0.346496;
    e_reproject[28].pts_i_[2] = 1.000000;
    e_reproject[28].pts_j_[0] = -1.094977;
    e_reproject[28].pts_j_[1] = -0.302561;
    e_reproject[28].pts_j_[2] = 1.000000;
    e_reproject[28].v_idx0 = 14;
    e_reproject[28].v_idx1 = 0;
    e_reproject[28].v_idx2 = 1;
    e_reproject[29].pts_i_[0] = -0.889074;
    e_reproject[29].pts_i_[1] = 0.346496;
    e_reproject[29].pts_i_[2] = 1.000000;
    e_reproject[29].pts_j_[0] = -1.430994;
    e_reproject[29].pts_j_[1] = -0.961874;
    e_reproject[29].pts_j_[2] = 1.000000;
    e_reproject[29].v_idx0 = 14;
    e_reproject[29].v_idx1 = 0;
    e_reproject[29].v_idx2 = 2;
    v_points[15].parameters[0] = 0.155015;
    e_reproject[30].pts_i_[0] = 0.595237;
    e_reproject[30].pts_i_[1] = 0.343305;
    e_reproject[30].pts_i_[2] = 1.000000;
    e_reproject[30].pts_j_[0] = 0.579474;
    e_reproject[30].pts_j_[1] = -0.905888;
    e_reproject[30].pts_j_[2] = 1.000000;
    e_reproject[30].v_idx0 = 15;
    e_reproject[30].v_idx1 = 0;
    e_reproject[30].v_idx2 = 1;
    e_reproject[31].pts_i_[0] = 0.595237;
    e_reproject[31].pts_i_[1] = 0.343305;
    e_reproject[31].pts_i_[2] = 1.000000;
    e_reproject[31].pts_j_[0] = -0.192231;
    e_reproject[31].pts_j_[1] = -1.977474;
    e_reproject[31].pts_j_[2] = 1.000000;
    e_reproject[31].v_idx0 = 15;
    e_reproject[31].v_idx1 = 0;
    e_reproject[31].v_idx2 = 2;
    v_points[16].parameters[0] = 0.163754;
    e_reproject[32].pts_i_[0] = 0.419542;
    e_reproject[32].pts_i_[1] = -0.532542;
    e_reproject[32].pts_i_[2] = 1.000000;
    e_reproject[32].pts_j_[0] = -0.112762;
    e_reproject[32].pts_j_[1] = -1.652898;
    e_reproject[32].pts_j_[2] = 1.000000;
    e_reproject[32].v_idx0 = 16;
    e_reproject[32].v_idx1 = 0;
    e_reproject[32].v_idx2 = 1;
    e_reproject[33].pts_i_[0] = 0.419542;
    e_reproject[33].pts_i_[1] = -0.532542;
    e_reproject[33].pts_i_[2] = 1.000000;
    e_reproject[33].pts_j_[0] = -1.153765;
    e_reproject[33].pts_j_[1] = -2.234775;
    e_reproject[33].pts_j_[2] = 1.000000;
    e_reproject[33].v_idx0 = 16;
    e_reproject[33].v_idx1 = 0;
    e_reproject[33].v_idx2 = 2;
    v_points[17].parameters[0] = 0.169557;
    e_reproject[34].pts_i_[0] = -0.353633;
    e_reproject[34].pts_i_[1] = 0.164353;
    e_reproject[34].pts_i_[2] = 1.000000;
    e_reproject[34].pts_j_[0] = -0.448508;
    e_reproject[34].pts_j_[1] = -0.334578;
    e_reproject[34].pts_j_[2] = 1.000000;
    e_reproject[34].v_idx0 = 17;
    e_reproject[34].v_idx1 = 0;
    e_reproject[34].v_idx2 = 1;
    e_reproject[35].pts_i_[0] = -0.353633;
    e_reproject[35].pts_i_[1] = 0.164353;
    e_reproject[35].pts_i_[2] = 1.000000;
    e_reproject[35].pts_j_[0] = -0.742471;
    e_reproject[35].pts_j_[1] = -0.767675;
    e_reproject[35].pts_j_[2] = 1.000000;
    e_reproject[35].v_idx0 = 17;
    e_reproject[35].v_idx1 = 0;
    e_reproject[35].v_idx2 = 2;
    v_points[18].parameters[0] = 0.119037;
    e_reproject[36].pts_i_[0] = 0.130487;
    e_reproject[36].pts_i_[1] = 0.090459;
    e_reproject[36].pts_i_[2] = 1.000000;
    e_reproject[36].pts_j_[0] = 0.023189;
    e_reproject[36].pts_j_[1] = -0.555872;
    e_reproject[36].pts_j_[2] = 1.000000;
    e_reproject[36].v_idx0 = 18;
    e_reproject[36].v_idx1 = 0;
    e_reproject[36].v_idx2 = 1;
    e_reproject[37].pts_i_[0] = 0.130487;
    e_reproject[37].pts_i_[1] = 0.090459;
    e_reproject[37].pts_i_[2] = 1.000000;
    e_reproject[37].pts_j_[0] = -0.411669;
    e_reproject[37].pts_j_[1] = -1.063408;
    e_reproject[37].pts_j_[2] = 1.000000;
    e_reproject[37].v_idx0 = 18;
    e_reproject[37].v_idx1 = 0;
    e_reproject[37].v_idx2 = 2;
    v_points[19].parameters[0] = 0.254659;
    e_reproject[38].pts_i_[0] = -0.307829;
    e_reproject[38].pts_i_[1] = 0.024030;
    e_reproject[38].pts_i_[2] = 1.000000;
    e_reproject[38].pts_j_[0] = -0.606387;
    e_reproject[38].pts_j_[1] = -0.871912;
    e_reproject[38].pts_j_[2] = 1.000000;
    e_reproject[38].v_idx0 = 19;
    e_reproject[38].v_idx1 = 0;
    e_reproject[38].v_idx2 = 1;
    e_reproject[39].pts_i_[0] = -0.307829;
    e_reproject[39].pts_i_[1] = 0.024030;
    e_reproject[39].pts_i_[2] = 1.000000;
    e_reproject[39].pts_j_[0] = -1.251924;
    e_reproject[39].pts_j_[1] = -1.545270;
    e_reproject[39].pts_j_[2] = 1.000000;
    e_reproject[39].v_idx0 = 19;
    e_reproject[39].v_idx1 = 0;
    e_reproject[39].v_idx2 = 2;
  }
  {
    for (int i = 0; i < 40; i++) {
      e_reproject[i].information_ =
          decltype(e_reproject[i].information_)::Identity();
    }
  }

  {
    Quaternion<float> qic;
    qic.x() = 0;
    qic.y() = 0;
    qic.z() = 0;
    qic.w() = 1;
    Matrix<float, 3, 1> tic;
    tic[0] = 0;
    tic[1] = 0;
    tic[2] = 0;
    for (int i = 0; i < e_reproject_size; i++) {
      e_reproject[i].qic = qic;
      e_reproject[i].tic = tic;
    }
  }

  BAProblem problem;
  problem.edge_reproject_size = e_reproject_size;
  problem.Solve();
  for (int i = 0; i < 20; i++) {
    std::cout << "i=" << i << "  " << v_points[i].parameters << std::endl;
  }
  // problem.MakeHessian();

  Quaternion<float> a;
  a.x() = 0.1;
  a.y() = 0.2;
  a.z() = 0.3;
  a.w() = std::sqrt(1 - 0.01 - 0.04 - 0.09);
  std::cout << a << std::endl;
  std::cout << a.inverse() << std::endl;
  std::cout << a.toRotationMatrix() << std::endl;
}
