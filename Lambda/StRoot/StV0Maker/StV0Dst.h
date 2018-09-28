#ifndef StV0Dst_def
#define StV0Dst_def

#define MAX_NUM_V0 3000
struct StV0Dst {
  /// 带电粒子的reference multiplicity.
  int nrefmult;
  /// 最初的碰撞顶点三坐标
  double primvertexX;
  double primvertexY;
  double primvertexZ;
  /// what is vpd ???
  double vpdVz;
  /// 磁场
  double magn;
  /// what is nv0 ???
  int nv0;
  /// 事件ID
  long eventid;
  /// Run ID
  long runid;

  // StV0Data v0[1000];
  double v0mass[MAX_NUM_V0];
  double v0pt[MAX_NUM_V0];
  double v0rapidity[MAX_NUM_V0];
  double v0eta[MAX_NUM_V0];
  /// 三坐标和三动量
  double v0x[MAX_NUM_V0];
  double v0y[MAX_NUM_V0];
  double v0z[MAX_NUM_V0];
  double v0px[MAX_NUM_V0];
  double v0py[MAX_NUM_V0];
  double v0pz[MAX_NUM_V0];
  /// what is v0path ???
  double v0xpath[MAX_NUM_V0];
  double v0ypath[MAX_NUM_V0];
  double v0zpath[MAX_NUM_V0];
  double v0pxpath[MAX_NUM_V0];
  double v0pypath[MAX_NUM_V0];
  double v0pzpath[MAX_NUM_V0];

  double v0declen[MAX_NUM_V0];
  double v0declenH[MAX_NUM_V0];
  /// what is v0dca ???
  double v0dca[MAX_NUM_V0];
  /// v0的r dot p
  double v0rdotp[MAX_NUM_V0];
  /// what is sinth ???
  double v0sinth[MAX_NUM_V0];
  /// v0 粒子的theta角
  double v0theta[MAX_NUM_V0];
  /// 粒子1和2之间的最短距离
  double dca1to2[MAX_NUM_V0];

  // dau1******************************************
  /// PID
  int dau1id[MAX_NUM_V0];
  /// what is dau1dca ???
  double dau1dca[MAX_NUM_V0];
  /// 质量
  double dau1mass[MAX_NUM_V0];
  /// eta
  double dau1eta[MAX_NUM_V0];
  /// ???
  int dau1nhitsfit[MAX_NUM_V0];
  int dau1nhitsdedx[MAX_NUM_V0];
  int dau1nhitsposs[MAX_NUM_V0];
  int dau1hftbits[MAX_NUM_V0];
  double dau1nsigma[MAX_NUM_V0];
  /// 横动量和三动量
  double dau1pt[MAX_NUM_V0];
  double dau1px[MAX_NUM_V0];
  double dau1py[MAX_NUM_V0];
  double dau1pz[MAX_NUM_V0];
  /// 能损
  double dau1dEdx[MAX_NUM_V0];
  /// Z的位置
  double dau1Z[MAX_NUM_V0];
  double dau1beta[MAX_NUM_V0];
  double dau1diffinvbeta[MAX_NUM_V0];
  double dau1toflocaly[MAX_NUM_V0];
  double dau1toflocalz[MAX_NUM_V0];
  int dau1PXL1[MAX_NUM_V0];
  int dau1PXL2[MAX_NUM_V0];
  int dau1PXL3[MAX_NUM_V0];
  int dau1IST1[MAX_NUM_V0];
  int dau1IST2[MAX_NUM_V0];
  int dau1SSD1[MAX_NUM_V0];
  int dau1SSD2[MAX_NUM_V0];

  // dau2********************************************
  int dau2id[MAX_NUM_V0];
  double dau2dca[MAX_NUM_V0];
  double dau2mass[MAX_NUM_V0];
  double dau2eta[MAX_NUM_V0];
  int dau2nhitsfit[MAX_NUM_V0];
  int dau2nhitsdedx[MAX_NUM_V0];
  int dau2nhitsposs[MAX_NUM_V0];
  int dau2hftbits[MAX_NUM_V0];
  double dau2nsigma[MAX_NUM_V0];
  double dau2pt[MAX_NUM_V0];
  double dau2px[MAX_NUM_V0];
  double dau2py[MAX_NUM_V0];
  double dau2pz[MAX_NUM_V0];
  double dau2dEdx[MAX_NUM_V0];
  double dau2Z[MAX_NUM_V0];
  double dau2beta[MAX_NUM_V0];
  double dau2diffinvbeta[MAX_NUM_V0];
  double dau2toflocaly[MAX_NUM_V0];
  double dau2toflocalz[MAX_NUM_V0];
  int dau2PXL1[MAX_NUM_V0];
  int dau2PXL2[MAX_NUM_V0];
  int dau2PXL3[MAX_NUM_V0];
  int dau2IST1[MAX_NUM_V0];
  int dau2IST2[MAX_NUM_V0];
  int dau2SSD1[MAX_NUM_V0];
  int dau2SSD2[MAX_NUM_V0];
};

#endif
