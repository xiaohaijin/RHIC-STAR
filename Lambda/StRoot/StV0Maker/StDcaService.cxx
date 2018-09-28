#include "StDcaService.h"
#include "SystemOfUnits.h"
#include "TMath.h"
#include "math.h"

// the following two functions will calculate approximately the dca of two
// helices. this is not a code for general purpose. use it when you understand
// the parameters inside the code. for specific, the current parameters are
// tuned for STAR TPC track. Original version from Hui LONG, improved by
// Xianglei ZHU in May, 2008

// master switch for the minimization algorithm used in the case of connecting
// helices: true -> StHelix method. slow, but accurate. just for checking
// purpose. false-> native minimization code. fast, and accurate enough.
bool kStHelixDca = false;

// switch for turn off the minimization in native code. this is OK for Lambda
// reconstruction. but minimization is essential for photon reconstruction.
// since adding minimization causes only minor lag in performance, we should
// keep it as default, even if for lambda reconstruction.
bool kMinimize = true;

// parameters to control the usage of minimization: for lambda recontruction
// 0.3, 0.3 is safe. for photon recontruction, 2, 2 is safe. By 'safe', I mean
// there is no signal loss compared to the sthelix dca method. For photons, this
// is essential to keep all signals (with the payment of performance), since we
// want to use same sign method to estimate the background. there will be
// systematic deviation if the signal loss (due to limited accuracy) is not the
// same for signals and background.
double kShiftConnect = 2;
double kShiftContain = 2;

extern bool kPeriod;
bool kPeriod = false;
// do not consider the period of helices. we should always ask pt>0.15 to avoid
// helices with multiple periods in TPC. those tracks can not be simulated in
// embedding, therefore useless.
// turn on this option will add more background!

double closestDistance(const StPhysicalHelixD& helix1,
                       const StPhysicalHelixD& helix2, double magnet,
                       const StThreeVectorF& pv, StThreeVectorF& xv0,
                       StThreeVectorF& op1, StThreeVectorF& op2) {
  // Hui LONG's function to calculate the Dca of two global tracks.
  // should use inner helix of tracks?
  double PI = TMath::Pi();
  double x[3], p1[3], p2[3], d;
  StTrackHelix protonobject, pionobject;
  StTrackHelix* proton = &protonobject;
  StTrackHelix* pion = &pionobject;

  // fill StTrackHelix
  StPhysicalHelixD helix = helix1;
  int hrot = helix.h();
  double ttheta = helix.phase();
  StThreeVectorF origin = helix.origin();
  double Xc = helix.xcenter();
  double Yc = helix.ycenter();
  double Zc = static_cast<double>(origin.z());
  double r = 1. / helix.curvature();
  double Vz = r * tan(helix.dipAngle());
  double pt = (helix.momentum(magnet * kilogauss)).perp();
  double Pz = (helix.momentum(magnet * kilogauss)).z();

  proton->Xc = Xc;
  proton->Yc = Yc;
  proton->Zc = Zc;
  proton->r = r;
  proton->theta = ttheta;
  proton->h = hrot;
  proton->Vz = Vz;
  proton->pt = pt;
  proton->Pz = Pz;
  proton->Flag = -1;  // max # of periods to scan forward

  helix = helix2;
  hrot = helix.h();
  ttheta = helix.phase();
  origin = helix.origin();
  Xc = helix.xcenter();
  Yc = helix.ycenter();
  Zc = static_cast<double>(origin.z());
  r = 1. / helix.curvature();
  Vz = r * tan(helix.dipAngle());
  pt = (helix.momentum(magnet * kilogauss)).perp();
  Pz = (helix.momentum(magnet * kilogauss)).z();

  pion->Xc = Xc;
  pion->Yc = Yc;
  pion->Zc = Zc;
  pion->r = r;
  pion->theta = ttheta;
  pion->h = hrot;
  pion->Vz = Vz;
  pion->pt = pt;
  pion->Pz = Pz;
  pion->Flag = -1;

  // then hand over to Hui Long's code
  double d_root1, d_root2, d_p_pi, theta, r_p, r_pi, zij, rec_zij, rec_zij2;
  double fi_root1_p, fi_root2_p, fi_root1_pi, fi_root2_pi;
  double reci, reci2, recj, recj2;
  double v0_root1[3], v0_root2[3];
  double fi_p, fi_pi, ti, tj;
  int n1, n2, m1, m2, i, j;
  double alfa = 0;

  int status;

  double R2_TPC = 30000.;  // hard-coded cut

  d_root2 = 500.;
  d_root1 = 500.;

  d_p_pi = sqrt((pion->Xc - proton->Xc) * (pion->Xc - proton->Xc) +
                (pion->Yc - proton->Yc) * (pion->Yc - proton->Yc));

  theta = atan2((pion->Yc - proton->Yc), (pion->Xc - proton->Xc));

  r_p = proton->r;  // radium of proton helix
  r_pi = pion->r;   // radium of pion helix
  if (d_p_pi > (r_p + r_pi) + 2.0 || d_p_pi < fabs(r_p - r_pi) - 2.0)
    return 500.;  // dca cut should be smaller
  /// XZHU: it will be rare that two tracks with d_p_pi ==
  /// r_p+r_pi come from a single particle decay. that
  /// means their momenta should be in the same direction
  /// or in opposite direction. it is normal for photon
  /// conversion to have two daughter electrons in the
  /// same direction. But it will be rare for two v0 decay
  /// daughters in the very opposite direction (that means
  /// the pt of v0 is 0). FIXME: please notice the hard
  /// cut 2.0cm above, you might change it for some
  /// purposes.
  if (d_p_pi >= (r_p + r_pi) || d_p_pi <= fabs(r_p - r_pi)) {
    fi_root1_p = 0.;
    fi_root2_p = 0.;
    if (d_p_pi >= (r_p + r_pi)) {
      fi_root1_pi = PI;
      fi_root2_pi = PI;
      alfa = PI;
      status = 0;
      // cout<<"two circles are seperated!"<<endl;
    } else {
      fi_root1_pi = 0;
      fi_root2_pi = 0;
      if (r_p < r_pi) {  // need to compensate the theta problem for r_p < r_pi
        theta = atan2((-pion->Yc + proton->Yc), (-pion->Xc + proton->Xc));
      }
      alfa = 0;
      status = 1;
      // cout<<"one circle is inside the other!"<<endl;
    }

    if (kStHelixDca) {
      pair<double, double> tmps = helix1.pathLengths(helix2);
      StThreeVectorF ox1 = helix1.at(tmps.first);
      StThreeVectorF ox2 = helix2.at(tmps.second);
      d = static_cast<double>((ox1 - ox2).mag());
      xv0 = (ox1 + ox2) / 2.0;
      op1 = helix1.momentumAt(tmps.first, magnet * kilogauss);
      op2 = helix2.momentumAt(tmps.second, magnet * kilogauss);

      if (tmps.first > 20 || tmps.second > 20) return 500;
      return d;
    }
  } else {
    status = 2;
    // cout<<"two circles are crossing!"<<endl;
    alfa = acos((+r_pi * r_pi - d_p_pi * d_p_pi + r_p * r_p) / 2. / r_p / r_pi);

    fi_root1_p =
        acos((-r_pi * r_pi + d_p_pi * d_p_pi + r_p * r_p) / 2. / r_p / d_p_pi);
    fi_root2_p = -fi_root1_p;
    fi_root1_pi = alfa + fi_root1_p;
    fi_root2_pi = -fi_root1_pi;

    if (kStHelixDca && (d_p_pi >= (r_p + r_pi) - kShiftConnect ||
                        d_p_pi <= fabs(r_p - r_pi) + kShiftContain)) {
      pair<double, double> tmps = helix1.pathLengths(helix2);
      StThreeVectorF ox1 = helix1.at(tmps.first);
      StThreeVectorF ox2 = helix2.at(tmps.second);
      d = static_cast<double>((ox1 - ox2).mag());
      xv0 = (ox1 + ox2) / 2.0;
      op1 = helix1.momentumAt(tmps.first, magnet * kilogauss);
      op2 = helix2.momentumAt(tmps.second, magnet * kilogauss);

      if (tmps.first > 20 || tmps.second > 20) return 500;
      return d;
    }
  }

  v0_root1[2] = -1000.;
  v0_root2[2] = 1000.;

  fi_root1_p = fi_root1_p + theta;
  fi_root2_p = fi_root2_p + theta;
  fi_root1_pi = fi_root1_pi + theta;
  fi_root2_pi = fi_root2_pi + theta;

  v0_root1[0] = proton->Xc + proton->r * cos(fi_root1_p);
  v0_root2[0] = proton->Xc + proton->r * cos(fi_root2_p);
  v0_root1[1] = proton->Yc + proton->r * sin(fi_root1_p);
  v0_root2[1] = proton->Yc + proton->r * sin(fi_root2_p);

  fi_root1_p = (fi_root1_p - proton->theta) / proton->h;
  fi_root2_p = (fi_root2_p - proton->theta) / proton->h;
  fi_root1_pi = (fi_root1_pi - pion->theta) / pion->h;
  fi_root2_pi = (fi_root2_pi - pion->theta) / pion->h;
  fi_root1_p =
      fi_root1_p -
      floor(fi_root1_p / 2 / PI) * 2 *
          PI;  // XZHU: try floor() here, instead of (int). since a wide range
               // of period is tried below, the result should be the same. in
               // fact, fi_root1_p should always be less than 0 (0,-2pi for 1
               // period), since we are looking backward for the position where
               // the helix originates. ADDED: it is not so simple. we should
               // keep the code as it is.
  fi_root2_p = fi_root2_p - floor(fi_root2_p / 2 / PI) * 2 * PI;
  fi_root1_pi = fi_root1_pi - floor(fi_root1_pi / 2 / PI) * 2 * PI;
  fi_root2_pi = fi_root2_pi - floor(fi_root2_pi / 2 / PI) * 2 * PI;

  if ((v0_root1[0] * v0_root1[0] + v0_root1[1] * v0_root1[1]) <
      R2_TPC) {  // XZHU: if the overlapping position is close (25 cm) to the
                 // outer layer of TPC. Even if the dca between two tracks is
                 // small, it is not interesting to us anymore. We find the
                 // right periods here, which is not necessary for particles
                 // with pt > 0.15GeV/c. in recHyperon, the period finding is
                 // not used. ADDED: lots of daughters has small pt. we do not
                 // understand how multiple periods tracks are recorded (split
                 // to multiple tracks or not?). if we scan periods, we could
                 // have multiple signals for a single lambda!!! hope this
                 // effect can be simulated in embedding analysis.
    rec_zij = 1000.;
    n1 = -(static_cast<int>(fabs((proton->Zc - static_cast<double>(pv.z())) /
                                 proton->Vz / 2. / PI)) +
           1);
    n2 = proton->Flag;  // XZHU: it is meanlingless to go forward in a helix to
                        // get the dca point where the particle (or helix)
                        // should orginates from! so i suggest to decrease this
                        // flag to 0, since fi_root1_pi is [-2pi,2pi]. ADDED: it
                        // is not so simple due to the TPC track momentum
                        // resolution. we should keep the code as it is.
    m1 = -(static_cast<int>(fabs((pion->Zc - static_cast<double>(pv.z())) /
                                 pion->Vz / 2. / PI)) +
           1);
    ;
    m2 = pion->Flag;
    if (n1 < -10) n1 = -10;
    if (m1 < -10) m1 = -10;

    if (!kPeriod) n1 = -1;
    if (!kPeriod) m1 = -1;
    if (proton->r * fi_root1_p < 20.) n2 = 0;
    if (pion->r * fi_root1_pi < 20.) m2 = 0;

    reci = fi_root1_p;
    recj = fi_root1_pi;
    for (i = n1; i <= n2; i++) {
      for (j = m1; j <= m2; j++) {
        ti = fi_root1_p + i * 2 * PI;
        tj = fi_root1_pi + j * 2 * PI;
        zij = -proton->Zc - proton->Vz * ti + pion->Zc + pion->Vz * tj;
        if (fabs(zij) < fabs(rec_zij)) {
          rec_zij = zij;
          reci = ti;
          recj = tj;
        }
      }
    }

    fi_root1_p = reci;
    fi_root1_pi = recj;
    d_root1 = rec_zij;

    dcaPToPi(&fi_root1_p, &fi_root1_pi, proton, pion, &d_root1, &v0_root1[0],
             alfa);
  }

  if (status == 2 &&
      (v0_root2[0] * v0_root2[0] + v0_root2[1] * v0_root2[1]) < R2_TPC) {
    rec_zij2 = 1000.;
    n1 = -(static_cast<int>(fabs((proton->Zc - static_cast<double>(pv.z())) /
                                 proton->Vz / 2. / PI)) +
           1);
    n2 = proton->Flag;
    m1 = -(static_cast<int>(fabs((pion->Zc - static_cast<double>(pv.z())) /
                                 pion->Vz / 2. / PI)) +
           1);
    m2 = pion->Flag;
    if (n1 < -10) n1 = -10;
    if (m1 < -10) m1 = -10;

    if (!kPeriod) n1 = -1;
    if (!kPeriod) m1 = -1;
    if (proton->r * fi_root2_p < 20.) n2 = 0;
    if (pion->r * fi_root2_pi < 20.) m2 = 0;

    reci2 = fi_root2_p;
    recj2 = fi_root2_pi;
    for (i = n1; i <= n2; i++) {
      for (j = m1; j <= m2; j++) {
        ti = fi_root2_p + i * 2 * PI;
        tj = fi_root2_pi + j * 2 * PI;
        zij = -proton->Zc - proton->Vz * ti + pion->Zc + pion->Vz * tj;
        if (fabs(zij) < fabs(rec_zij2)) {
          rec_zij2 = zij;
          reci2 = ti;
          recj2 = tj;
        }
      }
    }
    fi_root2_p = reci2;
    fi_root2_pi = recj2;
    d_root2 = rec_zij2;

    dcaPToPi(&fi_root2_p, &fi_root2_pi, proton, pion, &d_root2, &v0_root2[0],
             alfa);
  }

  if (status == 0 || status == 1) {
    fi_p = fi_root1_p;
    fi_pi = fi_root1_pi;
    *x = v0_root1[0];
    *(x + 1) = v0_root1[1];
    *(x + 2) = v0_root1[2];
    d = d_root1;
  } else {
    if (fabs(d_root1) < fabs(d_root2)) {
      fi_p = fi_root1_p;
      fi_pi = fi_root1_pi;
      *x = v0_root1[0];
      *(x + 1) = v0_root1[1];
      *(x + 2) = v0_root1[2];
      d = d_root1;
    } else {
      fi_p = fi_root2_p;
      fi_pi = fi_root2_pi;
      *x = v0_root2[0];
      *(x + 1) = v0_root2[1];
      *(x + 2) = v0_root2[2];
      d = d_root2;
    }
  }

  *p1 = proton->pt * cos(fi_p + proton->h * PI / 2.);
  *p2 = pion->pt * cos(fi_pi + pion->h * PI / 2.);

  *(p1 + 1) = proton->pt * sin(fi_p + proton->h * PI / 2.);
  *(p2 + 1) = pion->pt * sin(fi_pi + pion->h * PI / 2.);
  *(p1 + 2) = proton->Pz;
  *(p2 + 2) = pion->Pz;

  // Finally set the StThreeVectors
  xv0.setX(static_cast<float>(x[0]));
  xv0.setY(static_cast<float>(x[1]));
  xv0.setZ(static_cast<float>(x[2]));

  op1.setX(static_cast<float>(p1[0]));
  op1.setY(static_cast<float>(p1[1]));
  op1.setZ(static_cast<float>(p1[2]));

  op2.setX(static_cast<float>(p2[0]));
  op2.setY(static_cast<float>(p2[1]));
  op2.setZ(static_cast<float>(p2[2]));

  return d;
}

void dcaPToPi(double* fi_p, double* fi_pi, StTrackHelix* proton,
              StTrackHelix* pion, double* d_root, double* v0, double alfa) {
  double R2_TPC = 30000.;  // hard-coded cuts
  double Z_TPC = 180.;
  double t1, t2;

  double Co_a, Co_b, Co_c, Co_e, Co_f, dfi_p, dfi_pi;
  double new_d, x1, y1, z1, x2, y2, z2;
  double old_d, oldx1, oldy1, oldz1, oldx2, oldy2, oldz2;

  Co_a = proton->r * proton->r + proton->Vz * proton->Vz;
  Co_b = -proton->h * pion->h * proton->r * pion->r * cos(alfa) -
         proton->Vz * pion->Vz;  // XZHU: bug fix to consider also the case of
                                 // same charge helices
  Co_c = pion->r * pion->r + pion->Vz * pion->Vz;
  double Tp = proton->r * proton->h * sin(*fi_p * proton->h + proton->theta);
  double Kp = proton->r * proton->h * cos(*fi_p * proton->h + proton->theta);
  double Tn = pion->r * pion->h * sin(*fi_pi * pion->h + pion->theta);
  double Kn = pion->r * pion->h * cos(*fi_pi * pion->h + pion->theta);
  x1 = proton->Xc + proton->r * cos(*fi_p * proton->h + proton->theta);
  y1 = proton->Yc + proton->r * sin(*fi_p * proton->h + proton->theta);
  x2 = pion->Xc + pion->r * cos(*fi_pi * pion->h + pion->theta);
  y2 = pion->Yc + pion->r * sin(*fi_pi * pion->h + pion->theta);
  double xij = x1 - x2;
  double yij = y1 - y2;
  Co_e = (*d_root) * proton->Vz + xij * Tp - yij * Kp;
  Co_f = -(*d_root) * pion->Vz - xij * Tn + yij * Kn;
  dfi_p = 0.;
  dfi_pi = 0.;
  if (fabs(Co_a * Co_c - Co_b * Co_b) > 0.) {
    dfi_p = (Co_e * Co_c - Co_b * Co_f) / (Co_a * Co_c - Co_b * Co_b);

    dfi_pi = (Co_a * Co_f - Co_b * Co_e) / (Co_a * Co_c - Co_b * Co_b);
  }

  // do some minimization for the case of connecting helices.
  double d_p_pi = sqrt((pion->Xc - proton->Xc) * (pion->Xc - proton->Xc) +
                       (pion->Yc - proton->Yc) * (pion->Yc - proton->Yc));
  double r_p = proton->r;  // radium of proton helix
  double r_pi = pion->r;   // radium of pion helix

  if (kMinimize && (d_p_pi >= (r_p + r_pi) - kShiftConnect ||
                    d_p_pi <= fabs(r_p - r_pi) + kShiftContain)) {
    double dfi1, dfi2;

    double dfi1min, dfi2min;
    double dmin = 9999., dminsave;

    double dfi1low = 0;
    double dfi1high = dfi_p / fabs(dfi_p) * fmax(fabs(dfi_p), 0.01);
    double ddfi = (dfi1high - dfi1low) / 10;
    double minddfi = fabs(ddfi) / 1000;

    int nshiftl = 0;
    int nshiftr = 0;
    int ndig = 0;

    while (fabs(ddfi) > minddfi) {
      dminsave = dmin;

      dmin = 9999.;
      double dfi1last;
      for (dfi1 = dfi1low; (dfi1 - dfi1low) / ddfi <= 10.05; dfi1 += ddfi) {
        double dmin12;

        double tt1 = *fi_p + dfi1;
        double tx1, ty1, tz1;
        tx1 = proton->Xc + proton->r * cos(tt1 * proton->h + proton->theta);
        ty1 = proton->Yc + proton->r * sin(tt1 * proton->h + proton->theta);
        tz1 = proton->Zc + proton->Vz * tt1;

        // newton method to find the dca point in pion helix.
        dfi2 = 0;
        double dfi2save;
        int maxinteration = 50;
        int iter = 0;
        double M, M0;
        double tt2, phi2, cosphi2, sinphi2, tx2, ty2, tz2;
        do {
          dfi2save = dfi2;
          if (iter == 0) {
            tt2 = *fi_pi + dfi2;
            phi2 = tt2 * pion->h + pion->theta;
            cosphi2 = cos(phi2);
            sinphi2 = sin(phi2);
            tx2 = pion->Xc + pion->r * cosphi2;
            ty2 = pion->Yc + pion->r * sinphi2;
            tz2 = pion->Zc + pion->Vz * tt2;
            M = sqrt((tx1 - tx2) * (tx1 - tx2) + (ty1 - ty2) * (ty1 - ty2) +
                     (tz1 - tz2) * (tz1 - tz2));
          }
          double dtx2, dty2, dtz2;
          dtx2 = -pion->r * pion->h * sinphi2;
          dty2 = pion->r * pion->h * cosphi2;
          dtz2 = pion->Vz;
          double ddtx2, ddty2, ddtz2;
          ddtx2 = -pion->r * cosphi2;
          ddty2 = -pion->r * sinphi2;
          ddtz2 = 0;

          double fn =
              (tx1 - tx2) * dtx2 + (ty1 - ty2) * dty2 + (tz1 - tz2) * dtz2;
          double dfn = -dtx2 * dtx2 + (tx1 - tx2) * ddtx2 - dty2 * dty2 +
                       (ty1 - ty2) * ddty2 - dtz2 * dtz2 + (tz1 - tz2) * ddtz2;
          dfi2 = dfi2save - fn / dfn;
          M0 = M;
          tt2 = *fi_pi + dfi2;
          phi2 = tt2 * pion->h + pion->theta;
          cosphi2 = cos(phi2);
          sinphi2 = sin(phi2);
          tx2 = pion->Xc + pion->r * cosphi2;
          ty2 = pion->Yc + pion->r * sinphi2;
          tz2 = pion->Zc + pion->Vz * tt2;
          M = sqrt((tx1 - tx2) * (tx1 - tx2) + (ty1 - ty2) * (ty1 - ty2) +
                   (tz1 - tz2) * (tz1 - tz2));
          iter++;
        } while (fabs(M - M0) > 0.001 * micrometer && iter < maxinteration);
        dfi2min = dfi2;
        dmin12 = M;
        // cout<<"iter: "<<iter<<endl;

        if (dmin12 < dmin) {
          dfi1min = dfi1;
          dmin = dmin12;
        }

        dfi1last = dfi1;
      }

      // update ddfi,dfi1low,dfi1high
      if (fabs(static_cast<double>(dfi1min - dfi1low)) < 0.0000001) {
        // cout<<"shift left"<<endl;
        nshiftl++;
        dfi1low = dfi1min - 8 * ddfi;
        dfi1high = dfi1high - 8 * ddfi;
      } else if (fabs(static_cast<double>(dfi1min - dfi1last)) < 0.0000001) {
        // cout<<"shift right"<<endl;
        nshiftr++;
        dfi1low = dfi1low + 8 * ddfi;
        dfi1high = dfi1min + 8 * ddfi;
      } else {
        // cout<<"zooming"<<endl;
        ndig++;
        dfi1low = dfi1min - ddfi;
        dfi1high = dfi1min + ddfi;
      }
      ddfi = (dfi1high - dfi1low) / 10;
      if (nshiftl > 0 && nshiftr > 0) {
        // cout<<"We dropped a candidate"<<endl;
        // cout<<"warn"<<nshiftl<<" "<<nshiftr<<" "<<ndig<<" "<<status<<"
        // "<<10*ddfi<<" "<<dmin<<endl;
        break;
      }
    }
    // if(nshiftl>5||nshiftr>5)cout<<nshiftl<<" "<<nshiftr<<" "<<ndig<<"
    // "<<status<<endl;

    // the finer phase
    dfi_p = dfi1min;
    dfi_pi = dfi2min;
  }

  t1 = *fi_p + dfi_p;
  t2 = *fi_pi + dfi_pi;

  x1 = proton->Xc + proton->r * cos(t1 * proton->h + proton->theta);
  y1 = proton->Yc + proton->r * sin(t1 * proton->h + proton->theta);
  z1 = proton->Zc + proton->Vz * t1;

  x2 = pion->Xc + pion->r * cos(t2 * pion->h + pion->theta);
  y2 = pion->Yc + pion->r * sin(t2 * pion->h + pion->theta);
  z2 = pion->Zc + pion->Vz * t2;
  new_d = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) +
               (z1 - z2) * (z1 - z2));

  if (fabs(z1) > Z_TPC || fabs(z2) > Z_TPC) new_d = 500.;
  if ((x1 * x1 + y1 * y1) > R2_TPC || (x2 * x2 + y2 * y2) > R2_TPC)
    new_d = 500.;

  // added by xianglei, for the case of one (or no) crossing point
  oldx1 =
      proton->Xc + proton->r * cos((t1 - dfi_p) * proton->h + proton->theta);
  oldy1 =
      proton->Yc + proton->r * sin((t1 - dfi_p) * proton->h + proton->theta);
  oldz1 = proton->Zc + proton->Vz * (t1 - dfi_p);

  oldx2 = pion->Xc + pion->r * cos((t2 - dfi_pi) * pion->h + pion->theta);
  oldy2 = pion->Yc + pion->r * sin((t2 - dfi_pi) * pion->h + pion->theta);
  oldz2 = pion->Zc + pion->Vz * (t2 - dfi_pi);

  old_d = sqrt((oldx1 - oldx2) * (oldx1 - oldx2) +
               (oldy1 - oldy2) * (oldy1 - oldy2) +
               (oldz1 - oldz2) * (oldz1 - oldz2));

  if (new_d < old_d) {
    *fi_p = t1 * proton->h + proton->theta;
    *fi_pi = t2 * pion->h + pion->theta;

    *d_root = new_d;

    *v0 = (x2 + x1) / 2;
    *(v0 + 1) = (y2 + y1) / 2;
    *(v0 + 2) = (z2 + z1) / 2;
  } else {
    *fi_p = (t1 - dfi_p) * proton->h + proton->theta;
    *fi_pi = (t2 - dfi_pi) * pion->h + pion->theta;

    *d_root = old_d;

    *v0 = (oldx2 + oldx1) / 2;
    *(v0 + 1) = (oldy2 + oldy1) / 2;
    *(v0 + 2) = (oldz2 + oldz1) / 2;
  }
  return;
}

// the code below this wall is OBSOLETE, use it ONLY IF you have good reason.
////////////////////////////////////////////////////////////////////////////

double getDcaToPV(const StPhysicalHelixD& helix, const StThreeVectorF& pv) {
  // Hui LONG's function to calculate the Dca to PV of a global track.
  // should use inner helix of tracks. good for dca to PV.

  // first prepare the parameter of the helix.
  int hrot = helix.h();
  double psic = helix.phase();
  StThreeVectorF origin = helix.origin();
  double Xc = helix.xcenter() - static_cast<double>(pv.x());
  double Yc = helix.ycenter() - static_cast<double>(pv.y());
  double Zc = static_cast<double>(origin.z() - pv.z());
  double r = 1. / helix.curvature();
  double Vz = r * tan(helix.dipAngle());
  double PI = TMath::Pi();

  // then hand it over to Hui Long's func
  double f;
  double t1, t;
  double t2, a, b, c;
  if (fabs(Vz / r) < 0.2) {
    f = sqrt(Xc * Xc + Yc * Yc);
    t1 = -Xc * r / f;
    t2 = -Yc * r / f;
    a = r * cos(psic);
    b = r * sin(psic);
    c = acos((a * t1 + b * t2) / r / r);
    t1 = sqrt((f - r) * (f - r) + (Zc - c * Vz) * (Zc - c * Vz));
    c -= 2 * PI;
    t2 = sqrt((f - r) * (f - r) + (Zc + c * Vz) * (Zc + c * Vz));
    if (t1 < t2)
      return t1;
    else
      return t2;
  }

  // cout<<"test here: "<<endl;
  a = sqrt(Xc * Xc + Yc * Yc) * r / Vz / Vz;
  c = hrot * (psic - atan2(Yc, Xc));
  b = -Zc / Vz + c;
  int n1, n2;

  if (a > 0.) {
    if (b - a < 0.)
      n1 = static_cast<int>((b - a - PI / 2.) / PI);
    else
      n1 = static_cast<int>((b - a + PI / 2.) / PI);
    if (b + a < 0.)
      n2 = static_cast<int>((b + a - PI / 2.) / PI);
    else
      n2 = static_cast<int>((b + a + PI / 2.) / PI);
  } else {
    if (b + a < 0.)
      n1 = static_cast<int>((b + a - PI / 2.) / PI);
    else
      n1 = static_cast<int>((b + a + PI / 2.) / PI);
    if (b - a < 0.)
      n2 = static_cast<int>((b - a - PI / 2.) / PI);
    else
      n2 = static_cast<int>((b - a + PI / 2.) / PI);
  }
  double bound1, bound2, fb1, fb2, gb1, gb2;
  double save = 1.e+10;

  for (int n = n1; n <= n2; n++) {
    int i = 0;
    bound1 = n * PI - PI / 2.0;
    bound2 = n * PI + PI / 2.0;
    fb1 = a * sin(bound1);
    fb2 = a * sin(bound2);
    gb1 = bound1 - b;
    gb2 = bound2 - b;
    if ((gb1 > fb1 && gb2 > fb2) || (gb1 < fb1 && gb2 < fb2)) continue;
    if (gb1 > fb1 && gb2 < fb2) {
      t1 = bound1;
      t2 = bound2;
    } else {
      t1 = bound2;
      t2 = bound1;
    }
    while (fabs(t2 - t1) > 0.000001) {
      i++;
      t = (t1 + t2) / 2.0;
      if (a * sin(t) + b - t > 0.)
        t2 = t;
      else
        t1 = t;
      if (i > 1000) break;
    };
    f = t1 - c;
    a = sqrt((Xc + r * cos(psic + f * hrot)) * (Xc + r * cos(psic + f * hrot)) +
             (Yc + r * sin(psic + f * hrot)) * (Yc + r * sin(psic + f * hrot)) +
             (Zc + f * Vz) * (Zc + f * Vz));
    if (a < save) save = a;
  }
  return save;
}

double closestDistance(const StThreeVectorF& xv0, const StThreeVectorF& pv0,
                       const StPhysicalHelixD& helix, double magnet,
                       StThreeVectorF& xxi, StThreeVectorF& op) {
  double pxK, pyK, pzK, x, y, z, r, pL, pxL, pyL, pzL, pLK, pK, sign, xc, yc;
  // double P1=0, P2=0, P3=0;
  double a1, b1, c1;
  double tanK = -1;
  double x0, y0, z0;

  // gx gy gz are the decay vertex position of V0
  a1 = static_cast<double>(xv0.x());
  b1 = static_cast<double>(xv0.y());
  c1 = static_cast<double>(xv0.z());
  pxL = static_cast<double>(pv0.x());  // get momentum of V0
  pyL = static_cast<double>(pv0.y());
  pzL = static_cast<double>(pv0.z());

  StThreeVectorF origin = helix.origin();
  x = static_cast<double>(origin.x());
  y = static_cast<double>(origin.y());
  z = static_cast<double>(origin.z());
  StThreeVectorF pbach = helix.momentum(magnet * kilogauss);
  pK = static_cast<double>(pbach.mag());
  pxK = static_cast<double>(pbach.x());
  pyK = static_cast<double>(pbach.y());
  pzK = static_cast<double>(pbach.z());
  // double Xc = helix.xcenter();
  // double Yc = helix.ycenter();
  // double Zc = origin.z();
  r = 1. / helix.curvature();
  // double Vz = r*tan(helix.dipAngle());
  // double pt = (helix.momentum(magnet*kilogauss)).perp();
  // double Pz = (helix.momentum(magnet*kilogauss)).z();

  int hrot = helix.h();

  sign = -hrot;  // get the helicity of the helix, assigned to be
                 // field_reverse*charge
                 // get the x of the center of the pion
  // track //XZHU: why not gtrack->Xc ?
  xc = (x + r * sign * pyK / static_cast<double>(pbach.perp()));

  // get the Y of the center of the pion track
  yc = (y - r * sign * pxK / static_cast<double>(pbach.perp()));
  pL = sqrt(pxL * pxL + pyL * pyL + pzL * pzL);
  pLK = pxL * pxK + pyL * pyK + pzL * pzK;
  double xroot =
      pyL * (pyL * (a1 - xc) - (b1 - yc) * pxL) / (pxL * pxL + pyL * pyL);
  double yroot =
      pxL * (pxL * (b1 - yc) - (a1 - xc) * pyL) / (pxL * pxL + pyL * pyL);
  double dmax =
      fabs(pxL * (b1 - yc) - (a1 - xc) * pyL) / sqrt(pxL * pxL + pyL * pyL);
  // xroot and yroot are the coordinates of the crossing points or the closest
  // point between the pion circle and lambda line ( in 2D ), it is still on the
  // lambda line dmax is the distance from the center of the pion circle to the
  // Lambda
  if (dmax >= (r + 0.5)) return 999.0;
  double newx, newy, newz, sign0 = -1, dsave = 999.0;
  if (dmax >
      r) {  // no 2d crossing points but not too far either, then the newx and
            // newy are the coordinates of the 2D closest point to the Lambda
            // line but still on the circle. //XZHU: should be rare, since it
            // needs the momenta of lambda and K is in line with each other.
    newx = r / dmax * xroot + xc;
    newy = r / dmax * yroot + yc;
    if (((newx - x) * pxK + (newy - y) * pyK) > 0) sign0 = 1.0;
    newz =
        z + sign0 * r * 2 *
                asin(sqrt((newx - x) * (newx - x) + (newy - y) * (newy - y)) /
                     2.0 / r) *
                pzK / sqrt(pxK * pxK + pyK * pyK);
    double newzL = pzL / pyL * (newy - b1) +
                   c1;  // from new x, new y , calculate the newz on the helix
    // dsave=sqrt((newzL-newz)*(newzL-newz)+dmax*dmax); //XZHU: BUG! change dmax
    // -> dmax-r
    dsave = sqrt(
        (newzL - newz) * (newzL - newz) +
        (dmax - r) *
            (dmax - r));  // XZHU: fixed! in this case, we can not use the
                          // straight line approximation below. because these
                          // two lines always cross each other at somewhere far
                          // far away. therefore, you will have unphysical '0'
                          // dca, some of them will add up to the background.
    // tanK=-1;
    tanK = 0;
    // cout<<"I'm here"<<" "<<sign0<<endl;
    x0 = newx;
    y0 = newy;
    z0 = (newz + newzL) / 2.;
    pxK = -(yc - newy) / r * sign * sqrt(pxK * pxK + pyK * pyK);
    pyK = (xc - newx) / r * sign * sqrt(pxK * pxK + pyK * pyK);
  } else {  //  two crssoing points, fa, fb, fc  just intermidiate numbers in
            //  calculation without any geometric meaning. newx, new y has the
            //  same meaning in case 1
    double fa = 1 + pxL * pxL / pyL / pyL;
    double fb = 2 * pxL * ((a1 - xc) - (b1 - yc) * pxL / pyL) / pyL;
    double fc = ((a1 - xc) - (b1 - yc) * pxL / pyL) *
                    ((a1 - xc) - (b1 - yc) * pxL / pyL) -
                r * r;
    double newy1 = (-fb - sqrt(fb * fb - 4 * fa * fc)) * 0.5 / fa;
    double newx1 = pxL * (newy1 - (b1 - yc)) / pyL + (a1 - xc) + xc;
    newy1 = newy1 + yc;
    double newy2 = (-fb + sqrt(fb * fb - 4 * fa * fc)) * 0.5 / fa;
    double newx2 = pxL * (newy2 - (b1 - yc)) / pyL + (a1 - xc) + xc;
    newy2 = newy2 + yc;
    double ds1 = sqrt((newx1 - x) * (newx1 - x) + (newy1 - y) * (newy1 - y));
    double da1 = ds1 / r;
    sign0 = -1;
    dsave = 999;
    if ((newx1 - x) * pxK + (newy1 - y) * pyK > 0) sign0 = 1.0;
    double newz1 =
        z + sign0 * 2.0 * r * asin(da1 / 2.0) * pzK /
                sqrt(pxK * pxK +
                     pyK * pyK);  // XZHU: only work in one period. for
                                  // sign0=1, it is weird to check the
                                  // orign point in the forward direction
    double newzL1 = pzL / pyL * (newy1 - b1) + c1;
    double ds2 = sqrt((newx2 - x) * (newx2 - x) + (newy2 - y) * (newy2 - y));
    double da2 = ds2 / r;
    sign0 = -1;
    if ((newx2 - x) * pxK + (newy2 - y) * pyK > 0) sign0 = 1.0;
    double newz2 = z + sign0 * 2.0 * r * asin(da2 / 2.0) * pzK /
                           sqrt(pxK * pxK + pyK * pyK);
    double newzL2 = pzL / pyL * (newy2 - b1) + c1;
    if (fabs(newzL2 - newz2) > fabs(newzL1 - newz1)) {
      newx = newx1;
      newy = newy1;
      newz = newz1;
      dsave = (newzL1 - newz1);
      // pick up a better one from the two solutions
    } else {
      newx = newx2;
      newy = newy2;
      newz = newz2;
      dsave = (newzL2 - newz2);
    }

    double newpx = -(yc - newy) / r * sign * sqrt(pxK * pxK + pyK * pyK);
    double newpy = (xc - newx) / r * sign * sqrt(pxK * pxK + pyK * pyK);
    x = newx;
    y = newy;
    z = newz;
    pxK = newpx;
    pyK = newpy;  // // update the momentum of the pion track at the new postion

    // JUne 19, 2003   here calculate the distance between two straigt line, (
    // helix pion is approximated as straight line at the new position
    // calculated above. lambda is always a streight line.
    double aK = pxK / pK;
    double bK = pyK / pK;
    double cK = pzK / pK;
    double aL = pxL / pL;
    double bL = pyL / pL;
    double cL = pzL / pL;
    double eA = (aK * aK + bK * bK + cK * cK);
    double eB = -(aK * aL + bK * bL + cK * cL);
    double eC = cK * dsave;
    double eD = -eB;
    double eE = -(aL * aL + bL * bL + cL * cL);
    double eF = cL * dsave;
    double t1 = (eC * eE - eF * eB) / (eA * eE - eB * eD);
    double t2 = (eA * eF - eC * eD) / (eA * eE - eB * eD);
    double min_d =
        sqrt((aK * t1 - aL * t2) * (aK * t1 - aL * t2) +
             (bK * t1 - bL * t2) * (bK * t1 - bL * t2) +
             (cK * t1 - cL * t2 - dsave) * (cK * t1 - cL * t2 - dsave));
    // if(tanK==0)
    // cout<<min_d<<" "<<fabs(dsave)<<" "<<(cK*t1-cL*t2-dsave)<<"
    // "<<(bK*t1-bL*t2)<<" "<<(aK*t1-aL*t2)<<" "<<t1<<" "<<t2<<" "<<aK<<"
    // "<<bK<<" "<<cK<<" "<<aL<<" "<<bL<<" "<<cL<<" "<<dmax<<" "<<r<<"
    // "<<sqrt((newx+aK*t1-xc)*(newx+aK*t1-xc)+(newy+bK*t1-yc)*(newy+bK*t1-yc))<<endl;
    if (min_d < fabs(dsave)) {
      x0 = (aK * t1 + aL * t2) / 2.0 + newx;
      y0 = (bK * t1 + bL * t2) / 2.0 + newy;
      z0 = (cK * t1 + cL * t2 + dsave) / 2.0 + newz;
      dsave = min_d;  // update the x0,y0,z0 as the dca point between Lambda and
                      // pion

      // XZHU: update the Kaon px, py at dca point, should not be necessary.
      // just for consistent check!
      double newkx = aK * t1 + newx;
      double newky = bK * t1 + newy;
      double newr =
          sqrt((newkx - xc) * (newkx - xc) + (newky - yc) * (newky - yc));
      double ptK = sqrt(pxK * pxK + pyK * pyK);
      pxK = -(yc - newky) / newr * sign * ptK;
      pyK = (xc - newkx) / newr * sign * ptK;
    } else {
      x0 = newx;
      y0 = newy;
      z0 = dsave / 2.0 + newz;
    }
  }
  //  P1=pxL+pxK;
  //  P2=pyL+pyK;
  //  P3=pzL+pzK;// updated the hyperon P

  // double PT = sqrt(pow( P1,2) + pow( P2,2));

  // if(fabs(dsave)>sqrt(cutd2)) return 0;  //XZHU: dca of v0 to bachelor CuT
  xxi.setX(static_cast<float>(x0));
  xxi.setY(static_cast<float>(y0));
  xxi.setZ(static_cast<float>(z0));
  op.setX(static_cast<float>(pxK));
  op.setY(static_cast<float>(pyK));
  op.setZ(static_cast<float>(pzK));

  return fabs(dsave);
}
