#ifndef StDcaService_def
#define StDcaService_def

#include "StPhysicalHelixD.hh"
#include "StThreeVectorF.hh"

// this struct is for Hui LONG's dcaPToPi function ONLY.
// it is just the 'Helix' in his code.
// some useless items are commented out.
///////////////////////////////////////////////////////////////////////
/// was added by xiaohai (begin)
///////////////////////////////////////////////////////////////////////
/// The ROOT description: if a particle with charge q passes through a
/// point (x, y, z) with momentum (px, py, pz) with magnetic field B
/// along (nx, ny, nz), this can be constructed like:
/// THelix p(x, y, z, px, py, pz, q*B, nx, ny, nz);
/// (nx, ny, nz) defaults to (0, 0, 1)
///
/// A helix in its own frame can be defined with a pivotal point (x0, y0, z0),
/// the velocity at that point (vx0, vy0, vz0), and an angular frequency w.
/// Combining vx0 and vy0 to a transverse velocity vt0 one can parametrize the
/// helix as:
/// x(t) = x0 - vt0 / w*sin(- w*t + phi0)
/// y(t) = y0 + vt0 / w*cos(- w*t + phi0)
/// z(t) = z0 + vz0 * t
///////////////////////////////////////////////////////////////////////
/// was added by xiaohai (end)
///////////////////////////////////////////////////////////////////////
///
/// \brief The StTrackHelix struct
///        径迹的螺旋结构
struct StTrackHelix {
  // double  pid;
  // short  id;
  double Xc;  ///< xcenter in xy-plane measured from ring center
  double Yc;  ///< ycenter in xy-plane measured from ring center
  double Zc;
  double pt;
  // double X;
  // double Y;
  // double Z;
  // double Px;
  // double Py;
  double Pz;
  double r;
  double theta;
  double Vz;
  int h;
  // int    h2;
  int Flag;  ///< what is flag ?
  // int    q;
  // int    nhits; ///< the number of hit.
  // double  dca;
  // double  dedx;
};

extern bool kStHelixDca;
extern bool kMinimize;
extern double kShiftConnect;
extern double kShiftContain;

/// \brief closestDistance wrapper to Hui LONG's two helix dca code
/// \param helix1  first particle
/// \param helix2 second particle
/// \param magnet
/// \param pv
/// \param xv0
/// \param op1
/// \param op2
/// \return
///
double closestDistance(const StPhysicalHelixD& helix1,
                       const StPhysicalHelixD& helix2, double magnet,
                       const StThreeVectorF& pv, StThreeVectorF& xv0,
                       StThreeVectorF& op1, StThreeVectorF& op2);
/// P: is proton
/// Pi: is piplus or piminus
void dcaPToPi(double* fi_p, double* fi_pi, StTrackHelix* proton,
              StTrackHelix* pion, double* d_root, double* v0, double alfa);

/// get the distance from Dca to primary vertex.(OBSOLETE)
double getDcaToPV(const StPhysicalHelixD& helix, const StThreeVectorF& pv);

/// wrapper to Hui LONG's helix to straight line Dca code(for Xi or Omega
/// reconstruction) (OBSOLETE)
double closestDistance(const StThreeVectorF& xv0, const StThreeVectorF& pv0,
                       const StPhysicalHelixD& helix, double magnet,
                       StThreeVectorF& xxi, StThreeVectorF& op);

#endif
