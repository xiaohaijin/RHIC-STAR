// SpaceShower.h is a part of the PYTHIA event generator.
// Copyright (C) 2012 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the spacelike initial-state showers.
// SpaceDipoleEnd: radiating dipole end in ISR.
// SpaceShower: handles the showering description.

#ifndef Pythia8_SpaceShower_H
#define Pythia8_SpaceShower_H

#include "Basics.h"
#include "BeamParticle.h"
#include "Event.h"
#include "Info.h"
#include "ParticleData.h"
#include "PartonSystems.h"
#include "PythiaStdlib.h"
#include "Settings.h"
#include "StandardModel.h"
#include "UserHooks.h"

namespace Pythia8 {
 
//==========================================================================

// Data on radiating dipole ends, only used inside SpaceShower.

class SpaceDipoleEnd {
  
public:

  // Constructor.
  SpaceDipoleEnd( int systemIn = 0, int sideIn = 0, int iRadiatorIn = 0, 
    int iRecoilerIn = 0, double pTmaxIn = 0., int colTypeIn = 0, 
    int chgTypeIn = 0,  int MEtypeIn = 0, bool normalRecoilIn = true) : 
    system(systemIn), side(sideIn), iRadiator(iRadiatorIn), 
    iRecoiler(iRecoilerIn), pTmax(pTmaxIn), colType(colTypeIn), 
    chgType(chgTypeIn), MEtype(MEtypeIn), normalRecoil(normalRecoilIn), 
    nBranch(0), pT2Old(0.), zOld(0.5) { }
 
  // Store values for trial emission.
  void store( int idDaughterIn, int idMotherIn, int idSisterIn,   
    double x1In, double x2In, double m2DipIn, double pT2In, double zIn,
    double xMoIn, double Q2In, double mSisterIn, double m2SisterIn, 
    double pT2corrIn) {idDaughter = idDaughterIn; idMother = idMotherIn; 
    idSister = idSisterIn; x1 = x1In; x2 = x2In; m2Dip = m2DipIn; 
    pT2 = pT2In; z = zIn; xMo = xMoIn; Q2 = Q2In; mSister = mSisterIn; 
    m2Sister = m2SisterIn; pT2corr = pT2corrIn;}
 
  // Basic properties related to evolution and matrix element corrections.
  int    system, side, iRadiator, iRecoiler;
  double pTmax;
  int    colType, chgType, MEtype;
  bool   normalRecoil;
  
  // Properties specific to current trial emission.
  int    nBranch, idDaughter, idMother, idSister, iFinPol;  
  double x1, x2, m2Dip, pT2, z, xMo, Q2, mSister, m2Sister, pT2corr, 
         pT2Old, zOld, asymPol;

} ;
 
//==========================================================================

// The SpaceShower class does spacelike showers.

class SpaceShower {

public:

  // Constructor.
  SpaceShower() {}

  // Destructor.
  virtual ~SpaceShower() {}

  // Initialize various pointers.
  // (Separated from rest of init since not virtual.)
  void initPtr(Info* infoPtrIn, Settings* settingsPtrIn, 
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
    PartonSystems* partonSystemsPtrIn, UserHooks* userHooksPtrIn)  {
    infoPtr = infoPtrIn; settingsPtr = settingsPtrIn; 
    particleDataPtr = particleDataPtrIn; rndmPtr = rndmPtrIn; 
    partonSystemsPtr = partonSystemsPtrIn; userHooksPtr = userHooksPtrIn;}

  // Initialize generation. Possibility to force re-initialization by hand.
  virtual void init(BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn);

  // New beams possible for handling of hard diffraction. (Not virtual.)
  void reassignBeamPtrs( BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn) 
    {beamAPtr = beamAPtrIn; beamBPtr = beamBPtrIn;}

  // Find whether to limit maximum scale of emissions, and whether to dampen.
  virtual bool limitPTmax( Event& event, double Q2Fac = 0., 
    double Q2Ren = 0.);

  // Potential enhancement factor of pTmax scale for hardest emission.
  virtual double enhancePTmax() const {return pTmaxFudge;}

  // Prepare system for evolution; identify ME.
  virtual void prepare( int iSys, Event& event, bool limitPTmaxIn = true);

  // Update dipole list after each FSR emission. Currently superfluous.
  // Usage: update( iSys, event).  
  virtual void update( int , Event& ) {}

  // Select next pT in downwards evolution.
  virtual double pTnext( Event& event, double pTbegAll, double pTendAll,
    int nRadIn = -1);

  // ME corrections and kinematics that may give failure.
  virtual bool branch( Event& event); 

  // Tell which system was the last processed one.
  int system() const {return iSysSel;} 

  // Flag for failure in branch(...) that will force a retry of parton level.
  bool doRestart() const {return rescatterFail;}

  // Print dipole list; for debug mainly.
  virtual void list(ostream& os = cout) const;

protected:

  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to the settings database.
  Settings*      settingsPtr;

  // Pointer to the particle data table.
  ParticleData*  particleDataPtr;

  // Pointer to the random number generator.
  Rndm*          rndmPtr;

  // Pointers to the two incoming beams.
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;

  // Pointer to information on subcollision parton locations.
  PartonSystems* partonSystemsPtr;

  // Pointer to userHooks object for user interaction with program.
  UserHooks*     userHooksPtr;

  // Store properties to be returned by methods.
  bool   rescatterFail;
  int    iSysSel;
  double pTmaxFudge;

private: 

  // Constants: could only be changed in the code itself.
  static const bool   DEBUG;
  static const int    MAXLOOPTINYPDF;
  static const double CTHRESHOLD, BTHRESHOLD, EVALPDFSTEP, TINYPDF, 
         TINYKERNELPDF, TINYPT2, HEAVYPT2EVOL, HEAVYXEVOL, EXTRASPACEQ, 
         LAMBDA3MARGIN, LEPTONXMIN, LEPTONXMAX, LEPTONPT2MIN, LEPTONFUDGE;

  // Initialization data, normally only set once.
  bool   doQCDshower, doQEDshowerByQ, doQEDshowerByL, useSamePTasMPI,
         doMEcorrections, doMEafterFirst, doPhiPolAsym, doPhiIntAsym, 
         doRapidityOrder, canVetoEmission;
  int    pTmaxMatch, pTdampMatch, alphaSorder, alphaEMorder, nQuarkIn, 
         enhanceScreening;
  double pTdampFudge, mc, mb, m2c, m2b, alphaSvalue, alphaS2pi, 
         Lambda3flav, Lambda4flav, Lambda5flav, Lambda3flav2, Lambda4flav2, 
         Lambda5flav2, pT0Ref, ecmRef, ecmPow, pTmin, sCM, eCM, pT0, 
         pTminChgQ, pTminChgL, pT20, pT2min, pT2minChgQ, pT2minChgL, 
         pTmaxFudgeMPI, strengthIntAsym; 

  // alphaStrong and alphaEM calculations.
  AlphaStrong alphaS;
  AlphaEM alphaEM;

  // Some current values.
  bool   sideA, dopTdamp;
  int    iNow, iRec, idDaughter, nRad, idResFirst, idResSecond;
  double xDaughter, x1Now, x2Now, m2Dip, m2Rec, pT2damp, pTbegRef;

  // All dipole ends
  vector<SpaceDipoleEnd> dipEnd;

  // Pointers to the current and hardest (so far) dipole ends.
  int iDipNow, iSysNow;
  SpaceDipoleEnd* dipEndNow; 
  int iDipSel;
  SpaceDipoleEnd* dipEndSel; 
 
  // Evolve a QCD dipole end. 
  void pT2nextQCD( double pT2begDip, double pT2endDip);

  // Evolve a QCD dipole end near heavy quark threshold region. 
  void pT2nearQCDthreshold( BeamParticle& beam, double m2Massive, 
    double m2Threshold, double xMaxAbs, double zMinAbs, 
    double zMaxMassive);

  // Evolve a QED dipole end. 
  void pT2nextQED( double pT2begDip, double pT2endDip);

  // Find class of ME correction.
  int findMEtype( int iSys, Event& event);

  // Provide maximum of expected ME weight; for preweighting of evolution.
  double calcMEmax( int MEtype, int idMother, int idDaughterIn);

  // Provide actual ME weight for current branching.
  double calcMEcorr(int MEtype, int idMother, int idDaughterIn, double M2, 
    double z, double Q2); 

  // Find coefficient of azimuthal asymmetry from gluon polarization.
  void findAsymPol( Event& event, SpaceDipoleEnd* dip);

};
 
//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SpaceShower_H
