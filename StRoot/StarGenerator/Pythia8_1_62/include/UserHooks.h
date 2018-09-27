// UserHooks.h is a part of the PYTHIA event generator.
// Copyright (C) 2012 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file to allow user access to program at different stages.
// UserHooks: almost empty base class, with user to write the rela code. 
// MyUserHooks: derived class, only intended as an example.

#ifndef Pythia8_UserHooks_H
#define Pythia8_UserHooks_H

#include "Event.h"
#include "PartonSystems.h"
#include "PythiaStdlib.h"
#include "SigmaProcess.h"

namespace Pythia8 {

//==========================================================================

// Forward reference to the PhaseSpace class.
class PhaseSpace;

//==========================================================================

// UserHooks is base class for user access to program execution.

class UserHooks {

public:

  // Destructor.
  virtual ~UserHooks() {selBias = 1.;}

  // Initialize pointers and workEvent. Note: not virtual.
  void initPtr( Info* infoPtrIn, Settings* settingsPtrIn, 
    ParticleData* particleDataPtrIn,  Rndm* rndmPtrIn, 
    BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn, 
    BeamParticle* beamPomAPtrIn, BeamParticle* beamPomBPtrIn, 
    CoupSM* coupSMPtrIn, PartonSystems* partonSystemsPtrIn, 
    SigmaTotal* sigmaTotPtrIn) { infoPtr = infoPtrIn; 
    settingsPtr = settingsPtrIn; particleDataPtr = particleDataPtrIn;
    rndmPtr = rndmPtrIn; beamAPtr = beamAPtrIn; beamBPtr = beamBPtrIn;
    beamPomAPtr = beamPomAPtrIn; beamPomBPtr = beamPomBPtrIn;
    coupSMPtr = coupSMPtrIn; partonSystemsPtr = partonSystemsPtrIn;
    sigmaTotPtr = sigmaTotPtrIn;
    workEvent.init("(work event)", particleDataPtr);}

  // Possibility to modify cross section of process.
  virtual bool canModifySigma() {return false;}

  // Multiplicative factor modifying the cross section of a hard process.
  virtual double multiplySigmaBy(const SigmaProcess* sigmaProcessPtr,
    const PhaseSpace* phaseSpacePtr, bool inEvent);

  // Possibility to bias selection of events, compensated by a weight.
  virtual bool canBiasSelection() {return false;}

  // Multiplicative factor in the phase space selection of a hard process.
  virtual double biasSelectionBy(const SigmaProcess* sigmaProcessPtr,
    const PhaseSpace* phaseSpacePtr, bool inEvent);

  // Event weight to compensate for selection weight above.
  virtual double biasedSelectionWeight() {return 1./selBias;}  

  // Possibility to veto event after process-level selection.
  virtual bool canVetoProcessLevel() {return false;}

  // Decide whether to veto current process or not, based on process record.
  // Usage: doVetoProcessLevel( process).
  virtual bool doVetoProcessLevel(Event& ) {return false;}

  // Possibility to veto MPI + ISR + FSR evolution and kill event, 
  // making decision at a fixed pT scale. Useful for MLM-style matching.
  virtual bool canVetoPT() {return false;}  

  // Transverse-momentum scale for veto test. 
  virtual double scaleVetoPT() {return 0.;} 

  // Decide whether to veto current event or not, based on event record.
  // Usage: doVetoPT( iPos, event), where iPos = 0: no emissions so far;
  // iPos = 1/2/3 joint evolution, latest step was MPI/ISR/FSR;
  // iPos = 4: FSR only afterwards; iPos = 5: FSR in resonance decay. 
  virtual bool doVetoPT( int , const Event& ) {return false;} 

  // Possibility to veto MPI + ISR + FSR evolution and kill event, 
  // making decision after fixed number of ISR or FSR steps.
  virtual bool canVetoStep() {return false;}

  // Up to how many ISR + FSR steps of hardest interaction should be checked.
  virtual int numberVetoStep() {return 1;}

  // Decide whether to veto current event or not, based on event record.
  // Usage: doVetoStep( iPos, nISR, nFSR, event), where iPos as above,
  // nISR and nFSR number of emissions so far for hard interaction only.
  virtual bool doVetoStep( int , int , int , const Event& ) {return false;} 

  // Possibility to veto MPI + ISR + FSR evolution and kill event, 
  // making decision after fixed number of MPI steps.
  virtual bool canVetoMPIStep() {return false;}

  // Up to how many MPI steps should be checked.
  virtual int numberVetoMPIStep() {return 1;}

  // Decide whether to veto current event or not, based on event record.
  // Usage: doVetoMPIStep( nMPI, event), where nMPI is number of MPI's so far.
  virtual bool doVetoMPIStep( int , const Event& ) {return false;} 
   
  // Possibility to veto event after parton-level selection.
  virtual bool canVetoPartonLevel() {return false;}

  // Decide whether to veto current partons or not, based on event record.
  // Usage: doVetoPartonLevel( event).
  virtual bool doVetoPartonLevel( const Event& ) {return false;} 

  // Possibility to set initial scale in TimeShower for resonance decay.
  virtual bool canSetResonanceScale() {return false;}

  // Initial scale for TimeShower evolution. 
  // Usage: scaleResonance( iRes, event), where iRes is location
  // of decaying resonance in the event record. 
  virtual double scaleResonance( int, const Event& ) {return 0.;} 

  // Possibility to veto an emission in the ISR machinery.
  virtual bool canVetoISREmission() {return false;}

  // Decide whether to veto current emission or not, based on event record.
  // Usage: doVetoISREmission( sizeOld, event, iSys) where sizeOld is size
  // of event record before current emission-to-be-scrutinized was added,
  // and iSys is the system of the radiation (according to PartonSystems).
  virtual bool doVetoISREmission( int, const Event&, int ) {return false;} 

  // Possibility to veto an emission in the FSR machinery.
  virtual bool canVetoFSREmission() {return false;}

  // Decide whether to veto current emission or not, based on event record.
  // Usage: doVetoFSREmission( sizeOld, event, iSys, inResonance) where 
  // sizeOld is size of event record before current emission-to-be-scrutinized 
  // was added, iSys is the system of the radiation (according to 
  // PartonSystems), and inResonance is true if the emission takes place in a
  // resonance decay.
  virtual bool doVetoFSREmission( int, const Event&, int, bool = false )
      {return false;} 

  // Possibility to veto an MPI.
  virtual bool canVetoMPIEmission() { return false; }

  // Decide whether to veto an MPI based on event record.
  // Usage: doVetoMPIEmission( sizeOld, event) where sizeOld
  // is size of event record before the current MPI.
  virtual bool doVetoMPIEmission(int, const Event &) { return false; }

protected:

  // Constructor.
  UserHooks() {}

  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to the settings database.
  Settings*      settingsPtr;

  // Pointer to the particle data table.
  ParticleData*  particleDataPtr;

 // Pointer to the random number generator.
  Rndm*          rndmPtr;

  // Pointers to the two incoming beams and to Pomeron beam-inside-beam.
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;
  BeamParticle*  beamPomAPtr;
  BeamParticle*  beamPomBPtr;

  // Pointers to Standard Model couplings.
  CoupSM*        coupSMPtr;

  // Pointer to information on subcollision parton locations.
  PartonSystems* partonSystemsPtr;

  // Pointer to the total/elastic/diffractive cross sections.
  SigmaTotal*    sigmaTotPtr;

  // omitResonanceDecays omits resonance decay chains from process record.
  void omitResonanceDecays(const Event& process); 

  // subEvent extracts currently resolved partons in the hard process.
  void subEvent(const Event& event, bool isHardest = true); 

  // Have one event object around as work area.
  Event workEvent;

  // User-imposed selection bias.
  double selBias;

};

//==========================================================================

// SuppressSmallPT is a derived class for user access to program execution.
// It is a simple example, illustrating how to suppress the cross section
// of 2 -> 2 processes by a factor pT^4 / (pT0^2 + pT^2)^2, with pT0 input,
// and also modify alpha_strong scale similarly.

class SuppressSmallPT : public UserHooks {
 
public:

  // Constructor.
  SuppressSmallPT( double pT0timesMPIIn = 1., int numberAlphaSIn = 0, 
    bool useSameAlphaSasMPIIn = true) {isInit = false; 
    pT0timesMPI = pT0timesMPIIn; numberAlphaS = numberAlphaSIn; 
    useSameAlphaSasMPI = useSameAlphaSasMPIIn;}

  // Possibility to modify cross section of process.
  virtual bool canModifySigma() {return true;}

  // Multiplicative factor modifying the cross section of a hard process.
  // Usage: inEvent is true for event generation, false for initialization.
  virtual double multiplySigmaBy(const SigmaProcess* sigmaProcessPtr,
    const PhaseSpace* phaseSpacePtr, bool );

private:

  // Save input properties and the squared pT0 scale.
  bool   isInit, useSameAlphaSasMPI;
  int    numberAlphaS;
  double pT0timesMPI, pT20;

  // Alpha_strong calculation.
  AlphaStrong alphaS;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_UserHooks_H
