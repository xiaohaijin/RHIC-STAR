// ProcessContainer.cc is a part of the PYTHIA event generator.
// Copyright (C) 2012 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the 
// ProcessContainer and SetupContainers classes.

#include "ProcessContainer.h"

// Internal headers for special processes.
#include "SigmaCompositeness.h"
#include "SigmaEW.h"
#include "SigmaExtraDim.h"
#include "SigmaGeneric.h"
#include "SigmaHiggs.h"
#include "SigmaLeftRightSym.h"
#include "SigmaLeptoquark.h"
#include "SigmaNewGaugeBosons.h"
#include "SigmaOnia.h"
#include "SigmaQCD.h"
#include "SigmaSUSY.h"

namespace Pythia8 {

//==========================================================================

// ProcessContainer class.
// Information allowing the generation of a specific process.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of event tries to check maximization finding reliability.
const int ProcessContainer::N12SAMPLE = 100;

// Ditto, but increased for 2 -> 3 processes.
const int ProcessContainer::N3SAMPLE  = 1000;

//--------------------------------------------------------------------------

// Initialize phase space and counters.
// Argument isFirst distinguishes two hard processes in same event.

bool ProcessContainer::init(bool isFirst, Info* infoPtrIn, 
  Settings& settings, ParticleData* particleDataPtrIn, Rndm* rndmPtrIn, 
  BeamParticle* beamAPtr, BeamParticle* beamBPtr, Couplings* couplingsPtr, 
  SigmaTotal* sigmaTotPtr, ResonanceDecays* resDecaysPtrIn, 
  SusyLesHouches* slhaPtr, UserHooks* userHooksPtr) {

  // Extract info about current process from SigmaProcess object.
  isLHA       = sigmaProcessPtr->isLHA();
  isMinBias   = sigmaProcessPtr->isMinBias();
  isResolved  = sigmaProcessPtr->isResolved();
  isDiffA     = sigmaProcessPtr->isDiffA();
  isDiffB     = sigmaProcessPtr->isDiffB();
  isQCD3body  = sigmaProcessPtr->isQCD3body();
  int nFin    = sigmaProcessPtr->nFinal();
  lhaStrat    = (isLHA) ? lhaUpPtr->strategy() : 0;
  lhaStratAbs = abs(lhaStrat);
  allowNegSig = sigmaProcessPtr->allowNegativeSigma();

  // Flag for maximum violation handling.
  increaseMaximum = settings.flag("PhaseSpace:increaseMaximum");

  // Pick and create phase space generator. Send pointers where required.
  if      (isLHA)       phaseSpacePtr = new PhaseSpaceLHA();
  else if (isMinBias)   phaseSpacePtr = new PhaseSpace2to2minbias();
  else if (!isResolved && !isDiffA  && !isDiffB )
                        phaseSpacePtr = new PhaseSpace2to2elastic();
  else if (!isResolved) phaseSpacePtr = new PhaseSpace2to2diffractive( 
                                        isDiffA, isDiffB);
  else if (nFin == 1)   phaseSpacePtr = new PhaseSpace2to1tauy();
  else if (nFin == 2)   phaseSpacePtr = new PhaseSpace2to2tauyz();
  else if (isQCD3body)  phaseSpacePtr = new PhaseSpace2to3yyycyl();
  else                  phaseSpacePtr = new PhaseSpace2to3tauycyl();

  // Store pointers and perform simple initialization.
  infoPtr         = infoPtrIn;
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;
  resDecaysPtr    = resDecaysPtrIn;
  if (isLHA) {
    sigmaProcessPtr->setLHAPtr(lhaUpPtr);
    phaseSpacePtr->setLHAPtr(lhaUpPtr);
  }
  sigmaProcessPtr->init(infoPtr, &settings, particleDataPtr, rndmPtr, 
    beamAPtr, beamBPtr, couplingsPtr, sigmaTotPtr, slhaPtr);
  phaseSpacePtr->init( isFirst, sigmaProcessPtr, infoPtr, &settings,
    particleDataPtr, rndmPtr, beamAPtr,  beamBPtr, couplingsPtr, sigmaTotPtr, 
    userHooksPtr);

  // Reset cross section statistics.
  nTry      = 0;
  nSel      = 0;
  nAcc      = 0;
  nTryStat  = 0;
  sigmaMx   = 0.;
  sigmaSum  = 0.;
  sigma2Sum = 0.;
  sigmaNeg  = 0.;
  sigmaAvg  = 0.;
  sigmaFin  = 0.;
  deltaFin  = 0.;
  wtAccSum  = 0.;

  // Initialize process and allowed incoming partons.
  sigmaProcessPtr->initProc();
  if (!sigmaProcessPtr->initFlux()) return false;

  // Find maximum of differential cross section * phasespace.
  bool physical       = phaseSpacePtr->setupSampling();
  sigmaMx             = phaseSpacePtr->sigmaMax();
  double sigmaHalfWay = sigmaMx;

  // Separate signed maximum needed for LHA with negative weight.
  sigmaSgn            = phaseSpacePtr->sigmaSumSigned();

  // Check maximum by a few events, and extrapolate a further increase.
  if (physical & !isLHA) {
    int nSample = (nFin < 3) ? N12SAMPLE : N3SAMPLE;
    for (int iSample = 0; iSample < nSample; ++iSample) {
      bool test = false;
      while (!test) test = phaseSpacePtr->trialKin(false); 
      if (iSample == nSample/2) sigmaHalfWay = phaseSpacePtr->sigmaMax();
    }   
    sigmaMx = pow2(phaseSpacePtr->sigmaMax()) / sigmaHalfWay;
    phaseSpacePtr->setSigmaMax(sigmaMx);
  }

  // Done.
  return physical;
}

//--------------------------------------------------------------------------

// Generate a trial event; selected or not.
 
bool ProcessContainer::trialProcess() { 

  // Loop over tries only occurs for Les Houches strategy = +-2.
  for (int iTry = 0;  ; ++iTry) {

    // Generate a trial phase space point, if meaningful.
    if (sigmaMx == 0.) return false;
    infoPtr->setEndOfFile(false);
    bool repeatSame = (iTry > 0);
    bool physical = phaseSpacePtr->trialKin(true, repeatSame);

    // Possibly fail, e.g. if at end of Les Houches file, else cross section.
    if (isLHA && !physical) infoPtr->setEndOfFile(true);
    else ++nTry;
    if (!physical) return false;
    double sigmaNow = phaseSpacePtr->sigmaNow(); 

    // Tell if this event comes with weight from cross section.
    double sigmaWeight = 1.;
    if (!isLHA && !increaseMaximum && sigmaNow > sigmaMx) 
      sigmaWeight = sigmaNow / sigmaMx;
    if ( lhaStrat < 0 && sigmaNow < 0.) sigmaWeight = -1.;
    if ( lhaStratAbs == 4) sigmaWeight = sigmaNow;

    // Also compensating weight from biased phase-space selection. 
    double biasWeight = phaseSpacePtr->biasSelectionWeight();
    weightNow = sigmaWeight * biasWeight;
    infoPtr->setWeight( weightNow, lhaStrat);

    // Check that not negative cross section when not allowed.
    if (!allowNegSig) {
      if (sigmaNow < sigmaNeg) {
        infoPtr->errorMsg("Warning in ProcessContainer::trialProcess: neg"
          "ative cross section set 0", "for " +  sigmaProcessPtr->name() );
        sigmaNeg = sigmaNow;
      }
      if (sigmaNow < 0.) sigmaNow = 0.;
    }

    // Cross section of process may come with a weight. Update sum.
    double sigmaAdd = sigmaNow * biasWeight;
    if (lhaStratAbs == 2 || lhaStratAbs == 3) sigmaAdd = sigmaSgn;
    sigmaSum  += sigmaAdd;
    sigma2Sum += pow2(sigmaAdd);
 
    // Check if maximum violated.
    newSigmaMx = phaseSpacePtr->newSigmaMax();
    if (newSigmaMx) sigmaMx = phaseSpacePtr->sigmaMax();

    // Select or reject trial point.
    bool select = true;
    if (lhaStratAbs < 3) select 
      = (newSigmaMx || rndmPtr->flat() * abs(sigmaMx) < abs(sigmaNow)); 
    if (select) ++nSel;
    if (select || lhaStratAbs != 2) return select;
  }
 
}


//--------------------------------------------------------------------------
  
// Give the hard subprocess.

bool ProcessContainer::constructProcess( Event& process, bool isHardest) { 

  // Construct flavour and colours for selected event.
  if (isResolved && !isMinBias) sigmaProcessPtr->pickInState();
  sigmaProcessPtr->setIdColAcol();

  // Construct kinematics from selected phase space point.
  if (!phaseSpacePtr->finalKin()) return false;
  int nFin = sigmaProcessPtr->nFinal();

  // Basic info on process.
  if (isHardest) infoPtr->setType( name(), code(), nFin, isMinBias, 
    isResolved, isDiffA, isDiffB, isLHA);

  // Let hard process record begin with the event as a whole and
  // the two incoming beam particles.  
  process.append( 90, -11, 0, 0, 0, 0, 0, 0, 
    Vec4(0., 0., 0., infoPtr->eCM()), infoPtr->eCM(), 0. ); 
  process.append( infoPtr->idA(), -12, 0, 0, 0, 0, 0, 0, 
    Vec4(0., 0., infoPtr->pzA(), infoPtr->eA()), infoPtr->mA(), 0. ); 
  process.append( infoPtr->idB(), -12, 0, 0, 0, 0, 0, 0, 
    Vec4(0., 0., infoPtr->pzB(), infoPtr->eB()), infoPtr->mB(), 0. ); 

  // For minbias process no interaction selected so far, so done.
  if (isMinBias) return true;

  // Entries 3 and 4, now to be added, come from 1 and 2.
  process[1].daughter1(3);
  process[2].daughter1(4);
  double scale = 0.;

  // Insert the subprocess partons - resolved processes.
  int idRes = sigmaProcessPtr->idSChannel();
  if (isResolved && !isLHA) {

    // NOAM: Mothers and daughters without/with intermediate state.
    int m_M1 = 3; 
    int m_M2 = 4; 
    int m_D1 = 5; 
    int m_D2 = 4 + nFin;
    if (idRes != 0) { 
      m_M1   = 5;
      m_M2   = 0;
      m_D1   = 5;  
      m_D2   = 0; 
    }

    // Find scale from which to begin MPI/ISR/FSR evolution. 
    scale = sqrt(Q2Fac());
    process.scale( scale );

    // Loop over incoming and outgoing partons.
    int colOffset = process.lastColTag();
    for (int i = 1; i <= 2 + nFin; ++i) { 

      // Read out particle info from SigmaProcess object.
      int id        = sigmaProcessPtr->id(i);
      int status    = (i <= 2) ? -21 : 23;
      int mother1   = (i <= 2) ? i : m_M1 ;
      int mother2   = (i <= 2) ? 0 : m_M2 ;
      int daughter1 = (i <= 2) ? m_D1 : 0;
      int daughter2 = (i <= 2) ? m_D2 : 0;
      int col       = sigmaProcessPtr->col(i);
      if      (col > 0) col += colOffset;
      else if (col < 0) col -= colOffset;
      int acol      = sigmaProcessPtr->acol(i);
      if      (acol > 0) acol += colOffset;
      else if (acol < 0) acol -= colOffset;

      // Append to process record.
      int iNow = process.append( id, status, mother1, mother2, 
        daughter1, daughter2, col, acol, phaseSpacePtr->p(i), 
        phaseSpacePtr->m(i), scale);

      // NOAM: If there is an intermediate state, insert the it in 
      // the process record after the two incoming particles.
      if (i == 2 && idRes != 0) {

        // Sign of intermediate state: go by charge. 
        if (particleDataPtr->hasAnti(idRes)
          && process[3].chargeType() + process[4].chargeType() < 0) 
          idRes *= -1;

        // The colour configuration of the intermediate state has to be 
        // resolved separately.
        col         = 0;
        acol        = 0;
        int m_col1  = sigmaProcessPtr->col(1); 
        int m_acol1 = sigmaProcessPtr->acol(1);
        int m_col2  = sigmaProcessPtr->col(2); 
        int m_acol2 = sigmaProcessPtr->acol(2);
        if (m_col1 == m_acol2 && m_col2 != m_acol1) { 
          col       = m_col2; 
          acol       = m_acol1; 
        } else if (m_col2 == m_acol1 && m_col1 != m_acol2) { 
          col        = m_col1; 
          acol       = m_acol2; 
        }
        if      ( col > 0)  col += colOffset;
	else if ( col < 0)  col -= colOffset;
        if      (acol > 0) acol += colOffset;
	else if (acol < 0) acol -= colOffset;

        // Insert the intermediate state into the event record.
        Vec4 pIntMed = phaseSpacePtr->p(1) + phaseSpacePtr->p(2);
        process.append( idRes, -22, 3, 4,  6, 5 + nFin, col, acol, 
          pIntMed, pIntMed.mCalc(), scale);
      }
      
      // Pick lifetime where relevant, else not.
      if (process[iNow].tau0() > 0.) process[iNow].tau(
        process[iNow].tau0() * rndmPtr->exp() );
    }
  }

  // Insert the outgoing particles - unresolved processes.
  else if (!isLHA) {
    int id3     = sigmaProcessPtr->id(3);
    int status3 = (id3 == process[1].id()) ? 14 : 15;
    process.append( id3, status3, 1, 0, 0, 0, 0, 0, 
      phaseSpacePtr->p(3), phaseSpacePtr->m(3));
    int id4     = sigmaProcessPtr->id(4);
    int status4 = (id4 == process[2].id()) ? 14 : 15;
    process.append( id4, status4, 2, 0, 0, 0, 0, 0, 
      phaseSpacePtr->p(4), phaseSpacePtr->m(4));
  }

  // Insert the outgoing particles - Les Houches Accord processes.
  else {

    // Since LHA partons may be out of order, determine correct one.
    // (Recall that zeroth particle is empty.) 
    vector<int> newPos;
    newPos.reserve(lhaUpPtr->sizePart());
    newPos.push_back(0);
    for (int iNew = 0; iNew < lhaUpPtr->sizePart(); ++iNew) {
      // For iNew == 0 look for the two incoming partons, then for
      // partons having them as mothers, and so on layer by layer.
      for (int i = 1; i < lhaUpPtr->sizePart(); ++i)
        if (lhaUpPtr->mother1(i) == newPos[iNew]) newPos.push_back(i);
      if (int(newPos.size()) <= iNew) break;
    } 

    // Find scale from which to begin MPI/ISR/FSR evolution.
    scale = lhaUpPtr->scale();
    double scalePr = (scale < 0.) ? sqrt(Q2Fac()) : scale;
    process.scale( scalePr);

    // Copy over info from LHA event to process, in proper order.
    for (int i = 1; i < lhaUpPtr->sizePart(); ++i) {
      int iOld = newPos[i];
      int id = lhaUpPtr->id(iOld);

      // Translate from LHA status codes.
      int lhaStatus =  lhaUpPtr->status(iOld);
      int status = -21;
      if (lhaStatus == 2 || lhaStatus == 3) status = -22;
      if (lhaStatus == 1) status = 23;

      // Find where mothers have been moved by reordering.
      int mother1Old = lhaUpPtr->mother1(iOld);   
      int mother2Old = lhaUpPtr->mother2(iOld);   
      int mother1 = 0;
      int mother2 = 0; 
      for (int im = 1; im < i; ++im) {
        if (mother1Old == newPos[im]) mother1 = im + 2; 
        if (mother2Old == newPos[im]) mother2 = im + 2; 
      } 
      if (i <= 2) mother1 = i;

      // Ensure that second mother = 0 except for bona fide carbon copies.
      if (mother1 > 0 && mother2 == mother1) { 
        int sister1 = process[mother1].daughter1();
        int sister2 = process[mother1].daughter2();
        if (sister2 != sister1 && sister2 != 0) mother2 = 0;
      } 

      // Find daughters and where they have been moved by reordering. 
      // (Values shifted two steps to account for inserted beams.)
      int daughter1 = 0;
      int daughter2 = 0;
      for (int im = i + 1; im < lhaUpPtr->sizePart(); ++im) { 
        if (lhaUpPtr->mother1(newPos[im]) == iOld
          || lhaUpPtr->mother2(newPos[im]) == iOld) {
          if (daughter1 == 0 || im + 2 < daughter1) daughter1 = im + 2;
          if (daughter2 == 0 || im + 2 > daughter2) daughter2 = im + 2;
        }
      }
      // For 2 -> 1 hard scatterings reset second daughter to 0.
      if (daughter2 == daughter1) daughter2 = 0;

      // Colour trivial, except reset irrelevant colour indices.
      int colType = particleDataPtr->colType(id);
      int col1   = (colType == 1 || colType == 2 || abs(colType) == 3) 
                 ? lhaUpPtr->col1(iOld) : 0;   
      int col2   = (colType == -1 || colType == 2 || abs(colType) == 3) 
                 ?  lhaUpPtr->col2(iOld) : 0; 

      // Momentum trivial.
      double px  = lhaUpPtr->px(iOld);  
      double py  = lhaUpPtr->py(iOld);  
      double pz  = lhaUpPtr->pz(iOld);  
      double e   = lhaUpPtr->e(iOld);  
      double m   = lhaUpPtr->m(iOld);

      // Polarization
      double pol = lhaUpPtr->spin(iOld);

      // For resonance decay products use resonance mass as scale.
      double scaleNow = scalePr;
      if (mother1 > 4) scaleNow = process[mother1].m();

      // Store Les Houches Accord partons.
      int iNow = process.append( id, status, mother1, mother2, daughter1, 
        daughter2, col1, col2, Vec4(px, py, pz, e), m, scaleNow, pol);

      // Check if need to store lifetime.
      double tau = lhaUpPtr->tau(iOld);
      if (tau > 0.) process[iNow].tau(tau);
    }  
  }

  // Loop through decay chains and set secondary vertices when needed.
  for (int i = 3; i < process.size(); ++i) {
    int iMother  = process[i].mother1();
    
    // If sister to already assigned vertex then assign same.
    if ( process[i - 1].mother1() == iMother && process[i - 1].hasVertex() ) 
      process[i].vProd( process[i - 1].vProd() ); 

    // Else if mother already has vertex and/or lifetime then assign.
    else if ( process[iMother].hasVertex() || process[iMother].tau() > 0.)
      process[i].vProd( process[iMother].vDec() ); 
  }

  // Further info on process. Reset quantities that may or may not be known.
  int    id1Now  =  process[3].id(); 
  int    id2Now  =  process[4].id(); 
  double pdf1    = 0.;
  double pdf2    = 0.;
  double tHat    = 0.;
  double uHat    = 0.;
  double pTHatL  = 0.;
  double m3      = 0.;
  double m4      = 0.;
  double theta   = 0.;
  double phi     = 0.;
  double Q2FacNow, alphaEM, alphaS, Q2Ren, x1Now, x2Now, sHat;

  // Internally generated and stored information.
  if (!isLHA) {
    pdf1         = sigmaProcessPtr->pdf1();
    pdf2         = sigmaProcessPtr->pdf2();
    Q2FacNow     = sigmaProcessPtr->Q2Fac();
    alphaEM      = sigmaProcessPtr->alphaEMRen();
    alphaS       = sigmaProcessPtr->alphaSRen();
    Q2Ren        = sigmaProcessPtr->Q2Ren();
    x1Now        = phaseSpacePtr->x1();
    x2Now        = phaseSpacePtr->x2();
    sHat         = phaseSpacePtr->sHat();
    tHat         = phaseSpacePtr->tHat();
    uHat         = phaseSpacePtr->uHat();
    pTHatL       = phaseSpacePtr->pTHat();
    m3           = phaseSpacePtr->m(3);
    m4           = phaseSpacePtr->m(4);
    theta        = phaseSpacePtr->thetaHat();
    phi          = phaseSpacePtr->phiHat();
  }    

  // Les Houches Accord process partly available, partly to be constructed.
  else {
    Q2FacNow     = (scale < 0.) ? sigmaProcessPtr->Q2Fac() : pow2(scale);
    alphaEM      = lhaUpPtr->alphaQED();
    if (alphaEM < 0.001) alphaEM = sigmaProcessPtr->alphaEMRen();
    alphaS       = lhaUpPtr->alphaQCD();
    if (alphaS  < 0.001) alphaS  = sigmaProcessPtr->alphaSRen();
    Q2Ren        = (scale < 0.) ? sigmaProcessPtr->Q2Ren() : pow2(scale);
    x1Now        = 2. * process[3].e() / infoPtr->eCM();
    x2Now        = 2. * process[4].e() / infoPtr->eCM();
    Vec4 pSum    = process[3].p() + process[4].p();
    sHat         = pSum * pSum;

    // Read info on parton densities if provided.
    if (lhaUpPtr->pdfIsSet()) {
      pdf1       = lhaUpPtr->xpdf1();
      pdf2       = lhaUpPtr->xpdf2();
      Q2FacNow   = pow2(lhaUpPtr->scalePDF());
      x1Now      = lhaUpPtr->x1();
      x2Now      = lhaUpPtr->x2();
    }

    // Reconstruct kinematics of 2 -> 2 processes from momenta.
    if (nFin == 2) {
      Vec4 pDifT = process[3].p() - process[5].p();
      tHat       = pDifT * pDifT;    
      Vec4 pDifU = process[3].p() - process[6].p();
      uHat       = pDifU * pDifU;
      pTHatL     = process[5].pT();
      m3         = process[5].m();    
      m4         = process[6].m(); 
      Vec4 p5    = process[5].p();
      p5.bstback(pSum);
      theta      = p5.theta();   
      phi        = process[5].phi();   
    }
  }

  // Store information.
  if (isHardest) {
    infoPtr->setPDFalpha( id1Now, id2Now, pdf1, pdf2, Q2FacNow, 
      alphaEM, alphaS, Q2Ren);
    infoPtr->setKin( x1Now, x2Now, sHat, tHat, uHat, pTHatL, m3, m4, 
      theta, phi);
  }
  infoPtr->setTypeMPI( code(), pTHatL);

  // For Les Houches event store subprocess classification.
  if (isLHA) {
    int codeSub  = lhaUpPtr->idProcess();
    ostringstream nameSub;
    nameSub << "user process " << codeSub; 
    infoPtr->setSubType( nameSub.str(), codeSub, nFin);
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------
  
// Handle resonance decays.

bool ProcessContainer::decayResonances( Event& process) {

  // Save current event-record size.
  process.saveSize();
  bool physical    = true;
  bool newFlavours = false;

  // Do sequential chain of uncorrelated isotropic decays.
  do {
    physical = resDecaysPtr->next( process);
    if (!physical) return false;

    // Check whether flavours should be correlated.
    // (Currently only relevant for f fbar -> gamma*/Z0 gamma*/Z0.)
    newFlavours = ( sigmaProcessPtr->weightDecayFlav( process) 
                  < rndmPtr->flat() ); 

    // Reset the decay chains if have to redo.
    if (newFlavours) {
      process.restoreSize();
      for (int i = 5; i < process.size(); ++i) process[i].statusPos();
    } 

  // Loop back where required to generate new decays with new flavours.    
  } while (newFlavours);

  // Correct to nonisotropic decays.
  phaseSpacePtr->decayKinematics( process); 

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Reset event generation statistics; but NOT maximum of cross section.

void ProcessContainer::reset() {

  nTry      = 0;
  nSel      = 0;
  nAcc      = 0;
  nTryStat  = 0;
  sigmaSum  = 0.;
  sigma2Sum = 0.;
  sigmaNeg  = 0.;
  sigmaAvg  = 0.;
  sigmaFin  = 0.;
  deltaFin  = 0.;
  wtAccSum  = 0.;

}

//--------------------------------------------------------------------------

// Estimate integrated cross section and its uncertainty.

void ProcessContainer::sigmaDelta() {

  // Initial values. No analysis meaningful unless accepted events.
  nTryStat = nTry;
  sigmaAvg = 0.;
  sigmaFin = 0.;
  deltaFin = 0.;
  if (nAcc == 0) return;

  // Average value. No error analysis unless at least two events.
  double nTryInv  = 1. / nTry;
  double nSelInv  = 1. / nSel;
  double nAccInv  = 1. / nAcc;
  sigmaAvg        = sigmaSum * nTryInv;
  double fracAcc  = nAcc * nSelInv;
  sigmaFin        = sigmaAvg * fracAcc;
  deltaFin        = sigmaFin;
  if (nAcc == 1) return;

  // Estimated error. Quadratic sum of cross section term and
  // binomial from accept/reject step.
  double delta2Sig   = (sigma2Sum * nTryInv - pow2(sigmaAvg)) * nTryInv
    / pow2(sigmaAvg);
  double delta2Veto  = (nSel - nAcc) * nAccInv * nSelInv;
  double delta2Sum   = delta2Sig + delta2Veto;
  deltaFin           = sqrtpos(delta2Sum) * sigmaFin; 

}
 
//==========================================================================

// SetupContainer class.
// Turns list of user-desired processes into a vector of containers.

//--------------------------------------------------------------------------

// Main routine to initialize list of processes.

bool SetupContainers::init(vector<ProcessContainer*>& containerPtrs, 
  Settings& settings, ParticleData* particleDataPtr, Couplings* couplings) {

  // Reset process list, if filled in previous subrun.
  if (containerPtrs.size() > 0) {
    for (int i = 0; i < int(containerPtrs.size()); ++i) 
      delete containerPtrs[i];
    containerPtrs.clear(); 
  }
  SigmaProcess* sigmaPtr;

  // Set up requested objects for soft QCD processes.
  bool softQCD = settings.flag("SoftQCD:all");
  if (softQCD || settings.flag("SoftQCD:minBias")) {
    sigmaPtr = new Sigma0minBias;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (softQCD || settings.flag("SoftQCD:elastic")) {
    sigmaPtr = new Sigma0AB2AB;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (softQCD || settings.flag("SoftQCD:singleDiffractive")) {
    sigmaPtr = new Sigma0AB2XB;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma0AB2AX;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (softQCD || settings.flag("SoftQCD:doubleDiffractive")) {
    sigmaPtr = new Sigma0AB2XX;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for hard QCD processes.
  bool hardQCD = settings.flag("HardQCD:all");
  if (hardQCD || settings.flag("HardQCD:gg2gg")) {
    sigmaPtr = new Sigma2gg2gg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || settings.flag("HardQCD:gg2qqbar")) {
    sigmaPtr = new Sigma2gg2qqbar;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || settings.flag("HardQCD:qg2qg")) {
    sigmaPtr = new Sigma2qg2qg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || settings.flag("HardQCD:qq2qq")) {
    sigmaPtr = new Sigma2qq2qq;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || settings.flag("HardQCD:qqbar2gg")) {
    sigmaPtr = new Sigma2qqbar2gg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || settings.flag("HardQCD:qqbar2qqbarNew")) {
    sigmaPtr = new Sigma2qqbar2qqbarNew;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for c cbar and b bbar, also hard QCD.
  if (hardQCD || settings.flag("HardQCD:gg2ccbar")) {
    sigmaPtr = new Sigma2gg2QQbar(4, 121);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || settings.flag("HardQCD:qqbar2ccbar")) {
    sigmaPtr = new Sigma2qqbar2QQbar(4, 122);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || settings.flag("HardQCD:gg2bbbar")) {
    sigmaPtr = new Sigma2gg2QQbar(5, 123);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || settings.flag("HardQCD:qqbar2bbbar")) {
    sigmaPtr = new Sigma2qqbar2QQbar(5, 124);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
 
  // Set up requested objects for hard QCD 2 -> 3 processes.
  bool hardQCD3parton = settings.flag("HardQCD:3parton");
  if (hardQCD3parton || settings.flag("HardQCD:gg2ggg")) {
    sigmaPtr = new Sigma3gg2ggg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD3parton || settings.flag("HardQCD:qqbar2ggg")) {
    sigmaPtr = new Sigma3qqbar2ggg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD3parton || settings.flag("HardQCD:qg2qgg")) {
    sigmaPtr = new Sigma3qg2qgg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.flag("HardQCD:qq2qqgDiff")) {
    sigmaPtr = new Sigma3qq2qqgDiff;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.flag("HardQCD:qq2qqgSame")) {
    sigmaPtr = new Sigma3qq2qqgSame;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.flag("HardQCD:qqbar2qqbargDiff")) {
    sigmaPtr = new Sigma3qqbar2qqbargDiff;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.flag("HardQCD:qqbar2qqbargSame")) {
    sigmaPtr = new Sigma3qqbar2qqbargSame;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.flag("HardQCD:gg2qqbarg")) {
    sigmaPtr = new Sigma3gg2qqbarg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.flag("HardQCD:qg2qqqbarDiff")) {
    sigmaPtr = new Sigma3qg2qqqbarDiff;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (hardQCD3parton || settings.flag("HardQCD:qg2qqqbarSame")) {
    sigmaPtr = new Sigma3qg2qqqbarSame;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
 
  // Set up requested objects for prompt photon processes.
  bool promptPhotons = settings.flag("PromptPhoton:all");
  if (promptPhotons
    || settings.flag("PromptPhoton:qg2qgamma")) {
    sigmaPtr = new Sigma2qg2qgamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (promptPhotons 
    || settings.flag("PromptPhoton:qqbar2ggamma")) {
    sigmaPtr = new Sigma2qqbar2ggamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (promptPhotons
    || settings.flag("PromptPhoton:gg2ggamma")) {
    sigmaPtr = new Sigma2gg2ggamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (promptPhotons 
    || settings.flag("PromptPhoton:ffbar2gammagamma")) {
    sigmaPtr = new Sigma2ffbar2gammagamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (promptPhotons
    || settings.flag("PromptPhoton:gg2gammagamma")) {
    sigmaPtr = new Sigma2gg2gammagamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Set up requested objects for weak gauge boson t-channel exchange.
  bool weakBosonExchanges = settings.flag("WeakBosonExchange:all");
  if (weakBosonExchanges
    || settings.flag("WeakBosonExchange:ff2ff(t:gmZ)")) {
    sigmaPtr = new Sigma2ff2fftgmZ;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonExchanges
    || settings.flag("WeakBosonExchange:ff2ff(t:W)")) {
    sigmaPtr = new Sigma2ff2fftW;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Set up requested objects for weak gauge boson processes.
  bool weakSingleBosons = settings.flag("WeakSingleBoson:all");
  if (weakSingleBosons
    || settings.flag("WeakSingleBoson:ffbar2gmZ")) {
    sigmaPtr = new Sigma1ffbar2gmZ;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakSingleBosons
    || settings.flag("WeakSingleBoson:ffbar2W")) {
    sigmaPtr = new Sigma1ffbar2W;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Set up requested object for s-channel gamma exchange.
  // Subset of gamma*/Z0 above, intended for multiparton interactions.
  if (settings.flag("WeakSingleBoson:ffbar2ffbar(s:gm)")) {
    sigmaPtr = new Sigma2ffbar2ffbarsgm;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
   
  // Set up requested objects for weak gauge boson pair processes.
  bool weakDoubleBosons = settings.flag("WeakDoubleBoson:all");
  if (weakDoubleBosons
    || settings.flag("WeakDoubleBoson:ffbar2gmZgmZ")) {
    sigmaPtr = new Sigma2ffbar2gmZgmZ;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakDoubleBosons
    || settings.flag("WeakDoubleBoson:ffbar2ZW")) {
    sigmaPtr = new Sigma2ffbar2ZW;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakDoubleBosons
    || settings.flag("WeakDoubleBoson:ffbar2WW")) {
    sigmaPtr = new Sigma2ffbar2WW;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for weak gauge boson + parton processes.
  bool weakBosonAndPartons = settings.flag("WeakBosonAndParton:all");
  if (weakBosonAndPartons
    || settings.flag("WeakBosonAndParton:qqbar2gmZg")) {
    sigmaPtr = new Sigma2qqbar2gmZg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonAndPartons
    || settings.flag("WeakBosonAndParton:qg2gmZq")) {
    sigmaPtr = new Sigma2qg2gmZq;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonAndPartons
    || settings.flag("WeakBosonAndParton:ffbar2gmZgm")) {
    sigmaPtr = new Sigma2ffbar2gmZgm;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonAndPartons
    || settings.flag("WeakBosonAndParton:fgm2gmZf")) {
    sigmaPtr = new Sigma2fgm2gmZf;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonAndPartons
    || settings.flag("WeakBosonAndParton:qqbar2Wg")) {
    sigmaPtr = new Sigma2qqbar2Wg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonAndPartons
    || settings.flag("WeakBosonAndParton:qg2Wq")) {
    sigmaPtr = new Sigma2qg2Wq;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonAndPartons
    || settings.flag("WeakBosonAndParton:ffbar2Wgm")) {
    sigmaPtr = new Sigma2ffbar2Wgm;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonAndPartons
    || settings.flag("WeakBosonAndParton:fgm2Wf")) {
    sigmaPtr = new Sigma2fgm2Wf;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for photon collision processes.
  bool photonCollisions = settings.flag("PhotonCollision:all");
  if (photonCollisions || settings.flag("PhotonCollision:gmgm2qqbar")) {
    sigmaPtr = new Sigma2gmgm2ffbar(1, 261);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (photonCollisions || settings.flag("PhotonCollision:gmgm2ccbar")) {
    sigmaPtr = new Sigma2gmgm2ffbar(4, 262);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (photonCollisions || settings.flag("PhotonCollision:gmgm2bbbar")) {
    sigmaPtr = new Sigma2gmgm2ffbar(5, 263);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (photonCollisions || settings.flag("PhotonCollision:gmgm2ee")) {
    sigmaPtr = new Sigma2gmgm2ffbar(11, 264);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (photonCollisions || settings.flag("PhotonCollision:gmgm2mumu")) {
    sigmaPtr = new Sigma2gmgm2ffbar(13, 265);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (photonCollisions || settings.flag("PhotonCollision:gmgm2tautau")) {
    sigmaPtr = new Sigma2gmgm2ffbar(15, 266);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for charmonium production
  bool charmoniums = settings.flag("Charmonium:all");
  if (charmoniums || settings.flag("Charmonium:gg2QQbar[3S1(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3S11g(4, 401);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:gg2QQbar[3P0(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(4, 0, 402);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:gg2QQbar[3P1(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(4, 1, 403);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:gg2QQbar[3P2(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(4, 2, 404);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:qg2QQbar[3P0(1)]q")) {
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(4, 0, 405);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:qg2QQbar[3P1(1)]q")) {
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(4, 1, 406);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:qg2QQbar[3P2(1)]q")) {
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(4, 2, 407);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:qqbar2QQbar[3P0(1)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(4, 0, 408);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:qqbar2QQbar[3P1(1)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(4, 1, 409);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:qqbar2QQbar[3P2(1)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(4, 2, 410);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:gg2QQbar[3S1(8)]g")) {
    sigmaPtr = new Sigma2gg2QQbarX8g(4, 0, 411);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:gg2QQbar[1S0(8)]g")) {
    sigmaPtr = new Sigma2gg2QQbarX8g(4, 1, 412);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:gg2QQbar[3PJ(8)]g")) {
    sigmaPtr = new Sigma2gg2QQbarX8g(4, 2, 413);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:qg2QQbar[3S1(8)]q")) {
    sigmaPtr = new Sigma2qg2QQbarX8q(4, 0, 414);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:qg2QQbar[1S0(8)]q")) {
    sigmaPtr = new Sigma2qg2QQbarX8q(4, 1, 415);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:qg2QQbar[3PJ(8)]q")) {
    sigmaPtr = new Sigma2qg2QQbarX8q(4, 2, 416);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:qqbar2QQbar[3S1(8)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbarX8g(4, 0, 417);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:qqbar2QQbar[1S0(8)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbarX8g(4, 1, 418);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || settings.flag("Charmonium:qqbar2QQbar[3PJ(8)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbarX8g(4, 2, 419);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
    
  // Set up requested objects for bottomonium production
  bool bottomoniums = settings.flag("Bottomonium:all");
  if (bottomoniums || settings.flag("Bottomonium:gg2QQbar[3S1(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3S11g(5, 501);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:gg2QQbar[3P0(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(5, 0, 502);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:gg2QQbar[3P1(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(5, 1, 503);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:gg2QQbar[3P2(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(5, 2, 504);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:qg2QQbar[3P0(1)]q")) {
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(5, 0, 505);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:qg2QQbar[3P1(1)]q")) {
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(5, 1, 506);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:qg2QQbar[3P2(1)]q")) {
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(5, 2, 507);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:qqbar2QQbar[3P0(1)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(5, 0, 508);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:qqbar2QQbar[3P1(1)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(5, 1, 509);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:qqbar2QQbar[3P2(1)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(5, 2, 510);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:gg2QQbar[3S1(8)]g")) {
    sigmaPtr = new Sigma2gg2QQbarX8g(5, 0, 511);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:gg2QQbar[1S0(8)]g")) {
    sigmaPtr = new Sigma2gg2QQbarX8g(5, 1, 512);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:gg2QQbar[3PJ(8)]g")) {
    sigmaPtr = new Sigma2gg2QQbarX8g(5, 2, 513);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:qg2QQbar[3S1(8)]q")) {
    sigmaPtr = new Sigma2qg2QQbarX8q(5, 0, 514);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:qg2QQbar[1S0(8)]q")) {
    sigmaPtr = new Sigma2qg2QQbarX8q(5, 1, 515);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:qg2QQbar[3PJ(8)]q")) {
    sigmaPtr = new Sigma2qg2QQbarX8q(5, 2, 516);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:qqbar2QQbar[3S1(8)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbarX8g(5, 0, 517);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:qqbar2QQbar[1S0(8)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbarX8g(5, 1, 518);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || settings.flag("Bottomonium:qqbar2QQbar[3PJ(8)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbarX8g(5, 2, 519);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for top production
  bool tops = settings.flag("Top:all");
  if (tops || settings.flag("Top:gg2ttbar")) {
    sigmaPtr = new Sigma2gg2QQbar(6, 601);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (tops || settings.flag("Top:qqbar2ttbar")) {
    sigmaPtr = new Sigma2qqbar2QQbar(6, 602);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (tops || settings.flag("Top:qq2tq(t:W)")) {
    sigmaPtr = new Sigma2qq2QqtW(6, 603);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (tops || settings.flag("Top:ffbar2ttbar(s:gmZ)")) {
    sigmaPtr = new Sigma2ffbar2FFbarsgmZ(6, 604);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (tops || settings.flag("Top:ffbar2tqbar(s:W)")) {
    sigmaPtr = new Sigma2ffbar2FfbarsW(6, 0, 605);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (tops || settings.flag("Top:gmgm2ttbar")) {
    sigmaPtr = new Sigma2gmgm2ffbar(6, 606);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for fourth-generation b' production
  bool bPrimes = settings.flag("FourthBottom:all");
  if (bPrimes || settings.flag("FourthBottom:gg2bPrimebPrimebar")) {
    sigmaPtr = new Sigma2gg2QQbar(7, 801);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bPrimes || settings.flag("FourthBottom:qqbar2bPrimebPrimebar")) {
    sigmaPtr = new Sigma2qqbar2QQbar(7, 802);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bPrimes || settings.flag("FourthBottom:qq2bPrimeq(t:W)")) {
    sigmaPtr = new Sigma2qq2QqtW(7, 803);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bPrimes || settings.flag("FourthBottom:ffbar2bPrimebPrimebar(s:gmZ)")) {
    sigmaPtr = new Sigma2ffbar2FFbarsgmZ(7, 804);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bPrimes || settings.flag("FourthBottom:ffbar2bPrimeqbar(s:W)")) {
    sigmaPtr = new Sigma2ffbar2FfbarsW(7, 0, 805);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bPrimes || settings.flag("FourthBottom:ffbar2bPrimetbar(s:W)")) {
    sigmaPtr = new Sigma2ffbar2FfbarsW(7, 6, 806);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for fourth-generation t' production
  bool tPrimes = settings.flag("FourthTop:all");
  if (tPrimes || settings.flag("FourthTop:gg2tPrimetPrimebar")) {
    sigmaPtr = new Sigma2gg2QQbar(8, 821);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (tPrimes || settings.flag("FourthTop:qqbar2tPrimetPrimebar")) {
    sigmaPtr = new Sigma2qqbar2QQbar(8, 822);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (tPrimes || settings.flag("FourthTop:qq2tPrimeq(t:W)")) {
    sigmaPtr = new Sigma2qq2QqtW(8, 823);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (tPrimes || settings.flag("FourthTop:ffbar2tPrimetPrimebar(s:gmZ)")) {
    sigmaPtr = new Sigma2ffbar2FFbarsgmZ(8, 824);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (tPrimes || settings.flag("FourthTop:ffbar2tPrimeqbar(s:W)")) {
    sigmaPtr = new Sigma2ffbar2FfbarsW(8, 0, 825);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Set up requested objects for two different fourth-generation fermions.
  if (bPrimes || tPrimes 
    || settings.flag("FourthPair:ffbar2tPrimebPrimebar(s:W)")) {
    sigmaPtr = new Sigma2ffbar2FfbarsW(8, 7, 841);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (settings.flag("FourthPair:ffbar2tauPrimenuPrimebar(s:W)")) {
    sigmaPtr = new Sigma2ffbar2FfbarsW(17, 18, 842);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Flag for global choice between SM and BSM Higgses.
  bool useBSMHiggses = settings.flag("Higgs:useBSM");
  
  // Set up requested objects for Standard-Model Higgs production.
  if (!useBSMHiggses) {
    bool HiggsesSM = settings.flag("HiggsSM:all");
    if (HiggsesSM || settings.flag("HiggsSM:ffbar2H")) {
      sigmaPtr = new Sigma1ffbar2H(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (HiggsesSM || settings.flag("HiggsSM:gg2H")) {
     sigmaPtr = new Sigma1gg2H(0);
     containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (HiggsesSM || settings.flag("HiggsSM:gmgm2H")) {
      sigmaPtr = new Sigma1gmgm2H(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (HiggsesSM || settings.flag("HiggsSM:ffbar2HZ")) {
      sigmaPtr = new Sigma2ffbar2HZ(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (HiggsesSM || settings.flag("HiggsSM:ffbar2HW")) {
      sigmaPtr = new Sigma2ffbar2HW(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (HiggsesSM || settings.flag("HiggsSM:ff2Hff(t:ZZ)")) {
      sigmaPtr = new Sigma3ff2HfftZZ(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (HiggsesSM || settings.flag("HiggsSM:ff2Hff(t:WW)")) {
      sigmaPtr = new Sigma3ff2HfftWW(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (HiggsesSM || settings.flag("HiggsSM:gg2Httbar")) {
      sigmaPtr = new Sigma3gg2HQQbar(6,0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (HiggsesSM || settings.flag("HiggsSM:qqbar2Httbar")) {
      sigmaPtr = new Sigma3qqbar2HQQbar(6,0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 

    // Further Standard-Model Higgs processes, not included in "all".
    if (settings.flag("HiggsSM:qg2Hq")) {
      sigmaPtr = new Sigma2qg2Hq(4,0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
      sigmaPtr = new Sigma2qg2Hq(5,0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (settings.flag("HiggsSM:gg2Hbbbar")) {
      sigmaPtr = new Sigma3gg2HQQbar(5,0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (settings.flag("HiggsSM:qqbar2Hbbbar")) {
      sigmaPtr = new Sigma3qqbar2HQQbar(5,0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (settings.flag("HiggsSM:gg2Hg(l:t)")) {
      sigmaPtr = new Sigma2gg2Hglt(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (settings.flag("HiggsSM:qg2Hq(l:t)")) {
      sigmaPtr = new Sigma2qg2Hqlt(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (settings.flag("HiggsSM:qqbar2Hg(l:t)")) {
      sigmaPtr = new Sigma2qqbar2Hglt(0);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
  }

  // Common switch for the group of Higgs production BSM.
  if (useBSMHiggses) {
    bool HiggsesBSM = settings.flag("HiggsBSM:all");    

    // Set up requested objects for BSM H1 production.
    bool HiggsesH1 = settings.flag("HiggsBSM:allH1"); 
    if (HiggsesBSM || HiggsesH1 || settings.flag("HiggsBSM:ffbar2H1")) {
      sigmaPtr = new Sigma1ffbar2H(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );  
    }
    if (HiggsesBSM || HiggsesH1 || settings.flag("HiggsBSM:gg2H1")) {
      sigmaPtr = new Sigma1gg2H(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.flag("HiggsBSM:gmgm2H1")) {
      sigmaPtr = new Sigma1gmgm2H(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.flag("HiggsBSM:ffbar2H1Z")) {
      sigmaPtr = new Sigma2ffbar2HZ(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.flag("HiggsBSM:ffbar2H1W")) {
      sigmaPtr = new Sigma2ffbar2HW(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.flag("HiggsBSM:ff2H1ff(t:ZZ)")) {
      sigmaPtr = new Sigma3ff2HfftZZ(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.flag("HiggsBSM:ff2H1ff(t:WW)")) {
      sigmaPtr = new Sigma3ff2HfftWW(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.flag("HiggsBSM:gg2H1ttbar")) {
      sigmaPtr = new Sigma3gg2HQQbar(6,1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH1 || settings.flag("HiggsBSM:qqbar2H1ttbar")) {
      sigmaPtr = new Sigma3qqbar2HQQbar(6,1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }

    // Further BSM H1 processes, not included in "all".
    if (settings.flag("HiggsBSM:qg2H1q")) {
      sigmaPtr = new Sigma2qg2Hq(4,1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
      sigmaPtr = new Sigma2qg2Hq(5,1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:gg2H1bbbar")) {
      sigmaPtr = new Sigma3gg2HQQbar(5,1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:qqbar2H1bbbar")) {
      sigmaPtr = new Sigma3qqbar2HQQbar(5,1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:gg2H1g(l:t)")) {
      sigmaPtr = new Sigma2gg2Hglt(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:qg2H1q(l:t)")) {
      sigmaPtr = new Sigma2qg2Hqlt(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:qqbar2H1g(l:t)")) {
      sigmaPtr = new Sigma2qqbar2Hglt(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }   

    // Set up requested objects for BSM H2 production.
    bool HiggsesH2 = settings.flag("HiggsBSM:allH2");
    if (HiggsesBSM || HiggsesH2 || settings.flag("HiggsBSM:ffbar2H2")) {
      sigmaPtr = new Sigma1ffbar2H(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (HiggsesBSM || HiggsesH2 || settings.flag("HiggsBSM:gg2H2")) {
      sigmaPtr = new Sigma1gg2H(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.flag("HiggsBSM:gmgm2H2")) {
      sigmaPtr = new Sigma1gmgm2H(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.flag("HiggsBSM:ffbar2H2Z")) {
      sigmaPtr = new Sigma2ffbar2HZ(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.flag("HiggsBSM:ffbar2H2W")) {
      sigmaPtr = new Sigma2ffbar2HW(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.flag("HiggsBSM:ff2H2ff(t:ZZ)")) {
      sigmaPtr = new Sigma3ff2HfftZZ(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.flag("HiggsBSM:ff2H2ff(t:WW)")) {
      sigmaPtr = new Sigma3ff2HfftWW(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.flag("HiggsBSM:gg2H2ttbar")) {
      sigmaPtr = new Sigma3gg2HQQbar(6,2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesH2 || settings.flag("HiggsBSM:qqbar2H2ttbar")) {
      sigmaPtr = new Sigma3qqbar2HQQbar(6,2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }

    // Further BSM H2 processes, not included in "all".
   if (settings.flag("HiggsBSM:qg2H2q")) {
      sigmaPtr = new Sigma2qg2Hq(4,2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
      sigmaPtr = new Sigma2qg2Hq(5,2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:gg2H2bbbar")) {
      sigmaPtr = new Sigma3gg2HQQbar(5,2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:qqbar2H2bbbar")) {
      sigmaPtr = new Sigma3qqbar2HQQbar(5,2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:gg2H2g(l:t)")) {
      sigmaPtr = new Sigma2gg2Hglt(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:qg2H2q(l:t)")) {
      sigmaPtr = new Sigma2qg2Hqlt(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:qqbar2H2g(l:t)")) {
      sigmaPtr = new Sigma2qqbar2Hglt(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }

    // Set up requested objects for BSM A3 production.
    bool HiggsesA3 = settings.flag("HiggsBSM:allA3");
    if (HiggsesBSM || HiggsesA3 || settings.flag("HiggsBSM:ffbar2A3")) {
      sigmaPtr = new Sigma1ffbar2H(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.flag("HiggsBSM:gg2A3")) {
      sigmaPtr = new Sigma1gg2H(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.flag("HiggsBSM:gmgm2A3")) {
      sigmaPtr = new Sigma1gmgm2H(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.flag("HiggsBSM:ffbar2A3Z")) {
      sigmaPtr = new Sigma2ffbar2HZ(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.flag("HiggsBSM:ffbar2A3W")) {
      sigmaPtr = new Sigma2ffbar2HW(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.flag("HiggsBSM:ff2A3ff(t:ZZ)")) {
      sigmaPtr = new Sigma3ff2HfftZZ(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.flag("HiggsBSM:ff2A3ff(t:WW)")) {
      sigmaPtr = new Sigma3ff2HfftWW(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.flag("HiggsBSM:gg2A3ttbar")) {
      sigmaPtr = new Sigma3gg2HQQbar(6,3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (HiggsesBSM || HiggsesA3 || settings.flag("HiggsBSM:qqbar2A3ttbar")) {
      sigmaPtr = new Sigma3qqbar2HQQbar(6,3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }

    // Further BSM A3 processes, not included in "all".
    if (settings.flag("HiggsBSM:qg2A3q")) {
      sigmaPtr = new Sigma2qg2Hq(4,3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
      sigmaPtr = new Sigma2qg2Hq(5,3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:gg2A3bbbar")) {
      sigmaPtr = new Sigma3gg2HQQbar(5,3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:qqbar2A3bbbar")) {
      sigmaPtr = new Sigma3qqbar2HQQbar(5,3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:gg2A3g(l:t)")) {
      sigmaPtr = new Sigma2gg2Hglt(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:qg2A3q(l:t)")) {
      sigmaPtr = new Sigma2qg2Hqlt(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }
    if (settings.flag("HiggsBSM:qqbar2A3g(l:t)")) {
      sigmaPtr = new Sigma2qqbar2Hglt(3);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    }

    // Set up requested objects for Charged Higgs production
    bool HiggsesChg = settings.flag("HiggsBSM:allH+-");
    if (HiggsesBSM || HiggsesChg || settings.flag("HiggsBSM:ffbar2H+-")) {
      sigmaPtr = new Sigma1ffbar2Hchg;
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (HiggsesBSM || HiggsesChg || settings.flag("HiggsBSM:bg2H+-t")) {
      sigmaPtr = new Sigma2qg2Hchgq(6, 1062, "b g -> H+- t");
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 

    // Set up requested objects for Higgs pair-production
    bool HiggsesPairs = settings.flag("HiggsBSM:allHpair");
    if (HiggsesBSM || HiggsesPairs || settings.flag("HiggsBSM:ffbar2A3H1")) {
      sigmaPtr = new Sigma2ffbar2A3H12(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (HiggsesBSM || HiggsesPairs || settings.flag("HiggsBSM:ffbar2A3H2")) {
      sigmaPtr = new Sigma2ffbar2A3H12(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (HiggsesBSM || HiggsesPairs || settings.flag("HiggsBSM:ffbar2H+-H1")) {
      sigmaPtr = new Sigma2ffbar2HchgH12(1);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (HiggsesBSM || HiggsesPairs || settings.flag("HiggsBSM:ffbar2H+-H2")) {
      sigmaPtr = new Sigma2ffbar2HchgH12(2);
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
    if (HiggsesBSM || HiggsesPairs || settings.flag("HiggsBSM:ffbar2H+H-")) {
      sigmaPtr = new Sigma2ffbar2HposHneg();
      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    } 
  }

  // Set up requested objects for SUSY pair processes.
  if(couplings->isSUSY){ 
    CoupSUSY* coupSUSY = (CoupSUSY *) couplings;

    bool SUSYs = settings.flag("SUSY:all");
    bool nmssm = settings.flag("SLHA:NMSSM");

    // Preselected SUSY codes
    int codeA = max( abs(settings.mode("SUSY:idA")),
		     abs(settings.mode("SUSY:idB")));
    int codeB = min( abs(settings.mode("SUSY:idA")),
		     abs(settings.mode("SUSY:idB")));

    // MSSM: 4 neutralinos
    int nNeut = 4;
    if (nmssm) nNeut = 5;
    
    // Gluino-gluino
    if (SUSYs || settings.flag("SUSY:gg2gluinogluino")) {
      // Skip if specific codes not asked for 
      if (codeA == 0 || codeA == 1000021) {
	if (codeB == 0 || codeB == 1000021 ) {
	  sigmaPtr = new Sigma2gg2gluinogluino();
	  containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
	}
      }
    }
    if (SUSYs || settings.flag("SUSY:qqbar2gluinogluino")) {
      // Skip if specific codes not asked for 
      if (codeA == 0 || codeA == 1000021) {
	if (codeB == 0 || codeB == 1000021 ) {
	  sigmaPtr = new Sigma2qqbar2gluinogluino();
	  containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
	}
      }
    }
    
    // Gluino-squark  
    if (SUSYs || settings.flag("SUSY:qg2squarkgluino")) {
      int iproc = 1202; 
      for (int idx = 1; idx <= 6; ++idx) {
	for (int iso = 1; iso <= 2; ++iso) {
	  iproc++;
	  int id3 = iso + ((idx <= 3) 
                  ? 1000000+2*(idx-1) : 2000000+2*(idx-4));
	  int id4 = 1000021;
	  // Skip if specific codes not asked for
	  if (codeA != 0 && codeA != abs(id3) && codeA != abs(id4)) continue;
	  if (codeB != 0 && ( codeA != max(abs(id3),abs(id4)) 
			      || codeB != min(abs(id3),abs(id4)) ) ) continue;
	  sigmaPtr = new Sigma2qg2squarkgluino(id3,iproc);
	  containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
	}
      }
    }
    
    // Squark-antisquark (gg initiated)
    if (SUSYs || settings.flag("SUSY:gg2squarkantisquark")) {
      int iproc = 1214;
      for (int idx = 1; idx <= 6; ++idx) {
	for (int iso = 1; iso <= 2; ++iso) {
	  iproc++;
	  int id = iso + ((idx <= 3) 
                 ? 1000000+2*(idx-1) : 2000000+2*(idx-4)); 	  
	  // Skip if specific codes not asked for
	  if (codeA != 0 && codeA != abs(id)) continue;
	  if (codeA != 0 && codeB != 0 && codeB != abs(id)) continue;
	  sigmaPtr = new Sigma2gg2squarkantisquark(id,iproc);
	  containerPtrs.push_back( new ProcessContainer(sigmaPtr) ); 	
	}
      }
    }
    
    // Squark-antisquark (qqbar initiated)
    if (SUSYs || settings.flag("SUSY:qqbar2squarkantisquark")) {
      int iproc = 1230;
      for (int idx = 1; idx <= 6; ++idx) {
	for (int iso = 1; iso <= 2; ++iso) {
	  for (int jso = iso; jso >= 1; --jso) {
	    for (int jdx = 1; jdx <= 6; ++jdx) {
	      if (iso == jso && jdx < idx) continue;	    
	      int id1 = iso + ((idx <= 3) ? 1000000+2*(idx-1) 
			       : 2000000+2*(idx-4));
	      int id2 = jso + ((jdx <= 3) ? 1000000+2*(jdx-1) 
			       : 2000000+2*(jdx-4));
	      // Update process number counter (for ~q~q, +2 if not self-conj)
	      //if (iproc == 1302) iproc=1310;
	      iproc++;
	      if (iso == jso && id1 != id2) iproc++;
	      // Skip if specific codes not asked for
	      if (codeA != 0 && codeA != abs(id1) 
                && codeA != abs(id2)) continue;
	      if (codeB != 0 && ( codeA != max(abs(id1),abs(id2)) 
                || codeB != min(abs(id1),abs(id2)) ) ) continue;
	      if (iso == jso && id1 != id2) {
		sigmaPtr = new Sigma2qqbar2squarkantisquark(id1,-id2,iproc-1);
		containerPtrs.push_back( new ProcessContainer(sigmaPtr) );   
		sigmaPtr = new Sigma2qqbar2squarkantisquark(id2,-id1,iproc);
		containerPtrs.push_back( new ProcessContainer(sigmaPtr) ); 
	      } else {
		sigmaPtr = new Sigma2qqbar2squarkantisquark(id1,-id2,iproc);
		containerPtrs.push_back( new ProcessContainer(sigmaPtr) );   
	      }
	    }
	  }
	}
      }
    }
    
    // Squark-squark
    if (SUSYs || settings.flag("SUSY:qq2squarksquark")) {
      int iproc = 1350; 
      for (int idx = 1; idx <= 6; ++idx) {
	for (int iso = 1; iso <= 2; ++iso) {
	  for (int jso = iso; jso >= 1; jso--) {
	    for (int jdx = 1; jdx <= 6; ++jdx) {
	      if (iso == jso && jdx < idx) continue;
	      iproc++;
              int id1 = iso + ((idx <= 3) 
                      ? 1000000+2*(idx-1) : 2000000+2*(idx-4));
              int id2 = jso + ((jdx <= 3) 
                      ? 1000000+2*(jdx-1) : 2000000+2*(jdx-4));
	      // Skip if specific codes not asked for
              if (codeA != 0 && codeA != abs(id1) && codeA != abs(id2)) 
                continue;
              if (codeB != 0 && ( codeA != max(abs(id1),abs(id2)) 
                || codeB != min(abs(id1),abs(id2)) ) ) continue;
	      sigmaPtr = new Sigma2qq2squarksquark(id1,id2,iproc);
	      containerPtrs.push_back( new ProcessContainer(sigmaPtr) );    
	    }
	  }
	}
      }
    }
    
    // Neutralino + squark
    if (SUSYs || settings.flag("SUSY:qg2chi0squark")) {
      int iproc = 1430;
      for (int iNeut= 1; iNeut <= nNeut; iNeut++) {
	for (int idx = 1; idx <= 6; idx++) {
	  bool isUp = false;
	  for (int iso = 1; iso <= 2; iso++) {
	    if (iso == 2) isUp = true;
	    iproc++;
	    int id3 = coupSUSY->idNeut(iNeut);
            int id4 = iso + ((idx <= 3) 
                    ? 1000000+2*(idx-1) : 2000000+2*(idx-4));
	    // Skip if specific codes not asked for
	    if (codeA != 0 && codeA != abs(id3) && codeA != abs(id4)) continue;
	    if (codeB != 0 && codeB != min(abs(id3),abs(id4)) ) continue;
	    if (codeA != 0 && codeA == codeB) continue;
	    sigmaPtr = new Sigma2qg2chi0squark(iNeut,idx,isUp,iproc);
	    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
	  }
	}
      }
    }
    
    // Chargino + squark
    if (SUSYs || settings.flag("SUSY:qg2chi+-squark")) {
      int iproc = 1490;
      for (int iChar = 1; iChar <= 2; iChar++) {
	for (int idx = 1; idx <= 6; idx++) {
	  bool isUp = false;
	  for (int iso = 1; iso <= 2; iso++) {
	    if (iso == 2) isUp = true;
	    iproc++;
	    int id3 = coupSUSY->idChar(iChar);
	    int id4 = iso + ((idx <= 3) 
                    ? 1000000+2*(idx-1) : 2000000+2*(idx-4));
	    // Skip if specific codes not asked for
	    if (codeA != 0 && codeA != abs(id3) && codeA != abs(id4)) continue;
	    if (codeB != 0 && codeB != min(abs(id3),abs(id4)) ) continue;
	    if (codeA != 0 && codeA == codeB) continue;
	    sigmaPtr = new Sigma2qg2charsquark(iChar,idx,isUp,iproc);
	    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
	  }
	}
      }
    }
    
    // Neutralino pairs
    if (SUSYs || settings.flag("SUSY:qqbar2chi0chi0")) {
      int iproc = 1550;
      for (int iNeut2 = 1; iNeut2 <= nNeut; iNeut2++) {
	for (int iNeut1 = 1; iNeut1 <= iNeut2; iNeut1++) {
	  iproc++;
	  if (codeA != 0 && codeA != abs(coupSUSY->idNeut(iNeut1)) && 
	      codeA != abs(coupSUSY->idNeut(iNeut2))) continue;
	  if (codeB != 0 && (codeA != max(abs(coupSUSY->idNeut(iNeut1)),
			     abs(coupSUSY->idNeut(iNeut2))) 
                         ||  codeB != min(abs(coupSUSY->idNeut(iNeut1)),
			     abs(coupSUSY->idNeut(iNeut2)))) ) continue;
	  sigmaPtr = new Sigma2qqbar2chi0chi0(iNeut1, iNeut2,iproc);
	  containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
	}
      }
    }
    
    // Neutralino-Chargino
    if (SUSYs || settings.flag("SUSY:qqbar2chi+-chi0")) {
      int iproc = 1570;
      for (int iNeut = 1; iNeut <= nNeut; iNeut++) {
	for (int iChar = 1; iChar <= 2; ++iChar) {
	  iproc += 2;
	  if (codeA != 0 && codeA != coupSUSY->idNeut(iNeut)
	      && codeA != coupSUSY->idChar(iChar)) continue;
	  if (codeB != 0 
	    && ( codeA != max(coupSUSY->idNeut(iNeut),coupSUSY->idChar(iChar))
	      || codeB != min(coupSUSY->idNeut(iNeut),coupSUSY->idChar(iChar))
	      ) ) continue;
	  sigmaPtr = new Sigma2qqbar2charchi0( iChar, iNeut, iproc-1);
	  containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
	  sigmaPtr = new Sigma2qqbar2charchi0(-iChar, iNeut, iproc); 
	  containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
	}
      }
    } 
    
    // Chargino-Chargino
    if (SUSYs || settings.flag("SUSY:qqbar2chi+chi-")) {
      int iproc = 1590;
      for (int i = 1; i <= 2; ++i) {
	for (int j = 1; j <= 2; ++j) {
	  iproc++;	
	  if (codeA != 0 && codeA != abs(coupSUSY->idChar(i)) 
	      && codeA != abs(coupSUSY->idChar(j))) continue;
	  if (codeB != 0 
	      && ( codeA != max(coupSUSY->idChar(i),coupSUSY->idChar(j))
		   || codeB != min(coupSUSY->idChar(i),coupSUSY->idChar(j)) ) ) 
	    continue;
	  sigmaPtr = new Sigma2qqbar2charchar( i,-j, iproc);
	  containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
	}
      }
    } 
    
    // RPV squark production
    if(SUSYs || settings.flag("SUSY:qq2antisquark")) {
      for (int idx = 1; idx <= 6; ++idx) {
	for (int iso = 1; iso <= 2; ++iso) {
	  int id1 = iso + ((idx <= 3) ? 1000000+2*(idx-1) : 2000000+2*(idx-4));
	  if(codeA !=0 && codeA != abs(id1)) continue;
	  sigmaPtr = new Sigma1qq2antisquark(id1);
	  containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
	}
      }
    }
    
  }

  // Set up requested objects for New-Gauge-Boson processes.
  if (settings.flag("NewGaugeBoson:ffbar2gmZZprime")) {
    sigmaPtr = new Sigma1ffbar2gmZZprime();
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (settings.flag("NewGaugeBoson:ffbar2Wprime")) {
    sigmaPtr = new Sigma1ffbar2Wprime();
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (settings.flag("NewGaugeBoson:ffbar2R0")) {
    sigmaPtr = new Sigma1ffbar2Rhorizontal();
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
   
  // Set up requested objects for Left-Right-Symmetry processes.
  bool leftrights = settings.flag("LeftRightSymmmetry:all");
  if (leftrights || settings.flag("LeftRightSymmmetry:ffbar2ZR")) {
    sigmaPtr = new Sigma1ffbar2ZRight();
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leftrights || settings.flag("LeftRightSymmmetry:ffbar2WR")) {
    sigmaPtr = new Sigma1ffbar2WRight();
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leftrights || settings.flag("LeftRightSymmmetry:ll2HL")) {
    sigmaPtr = new Sigma1ll2Hchgchg(1);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leftrights || settings.flag("LeftRightSymmmetry:lgm2HLe")) {
    sigmaPtr = new Sigma2lgm2Hchgchgl(1, 11);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leftrights || settings.flag("LeftRightSymmmetry:lgm2HLmu")) {
    sigmaPtr = new Sigma2lgm2Hchgchgl(1, 13);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leftrights || settings.flag("LeftRightSymmmetry:lgm2HLtau")) {
    sigmaPtr = new Sigma2lgm2Hchgchgl(1, 15);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leftrights || settings.flag("LeftRightSymmmetry:ff2HLff")) {
    sigmaPtr = new Sigma3ff2HchgchgfftWW(1);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leftrights || settings.flag("LeftRightSymmmetry:ffbar2HLHL")) {
    sigmaPtr = new Sigma2ffbar2HchgchgHchgchg(1);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leftrights || settings.flag("LeftRightSymmmetry:ll2HR")) {
    sigmaPtr = new Sigma1ll2Hchgchg(2);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leftrights || settings.flag("LeftRightSymmmetry:lgm2HRe")) {
    sigmaPtr = new Sigma2lgm2Hchgchgl(2, 11);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leftrights || settings.flag("LeftRightSymmmetry:lgm2HRmu")) {
    sigmaPtr = new Sigma2lgm2Hchgchgl(2, 13);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leftrights || settings.flag("LeftRightSymmmetry:lgm2HRtau")) {
    sigmaPtr = new Sigma2lgm2Hchgchgl(2, 15);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leftrights || settings.flag("LeftRightSymmmetry:ff2HRff")) {
    sigmaPtr = new Sigma3ff2HchgchgfftWW(2);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leftrights || settings.flag("LeftRightSymmmetry:ffbar2HRHR")) {
    sigmaPtr = new Sigma2ffbar2HchgchgHchgchg(2);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for leptoquark LQ processes.
  bool leptoquarks = settings.flag("LeptoQuark:all");
  if (leptoquarks || settings.flag("LeptoQuark:ql2LQ")) {
    sigmaPtr = new Sigma1ql2LeptoQuark;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leptoquarks || settings.flag("LeptoQuark:qg2LQl")) {
    sigmaPtr = new Sigma2qg2LeptoQuarkl;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leptoquarks || settings.flag("LeptoQuark:gg2LQLQbar")) {
    sigmaPtr = new Sigma2gg2LQLQbar;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (leptoquarks || settings.flag("LeptoQuark:qqbar2LQLQbar")) {
    sigmaPtr = new Sigma2qqbar2LQLQbar;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for excited-fermion processes.
  bool excitedfermions = settings.flag("ExcitedFermion:all");
  if (excitedfermions || settings.flag("ExcitedFermion:dg2dStar")) {
    sigmaPtr = new Sigma1qg2qStar(1);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:ug2uStar")) {
    sigmaPtr = new Sigma1qg2qStar(2);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:sg2sStar")) {
    sigmaPtr = new Sigma1qg2qStar(3);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:cg2cStar")) {
    sigmaPtr = new Sigma1qg2qStar(4);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:bg2bStar")) {
    sigmaPtr = new Sigma1qg2qStar(5);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:egm2eStar")) {
    sigmaPtr = new Sigma1lgm2lStar(11);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:mugm2muStar")) {
    sigmaPtr = new Sigma1lgm2lStar(13);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:taugm2tauStar")) {
    sigmaPtr = new Sigma1lgm2lStar(15);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:qq2dStarq")) {
    sigmaPtr = new Sigma2qq2qStarq(1);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:qq2uStarq")) {
    sigmaPtr = new Sigma2qq2qStarq(2);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:qq2sStarq")) {
    sigmaPtr = new Sigma2qq2qStarq(3);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:qq2cStarq")) {
    sigmaPtr = new Sigma2qq2qStarq(4);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:qq2bStarq")) {
    sigmaPtr = new Sigma2qq2qStarq(5);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:qqbar2eStare")) {
    sigmaPtr = new Sigma2qqbar2lStarlbar(11);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:qqbar2nueStarnue")) {
    sigmaPtr = new Sigma2qqbar2lStarlbar(12);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:qqbar2muStarmu")) {
    sigmaPtr = new Sigma2qqbar2lStarlbar(13);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:qqbar2numuStarnumu")) {
    sigmaPtr = new Sigma2qqbar2lStarlbar(14);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions || settings.flag("ExcitedFermion:qqbar2tauStartau")) {
    sigmaPtr = new Sigma2qqbar2lStarlbar(15);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (excitedfermions 
    || settings.flag("ExcitedFermion:qqbar2nutauStarnutau")) {
    sigmaPtr = new Sigma2qqbar2lStarlbar(16);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Set up requested objects for contact interaction processes.
  if (settings.flag("ContactInteractions:QCqq2qq")) {
    sigmaPtr = new Sigma2QCqq2qq();
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (settings.flag("ContactInteractions:QCqqbar2qqbar")) {
    sigmaPtr = new Sigma2QCqqbar2qqbar();
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (settings.flag("ContactInteractions:QCffbar2eebar")) {
    sigmaPtr = new Sigma2QCffbar2llbar(11, 4003);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ContactInteractions:QCffbar2mumubar")) {
    sigmaPtr = new Sigma2QCffbar2llbar(13, 4004);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ContactInteractions:QCffbar2tautaubar")) {
    sigmaPtr = new Sigma2QCffbar2llbar(15, 4005);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Set spin of particles in the Hidden Valley scenario.
  int spinFv = settings.mode("HiddenValley:spinFv");
  for (int i = 1; i < 7; ++i) {
    if (particleDataPtr->spinType( 4900000 + i) != spinFv + 1) 
        particleDataPtr->spinType( 4900000 + i,    spinFv + 1);
    if (particleDataPtr->spinType( 4900010 + i) != spinFv + 1) 
        particleDataPtr->spinType( 4900010 + i,    spinFv + 1);
  }
  if (spinFv != 1) {
    if (particleDataPtr->spinType( 4900101) != 2) 
       particleDataPtr->spinType( 4900101, 2);
  } else {
    int spinqv = settings.mode("HiddenValley:spinqv");
    if (particleDataPtr->spinType( 4900101) != 2 * spinqv + 1) 
        particleDataPtr->spinType( 4900101,    2 * spinqv + 1);
  }
  
  // Set up requested objects for HiddenValley processes.
  bool hiddenvalleys = settings.flag("HiddenValley:all");
  if (hiddenvalleys || settings.flag("HiddenValley:gg2DvDvbar")) {
    sigmaPtr = new Sigma2gg2qGqGbar( 4900001, 4901, spinFv, 
      "g g -> Dv Dvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:gg2UvUvbar")) {
    sigmaPtr = new Sigma2gg2qGqGbar( 4900002, 4902, spinFv, 
      "g g -> Uv Uvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:gg2SvSvbar")) {
    sigmaPtr = new Sigma2gg2qGqGbar( 4900003, 4903, spinFv, 
      "g g -> Sv Svbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:gg2CvCvbar")) {
    sigmaPtr = new Sigma2gg2qGqGbar( 4900004, 4904, spinFv, 
      "g g -> Cv Cvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:gg2BvBvbar")) {
    sigmaPtr = new Sigma2gg2qGqGbar( 4900005, 4905, spinFv, 
      "g g -> Bv Bvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:gg2TvTvbar")) {
    sigmaPtr = new Sigma2gg2qGqGbar( 4900006, 4906, spinFv, 
      "g g -> Tv Tvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:qqbar2DvDvbar")) {
    sigmaPtr = new Sigma2qqbar2qGqGbar( 4900001, 4911, spinFv, 
      "q qbar -> Dv Dvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:qqbar2UvUvbar")) {
    sigmaPtr = new Sigma2qqbar2qGqGbar( 4900002, 4912, spinFv, 
      "q qbar -> Uv Uvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:qqbar2SvSvbar")) {
    sigmaPtr = new Sigma2qqbar2qGqGbar( 4900003, 4913, spinFv, 
      "q qbar -> Sv Svbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:qqbar2CvCvbar")) {
    sigmaPtr = new Sigma2qqbar2qGqGbar( 4900004, 4914, spinFv, 
      "q qbar -> Cv Cvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:qqbar2BvBvbar")) {
    sigmaPtr = new Sigma2qqbar2qGqGbar( 4900005, 4915, spinFv, 
      "q qbar -> Bv Bvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:qqbar2TvTvbar")) {
    sigmaPtr = new Sigma2qqbar2qGqGbar( 4900006, 4916, spinFv, 
      "q qbar -> Tv Tvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:ffbar2DvDvbar")) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900001, 4921, spinFv, 
      "f fbar -> Dv Dvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:ffbar2UvUvbar")) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900002, 4922, spinFv, 
      "f fbar -> Uv Uvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:ffbar2SvSvbar")) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900003, 4923, spinFv, 
      "f fbar -> Sv Svbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:ffbar2CvCvbar")) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900004, 4924, spinFv, 
      "f fbar -> Cv Cvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:ffbar2BvBvbar")) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900005, 4925, spinFv, 
      "f fbar -> Bv Bvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:ffbar2TvTvbar")) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900006, 4926, spinFv, 
      "f fbar -> Tv Tvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:ffbar2EvEvbar")) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900011, 4931, spinFv, 
      "f fbar -> Ev Evbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:ffbar2nuEvnuEvbar")) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900012, 4932, spinFv, 
      "f fbar -> nuEv nuEvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:ffbar2MUvMUvbar")) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900013, 4933, spinFv, 
      "f fbar -> MUv MUvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:ffbar2nuMUvnuMUvbar")) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900014, 4934, spinFv, 
      "f fbar -> nuMUv nuMUvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:ffbar2TAUvTAUvbar")) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900015, 4935, spinFv, 
      "f fbar -> TAUv TAUvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:ffbar2nuTAUvnuTAUvbar")) {
    sigmaPtr = new Sigma2ffbar2fGfGbar( 4900016, 4936, spinFv, 
      "f fbar -> nuTAUv nuTAUvbar");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hiddenvalleys || settings.flag("HiddenValley:ffbar2Zv")) {
    sigmaPtr = new Sigma1ffbar2Zv();
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for RS extra-dimensional G* processes.
  bool extraDimGstars = settings.flag("ExtraDimensionsG*:all");
  if (extraDimGstars || settings.flag("ExtraDimensionsG*:gg2G*")) {
    sigmaPtr = new Sigma1gg2GravitonStar;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (extraDimGstars || settings.flag("ExtraDimensionsG*:ffbar2G*")) {
    sigmaPtr = new Sigma1ffbar2GravitonStar;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (settings.flag("ExtraDimensionsG*:gg2G*g")) {
    sigmaPtr = new Sigma2gg2GravitonStarg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (settings.flag("ExtraDimensionsG*:qg2G*q")) {
    sigmaPtr = new Sigma2qg2GravitonStarq;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (settings.flag("ExtraDimensionsG*:qqbar2G*g")) {
    sigmaPtr = new Sigma2qqbar2GravitonStarg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  //  Set up requested objects for RS extra-dimensional KKgluon processes.
  if (settings.flag("ExtraDimensionsG*:qqbar2KKgluon*")) {
    sigmaPtr = new Sigma1qqbar2KKgluonStar;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // NOAM: Set up requested objects for TEV extra-dimensional processes.
  if (settings.flag("ExtraDimensionsTEV:ffbar2ddbar")) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(1, 5061); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }  
  if (settings.flag("ExtraDimensionsTEV:ffbar2uubar")) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(2, 5062); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }  
  if (settings.flag("ExtraDimensionsTEV:ffbar2ssbar")) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(3, 5063); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }  
  if (settings.flag("ExtraDimensionsTEV:ffbar2ccbar")) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(4, 5064); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }  
  if (settings.flag("ExtraDimensionsTEV:ffbar2bbbar")) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(5, 5065); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }  
  if (settings.flag("ExtraDimensionsTEV:ffbar2ttbar")) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(6, 5066); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }  
  if (settings.flag("ExtraDimensionsTEV:ffbar2e+e-")) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(11, 5071); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsTEV:ffbar2nuenuebar")) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(12, 5072); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsTEV:ffbar2mu+mu-")) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(13, 5073); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsTEV:ffbar2numunumubar")) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(14, 5074); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsTEV:ffbar2tau+tau-")) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(15, 5075); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsTEV:ffbar2nutaunutaubar")) {
    sigmaPtr = new Sigma2ffbar2TEVffbar(16, 5076); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for large extra-dimensional G processes.
  bool extraDimLEDmono = settings.flag("ExtraDimensionsLED:monojet");
  if (extraDimLEDmono || settings.flag("ExtraDimensionsLED:gg2Gg")) {
    sigmaPtr = new Sigma2gg2LEDUnparticleg( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimLEDmono || settings.flag("ExtraDimensionsLED:qg2Gq")) {
    sigmaPtr = new Sigma2qg2LEDUnparticleq( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimLEDmono || settings.flag("ExtraDimensionsLED:qqbar2Gg")) {
    sigmaPtr = new Sigma2qqbar2LEDUnparticleg( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsLED:ffbar2GZ")) {
    sigmaPtr = new Sigma2ffbar2LEDUnparticleZ( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsLED:ffbar2Ggamma")) {
    sigmaPtr = new Sigma2ffbar2LEDUnparticlegamma( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsLED:ffbar2gammagamma")) {
    sigmaPtr = new Sigma2ffbar2LEDgammagamma( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsLED:gg2gammagamma")) {
    sigmaPtr = new Sigma2gg2LEDgammagamma( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsLED:ffbar2llbar")) {
    sigmaPtr = new Sigma2ffbar2LEDllbar( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsLED:gg2llbar")) {
    sigmaPtr = new Sigma2gg2LEDllbar( true );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  bool extraDimLEDdij = settings.flag("ExtraDimensionsLED:dijets");
  if (extraDimLEDdij || settings.flag("ExtraDimensionsLED:gg2DJgg")) {
    sigmaPtr = new Sigma2gg2LEDgg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimLEDdij || settings.flag("ExtraDimensionsLED:gg2DJqqbar")) {
    sigmaPtr = new Sigma2gg2LEDqqbar;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimLEDdij || settings.flag("ExtraDimensionsLED:qg2DJqg")) {
    sigmaPtr = new Sigma2qg2LEDqg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimLEDdij || settings.flag("ExtraDimensionsLED:qq2DJqq")) {
    sigmaPtr = new Sigma2qq2LEDqq;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimLEDdij || settings.flag("ExtraDimensionsLED:qqbar2DJgg")) {
    sigmaPtr = new Sigma2qqbar2LEDgg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimLEDdij || settings.flag("ExtraDimensionsLED:qqbar2DJqqbarNew")) {
    sigmaPtr = new Sigma2qqbar2LEDqqbarNew;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Set up requested objects for unparticle processes.
  bool extraDimUnpartmono = settings.flag("ExtraDimensionsUnpart:monojet");
  if (extraDimUnpartmono || settings.flag("ExtraDimensionsUnpart:gg2Ug")) {
    sigmaPtr = new Sigma2gg2LEDUnparticleg( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimUnpartmono || settings.flag("ExtraDimensionsUnpart:qg2Uq")) {
    sigmaPtr = new Sigma2qg2LEDUnparticleq( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (extraDimUnpartmono || settings.flag("ExtraDimensionsUnpart:qqbar2Ug")) {
    sigmaPtr = new Sigma2qqbar2LEDUnparticleg( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsUnpart:ffbar2UZ")) {
    sigmaPtr = new Sigma2ffbar2LEDUnparticleZ( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsUnpart:ffbar2Ugamma")) {
    sigmaPtr = new Sigma2ffbar2LEDUnparticlegamma( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsUnpart:ffbar2gammagamma")) {
    sigmaPtr = new Sigma2ffbar2LEDgammagamma( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsUnpart:gg2gammagamma")) {
    sigmaPtr = new Sigma2gg2LEDgammagamma( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsUnpart:ffbar2llbar")) {
    sigmaPtr = new Sigma2ffbar2LEDllbar( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }
  if (settings.flag("ExtraDimensionsUnpart:gg2llbar")) {
    sigmaPtr = new Sigma2gg2LEDllbar( false );
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Done. 
  return true;

}

//--------------------------------------------------------------------------

// Routine to initialize list of second hard processes.

bool SetupContainers::init2(vector<ProcessContainer*>& container2Ptrs, 
  Settings& settings) {

  // Reset process list, if filled in previous subrun.
  if (container2Ptrs.size() > 0) {
    for (int i = 0; i < int(container2Ptrs.size()); ++i) 
      delete container2Ptrs[i];
    container2Ptrs.clear(); 
  }
  SigmaProcess* sigmaPtr;

  // Two hard QCD jets.
  if (settings.flag("SecondHard:TwoJets")) {
    sigmaPtr = new Sigma2gg2gg;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2qqbar;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2qg;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qq2qq;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2gg;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2qqbarNew;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbar(4, 121);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar(4, 122);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbar(5, 123);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar(5, 124);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // A prompt photon and a hard jet.
  if (settings.flag("SecondHard:PhotonAndJet")) {
    sigmaPtr = new Sigma2qg2qgamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2ggamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2ggamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Two prompt photons.
  if (settings.flag("SecondHard:TwoPhotons")) {
    sigmaPtr = new Sigma2ffbar2gammagamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2gammagamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
 
  // Charmonium production.
  if (settings.flag("SecondHard:Charmonium")) {
    sigmaPtr = new Sigma2gg2QQbar3S11g(4, 401);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(4, 0, 402);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(4, 1, 403);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(4, 2, 404);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(4, 0, 405);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(4, 1, 406);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(4, 2, 407);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(4, 0, 408);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(4, 1, 409);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(4, 2, 410);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbarX8g(4, 0, 411);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbarX8g(4, 1, 412);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbarX8g(4, 2, 413);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2QQbarX8q(4, 0, 414);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2QQbarX8q(4, 1, 415);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2QQbarX8q(4, 2, 416);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbarX8g(4, 0, 417);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbarX8g(4, 1, 418);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbarX8g(4, 2, 419);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
    
  // Bottomonium production.
  if (settings.flag("SecondHard:Bottomonium")) {
    sigmaPtr = new Sigma2gg2QQbar3S11g(5, 501);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(5, 0, 502);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(5, 1, 503);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(5, 2, 504);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(5, 0, 505);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(5, 1, 506);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(5, 2, 507);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(5, 0, 508);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(5, 1, 509);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(5, 2, 510);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbarX8g(5, 0, 511);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbarX8g(5, 1, 512);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbarX8g(5, 2, 513);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2QQbarX8q(5, 0, 514);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2QQbarX8q(5, 1, 515);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2QQbarX8q(5, 2, 516);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbarX8g(5, 0, 517);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbarX8g(5, 1, 518);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbarX8g(5, 2, 519);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // A single gamma*/Z0.
  if (settings.flag("SecondHard:SingleGmZ")) {
    sigmaPtr = new Sigma1ffbar2gmZ;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // A single W+-.
  if (settings.flag("SecondHard:SingleW")) {
    sigmaPtr = new Sigma1ffbar2W;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // A gamma*/Z0 and a hard jet.
  if (settings.flag("SecondHard:GmZAndJet")) {
    sigmaPtr = new Sigma2qqbar2gmZg;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2gmZq;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // A W+- and a hard jet.
  if (settings.flag("SecondHard:WAndJet")) {
    sigmaPtr = new Sigma2qqbar2Wg;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2Wq;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Top pair production.
  if (settings.flag("SecondHard:TopPair")) {
    sigmaPtr = new Sigma2gg2QQbar(6, 601);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar(6, 602);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2ffbar2FFbarsgmZ(6, 604);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Single top production.
  if (settings.flag("SecondHard:SingleTop")) {
    sigmaPtr = new Sigma2qq2QqtW(6, 603);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2ffbar2FfbarsW(6, 0, 605);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Two b jets - already part of TwoJets sample above.
  if (settings.flag("SecondHard:TwoBJets")) {
    sigmaPtr = new Sigma2gg2QQbar(5, 123);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar(5, 124);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Done. 
  return true;

}

//==========================================================================

} // end namespace Pythia8
