// BeamParticle.cc is a part of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// BeamParticle class.

#include "Pythia8/BeamParticle.h"

namespace Pythia8 {
 
//==========================================================================

// The BeamParticle class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// A lepton that takes (almost) the full beam energy does not leave a remnant.
const double BeamParticle::XMINUNRESOLVED = 1. - 1e-10;

//--------------------------------------------------------------------------

// Initialize data on a beam particle and save pointers.

void BeamParticle::init( int idIn, double pzIn, double eIn, double mIn,
  Info* infoPtrIn, Settings& settings, ParticleData* particleDataPtrIn,
  Rndm* rndmPtrIn,PDF* pdfInPtr, PDF* pdfHardInPtr, bool isUnresolvedIn,
  StringFlav* flavSelPtrIn) {

  // Store input pointers (and one bool) for future use.
  infoPtr           = infoPtrIn;
  particleDataPtr   = particleDataPtrIn;
  rndmPtr           = rndmPtrIn;
  pdfBeamPtr        = pdfInPtr;
  pdfHardBeamPtr    = pdfHardInPtr;
  isUnresolvedBeam  = isUnresolvedIn;
  flavSelPtr        = flavSelPtrIn;

  // Maximum quark kind in allowed incoming beam hadrons.
  maxValQuark       = settings.mode("BeamRemnants:maxValQuark");

  // Power of (1-x)^power/sqrt(x) for remnant valence quark distribution.
  valencePowerMeson = settings.parm("BeamRemnants:valencePowerMeson");
  valencePowerUinP  = settings.parm("BeamRemnants:valencePowerUinP");
  valencePowerDinP  = settings.parm("BeamRemnants:valencePowerDinP");

  // Enhancement factor of x of diquark.
  valenceDiqEnhance = settings.parm("BeamRemnants:valenceDiqEnhance");

  // Assume g(x) ~ (1-x)^power/x to constrain companion to sea quark.
  companionPower    = settings.mode("BeamRemnants:companionPower");

  // Assume g(x) ~ (1-x)^power/x to constrain companion to sea quark.
  companionPower    = settings.mode("BeamRemnants:companionPower");

  // Allow or not more than one valence quark to be kicked out.
  allowJunction     = settings.flag("BeamRemnants:allowJunction");

  // For low-mass diffractive system kick out q/g = norm / mass^power.
  pickQuarkNorm     = settings.parm("Diffraction:pickQuarkNorm");
  pickQuarkPower    = settings.parm("Diffraction:pickQuarkPower");

  // Width of primordial kT distribution in low-mass diffractive systems.
  diffPrimKTwidth   = settings.parm("Diffraction:primKTwidth");

  // Suppress large masses of beam remnant in low-mass diffractive systems.
  diffLargeMassSuppress = settings.parm("Diffraction:largeMassSuppress");

  // Store info on the incoming beam.
  idBeam            = idIn;
  initBeamKind();
  pBeam             = Vec4( 0., 0., pzIn, eIn);
  mBeam             = mIn;

}

//--------------------------------------------------------------------------

// Initialize kind and valence flavour content of incoming beam.
// For recognized hadrons one can generate multiparton interactions.
// Dynamic choice of meson valence flavours in newValenceContent below.

void BeamParticle::initBeamKind() {

  // Reset.
  idBeamAbs    = abs(idBeam);
  isLeptonBeam = false;
  isHadronBeam = false;
  isMesonBeam  = false;
  isBaryonBeam = false;
  nValKinds    = 0;

  // Check for leptons.
  if (idBeamAbs > 10 && idBeamAbs < 17) {
    nValKinds = 1;
    nVal[0]   = 1;
    idVal[0]  = idBeam;
    isLeptonBeam = true;
  }

  //  Done if cannot be lowest-lying hadron state.
  if (idBeamAbs < 101 || idBeamAbs > 9999) return;

  // Resolve valence content for assumed Pomeron.
  if (idBeamAbs == 990) {
    isMesonBeam = true;
    nValKinds = 2;
    nVal[0]   = 1 ;
    nVal[1]   = 1 ;
    newValenceContent();
  
  // Resolve valence content for assumed meson. Flunk unallowed codes.
  } else if (idBeamAbs < 1000) {
    int id1 = idBeamAbs/100;
    int id2 = (idBeamAbs/10)%10;
    if ( id1 < 1 || id1 > maxValQuark
      || id2 < 1 || id2 > maxValQuark ) return;
    isMesonBeam = true;
    
    // Store valence content of a confirmed meson.
    nValKinds = 2;
    nVal[0]   = 1 ;
    nVal[1]   = 1;
    if (id1%2 == 0) {
      idVal[0] = id1;
      idVal[1] = -id2;
    } else {
      idVal[0] = id2;
      idVal[1] = -id1;
    }
    newValenceContent();
  
  // Resolve valence content for assumed baryon. Flunk unallowed codes.
  } else {
    int id1 = idBeamAbs/1000;
    int id2 = (idBeamAbs/100)%10;
    int id3 = (idBeamAbs/10)%10;
    if ( id1 < 1 || id1 > maxValQuark || id2 < 1 || id2 > maxValQuark
      || id3 < 1 || id3 > maxValQuark) return;
    if (id2 > id1 || id3 > id1) return;
    isBaryonBeam = true;

    // Store valence content of a confirmed baryon.
    nValKinds = 1; idVal[0] = id1; nVal[0] = 1;
    if (id2 == id1) ++nVal[0];
    else {
      nValKinds = 2;
      idVal[1]  = id2;
      nVal[1]   = 1;
    }
    if (id3 == id1) ++nVal[0];
    else if (id3 == id2) ++nVal[1];
    else {
      idVal[nValKinds] = id3;
      nVal[nValKinds]  = 1;
      ++nValKinds;
    }
  }
  
  // Flip flavours for antimeson or antibaryon, and then done.
  if (idBeam < 0) for (int i = 0; i < nValKinds; ++i) idVal[i] = -idVal[i];
  isHadronBeam = true;
  Q2ValFracSav = -1.;

}

//--------------------------------------------------------------------------


// Dynamic choice of meson valence flavours for pi0, K0S, K0L, Pomeron.

void BeamParticle::newValenceContent() {

  // A pi0 oscillates between d dbar and u ubar.
  if (idBeam == 111) {
    idVal[0] = (rndmPtr->flat() < 0.5) ? 1 : 2;
    idVal[1] = -idVal[0];

  // A K0S or K0L oscillates between d sbar and s dbar.
  } else if (idBeam == 130 || idBeam == 310) {
    idVal[0] = (rndmPtr->flat() < 0.5) ?  1 :  3;
    idVal[1] = (idVal[0] == 1)      ? -3 : -1;

  // For a Pomeron split gluon remnant into d dbar or u ubar.
  } else if (idBeam == 990) {
    idVal[0] = (rndmPtr->flat() < 0.5) ? 1 : 2;
    idVal[1] = -idVal[0];

  // Other hadrons so far do not require any event-by-event change.
  } else return;

  // Propagate change to PDF routine(s).
  pdfBeamPtr->newValenceContent( idVal[0], idVal[1]);
  if (pdfHardBeamPtr != pdfBeamPtr && pdfHardBeamPtr != 0)
    pdfHardBeamPtr->newValenceContent( idVal[0], idVal[1]);

}

//--------------------------------------------------------------------------

double BeamParticle::xMax(int iSkip) {

  // Minimum requirement on remaining energy > nominal mass for hadron.
  double xLeft = 1.;
  if (isHadron()) xLeft -= m() / e();
  if (size() == 0) return xLeft;

  // Subtract what was carried away by initiators (to date).
  for (int i = 0; i < size(); ++i)
    if (i != iSkip && resolved[i].isFromBeam()) xLeft -= resolved[i].x();
  return xLeft;

}

//--------------------------------------------------------------------------

// Parton distributions, reshaped to take into account previous
// multiparton interactions. By picking a non-negative iSkip value,
// one particular interaction is skipped, as needed for ISR

double BeamParticle::xfModified(int iSkip, int idIn, double x, double Q2) {

  // Initial values.
  idSave    = idIn;
  iSkipSave = iSkip;
  xqVal     = 0.;
  xqgSea    = 0.;
  xqCompSum = 0.;

  // Fast procedure for first interaction.
  if (size() == 0) {
    if (x >= 1.) return 0.;
    bool canBeVal = false;
    for (int i = 0; i < nValKinds; ++i)
      if (idIn == idVal[i]) canBeVal = true;
    if (canBeVal) {
      xqVal     = xfVal( idIn, x, Q2);
      xqgSea    = xfSea( idIn, x, Q2);
    }
    else xqgSea = xf( idIn, x, Q2);

  // More complicated procedure for non-first interaction.
  } else {

    // Sum up the x already removed, and check that remaining x is enough.
    double xUsed = 0.;
    for (int i = 0; i < size(); ++i)
      if (i != iSkip) xUsed += resolved[i].x();
    double xLeft = 1. - xUsed;
    if (x >= xLeft) return 0.;
    double xRescaled = x / xLeft;

    // Calculate total and remaining amount of x carried by valence quarks.
    double xValTot = 0.;
    double xValLeft = 0.;
    for (int i = 0; i < nValKinds; ++i) {
      nValLeft[i] = nVal[i];
      for (int j = 0; j < size(); ++j)
      if (j != iSkip && resolved[j].isValence()
        && resolved[j].id() == idVal[i]) --nValLeft[i];
      double xValNow =  xValFrac(i, Q2);
      xValTot += nVal[i] * xValNow;
      xValLeft += nValLeft[i] * xValNow;
    }

    // Calculate total amount of x carried by unmatched companion quarks.
    double xCompAdded = 0.;
    for (int i = 0; i < size(); ++i)
    if (i != iSkip && resolved[i].isUnmatched()) xCompAdded
      += xCompFrac( resolved[i].x() / (xLeft + resolved[i].x()) )
      // Typo warning: extrafactor missing in Skands&Sjostrand article;
      // <x> for companion refers to fraction of x left INCLUDING sea quark.
      // To be modified further??
      * (1. + resolved[i].x() / xLeft);
  
    // Calculate total rescaling factor and pdf for sea and gluon.
    double rescaleGS = max( 0., (1. - xValLeft - xCompAdded)
      / (1. - xValTot) );
    xqgSea = rescaleGS * xfSea( idIn, xRescaled, Q2);

    // Find valence part and rescale it to remaining number of quarks.
    for (int i = 0; i < nValKinds; ++i)
    if (idIn == idVal[i] && nValLeft[i] > 0)
      xqVal = xfVal( idIn, xRescaled, Q2)
      * double(nValLeft[i]) / double(nVal[i]);
                                                                               
    // Find companion part, by adding all companion contributions.
    for (int i = 0; i < size(); ++i)
    if (i != iSkip && resolved[i].id() == -idIn
      && resolved[i].isUnmatched()) {
      double xsRescaled = resolved[i].x() / (xLeft + resolved[i].x());
      double xcRescaled = x / (xLeft + resolved[i].x());
      double xqCompNow = xCompDist( xcRescaled, xsRescaled);
      resolved[i].xqCompanion( xqCompNow);
      xqCompSum += xqCompNow;
    }
  }

  // Add total, but only return relevant part for ISR. More cases??
  // Watch out, e.g. g can come from either kind of quark.??
  xqgTot = xqVal + xqgSea + xqCompSum;
  if (iSkip >= 0) {
    if (resolved[iSkip].isValence()) return xqVal;
    if (resolved[iSkip].isUnmatched()) return xqgSea + xqCompSum;
  }
  return xqgTot;
  
}

//--------------------------------------------------------------------------

// Decide whether a quark extracted from the beam is of valence, sea or
// companion kind; in the latter case also pick its companion.
// Assumes xfModified has already been called.

int BeamParticle::pickValSeaComp() {

  // If parton already has a companion than reset code for this.
  int oldCompanion = resolved[iSkipSave].companion();
  if (oldCompanion >= 0) resolved[oldCompanion].companion(-2);

  // Default assignment is sea.
  int vsc = -2;

  // For gluons or photons no sense of valence or sea.
  if (idSave == 21 || idSave == 22) vsc = -1;

  // For lepton beam assume same-kind lepton inside is valence.
  else if (isLeptonBeam && idSave == idBeam) vsc = -3;

  // Decide if valence or sea quark.
  else {
    double xqRndm = xqgTot * rndmPtr->flat();
    if (xqRndm < xqVal) vsc = -3;
    else if (xqRndm < xqVal + xqgSea) vsc = -2;
 
    // If not either, loop over all possible companion quarks.
    else {
      xqRndm -= xqVal + xqgSea;
      for (int i = 0; i < size(); ++i)
      if (i != iSkipSave && resolved[i].id() == -idSave
        && resolved[i].isUnmatched()) {
        xqRndm -= resolved[i].xqCompanion();
        if (xqRndm < 0.) vsc = i;
        break;
      }
    }
  }

  // Bookkeep assignment; for sea--companion pair both ways.
  resolved[iSkipSave].companion(vsc);
  if (vsc >= 0) resolved[vsc].companion(iSkipSave);
  
  // Done; return code for choice (to distinguish valence/sea in Info).
  return vsc;

}

//--------------------------------------------------------------------------

// Fraction of hadron momentum sitting in a valence quark distribution.
// Based on hardcoded parametrizations of CTEQ 5L numbers.

double BeamParticle::xValFrac(int j, double Q2) {

  // Only recalculate when required.
  if (Q2 != Q2ValFracSav) {
    Q2ValFracSav = Q2;
     
    // Q2-dependence of log-log form; assume fixed Lambda = 0.2.
    double llQ2 = log( log( max( 1., Q2) / 0.04 ));

    // Fractions carried by u and d in proton.
    uValInt =  0.48 / (1. + 1.56 * llQ2);
    dValInt = 0.385 / (1. + 1.60 * llQ2);
  }

  // Baryon with three different quark kinds: (2 * u + d) / 3 of proton.
  if (isBaryonBeam && nValKinds == 3) return (2. * uValInt + dValInt) / 3.;

  // Baryon with one or two identical: like d or u of proton.
  if (isBaryonBeam && nVal[j] == 1) return dValInt;
  if (isBaryonBeam && nVal[j] == 2) return uValInt;

  // Meson: (2 * u + d) / 2 of proton so same total valence quark fraction.
    return 0.5 * (2. * uValInt + dValInt);

}

//--------------------------------------------------------------------------

// The momentum integral of a companion quark, with its partner at x_s,
// using an approximate gluon density like (1 - x_g)^power / x_g.
// The value corresponds to an unrescaled range between 0 and 1 - x_s.

double BeamParticle::xCompFrac(double xs) {

  // Select case by power of gluon (1-x_g) shape.
  switch (companionPower) {

    case 0:
       return xs * ( 5. + xs * (-9. - 2. * xs * (-3. + xs)) + 3. * log(xs) )
         / ( (-1. + xs) * (2. + xs * (-1. + 2. * xs)) );

    case 1:
       return -1. -3. * xs + ( 2. * pow2(-1. + xs) * (1. + xs + xs*xs))
         / ( 2. + xs*xs * (xs - 3.) + 3. * xs * log(xs) );

    case 2:
       return xs * ( (1. - xs) * (19. + xs * (43. + 4. * xs))
         + 6. * log(xs) * (1. + 6. * xs + 4.*xs*xs) ) /
        ( 4. * ( (xs - 1.) * (1. + xs * (4. + xs) )
        - 3. * xs * log(xs) * (1 + xs) ) );

    case 3:
      return 3. * xs * ( (xs - 1.) * (7. + xs * (28. + 13. * xs))
        - 2. * log(xs) * (1. + xs * (9. + 2. * xs * (6. + xs))) )
        / ( 4. + 27. * xs - 31. * pow3(xs)
        + 6. * xs * log(xs) * (3. + 2. * xs * (3.+xs)) );

    default:
      return ( -9. * xs * (xs*xs - 1.) * (5. + xs * (24. + xs)) + 12. * xs
        * log(xs) * (1. + 2. * xs) * (1. + 2. * xs * (5. + 2. * xs)) )
        / ( 8. * (1. + 2. * xs) * ((xs - 1.) * (1. + xs * (10. + xs))
        - 6. * xs * log(xs) * (1. + xs)) );

  }
}

//--------------------------------------------------------------------------

// The x*f pdf of a companion quark at x_c, with its sea partner at x_s,
// using an approximate gluon density like (1 - x_g)^power / x_g.
// The value corresponds to an unrescaled range between 0 and 1 - x_s.

double BeamParticle::xCompDist(double xc, double xs) {

  // Mother gluon momentum fraction. Check physical limit.
  double xg = xc + xs;
  if (xg > 1.) return 0.;

  // Common factor, including splitting kernel and part of gluon density
  // (and that it is x_c * f that is coded).
  double fac = 3. * xc * xs * (xc*xc + xs*xs) / pow4(xg);

  // Select case by power of gluon (1-x_g) shape.
  switch (companionPower) {

    case 0:
      return fac / ( 2. - xs * (3. - xs * (3. - 2. * xs)) );

    case 1:
      return fac * (1. - xg) / ( 2. + xs*xs * (-3. + xs) + 3. * xs * log(xs) );

    case 2:
      return fac * pow2(1. - xg) / ( 2. * ((1. - xs) * (1. + xs * (4. + xs))
        + 3. * xs * (1. + xs) * log(xs)) );

    case 3:
      return fac * pow3(1. - xg) * 2. / ( 4. + 27. * xs - 31. * pow3(xs)
        + 6. * xs * log(xs) * (3. + 2. * xs * (3. + xs)) );

    default:
       return fac * pow4(1. - xg) / ( 2. * (1. + 2. * xs) * ((1. - xs)
         * (1. + xs * (10. + xs)) + 6. * xs * log(xs) * (1. + xs)) );

  }
}

//--------------------------------------------------------------------------

// Add required extra remnant flavour content. Also initial colours.

bool BeamParticle::remnantFlavours(Event& event) {

  // A baryon will have a junction, unless a diquark is formed later.
  hasJunctionBeam = (isBaryon());

  // Store how many hard-scattering partons were removed from beam.
  nInit = size();

  // Find remaining valence quarks.
  for (int i = 0; i < nValKinds; ++i) {
    nValLeft[i] = nVal[i];
    for (int j = 0; j < nInit; ++j) if (resolved[j].isValence()
      && resolved[j].id() == idVal[i]) --nValLeft[i];
    // Add remaining valence quarks to record. Partly temporary values.
    for (int k = 0; k < nValLeft[i]; ++k) append(0, idVal[i], 0., -3);
  }

  // If at least two valence quarks left in baryon remnant then form diquark.
  int nInitPlusVal = size();
  if (isBaryon() && nInitPlusVal - nInit >= 2) {

    // If three, pick two at random to form diquark, else trivial.
    int iQ1 = nInit;
    int iQ2 = nInit + 1;
    if (nInitPlusVal - nInit == 3) {
      double pickDq = 3. * rndmPtr->flat();
      if (pickDq > 1.) iQ2 = nInit + 2;
      if (pickDq > 2.) iQ1 = nInit + 1;
    }

    // Pick spin 0 or 1 according to SU(6) wave function factors.
    int idDq = flavSelPtr->makeDiquark( resolved[iQ1].id(),
      resolved[iQ2].id(), idBeam);

    // Overwrite with diquark flavour and remove one slot. No more junction.
    resolved[iQ1].id(idDq);
    if (nInitPlusVal - nInit == 3 && iQ2 == nInit + 1)
      resolved[nInit + 1].id( resolved[nInit + 2].id() );
    resolved.pop_back();
    hasJunctionBeam = false;
  }

  // Find companion quarks to unmatched sea quarks.
  for (int i = 0; i < nInit; ++i)
  if (resolved[i].isUnmatched()) {

    // Add companion quark to record; and bookkeep both ways.
    append(0, -resolved[i].id(), 0., i);
    resolved[i].companion(size() - 1);
  }

  // If no other remnants found, add a gluon or photon to carry momentum.
  if (size() == nInit) {
    int    idRemnant = (isHadronBeam) ? 21 : 22;
    append(0, idRemnant, 0., -1);
  }

  // Set initiator and remnant masses.
  for (int i = 0; i < size(); ++i) {
    if (i < nInit) resolved[i].m(0.);
    else resolved[i].m( particleDataPtr->m0( resolved[i].id() ) );
  }

  // For debug purposes: reject beams with resolved junction topology.
  if (hasJunctionBeam && !allowJunction) return false;

  // Pick initial colours for remnants.
  for (int i = nInit; i < size(); ++i) {
    int colType = particleDataPtr->colType( resolved[i].id() );
    int col = (colType == 1 || colType == 2) ? event.nextColTag() : 0;
    int acol = (colType == -1 || colType == 2) ? event.nextColTag() : 0;
    resolved[i].cols( col, acol);
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Correlate all initiators and remnants to make a colour singlet.

bool BeamParticle::remnantColours(Event& event, vector<int>& colFrom,
  vector<int>& colTo) {

  // No colours in lepton beams so no need to do anything.
  if (isLeptonBeam) return true;

  // Copy initiator colour info from the event record to the beam.
  for (int i = 0; i < size(); ++i) {
    int j =  resolved[i].iPos();
    resolved[i].cols( event[j].col(), event[j].acol());
  }

  // Find number and position of valence quarks, of gluons, and
  // of sea-companion pairs (counted as gluons) in the beam remnants.
  // Skip gluons with same colour as anticolour and rescattering partons.
  vector<int> iVal;
  vector<int> iGlu;
  for (int i = 0; i < size(); ++i)
  if (resolved[i].isFromBeam()) {
    if ( resolved[i].isValence() ) iVal.push_back(i);
    else if ( resolved[i].isCompanion() && resolved[i].companion() > i )
      iGlu.push_back(i);
    else if ( resolved[i].id() == 21
      && resolved[i].col() != resolved[i].acol() ) iGlu.push_back(i);
  }
      
  // Pick a valence quark to which gluons are attached.
  // Do not resolve quarks in diquark. (More sophisticated??)
  int iValSel= iVal[0];
  if (iVal.size() == 2) {
    if ( abs(resolved[iValSel].id()) > 10 ) iValSel = iVal[1];
  } else {
    double rndmValSel = 3. * rndmPtr->flat();
    if (rndmValSel > 1.) iValSel= iVal[1];
    if (rndmValSel > 2.) iValSel= iVal[2];
  }

  // This valence quark defines initial (anti)colour.
  int iBeg = iValSel;
  bool hasCol = (resolved[iBeg].col() > 0);
  int begCol = (hasCol) ? resolved[iBeg].col() : resolved[iBeg].acol();

  // Do random stepping through gluon/(sea+companion) list.
  vector<int> iGluRndm;
  for (int i = 0; i < int(iGlu.size()); ++i)
    iGluRndm.push_back( iGlu[i] );
  for (int iOrder = 0; iOrder < int(iGlu.size()); ++iOrder) {
    int iRndm = int( double(iGluRndm.size()) * rndmPtr->flat());
    int iGluSel = iGluRndm[iRndm];
    iGluRndm[iRndm] = iGluRndm[iGluRndm.size() - 1];
    iGluRndm.pop_back();

    // Find matching anticolour/colour to current colour/anticolour.
    int iEnd = iGluSel;
    int endCol = (hasCol) ? resolved[iEnd].acol() : resolved[iEnd].col();
    // Not gluon but sea+companion pair: go to other.
    if (endCol == 0) {
      iEnd = resolved[iEnd].companion();
      endCol = (hasCol) ? resolved[iEnd].acol() : resolved[iEnd].col();
    }

    // Collapse this colour-anticolour pair to the lowest one.
    if (begCol < endCol) {
      if (hasCol) resolved[iEnd].acol(begCol);
      else resolved[iEnd].col(begCol);
      colFrom.push_back(endCol);
      colTo.push_back(begCol);
    } else {
      if (hasCol) resolved[iBeg].col(endCol);
      else resolved[iBeg].acol(endCol);
      colFrom.push_back(begCol);
      colTo.push_back(endCol);
    }

    // Pick up the other colour of the recent gluon and repeat.
    iBeg = iEnd;
    begCol = (hasCol) ? resolved[iBeg].col() : resolved[iBeg].acol();
    // Not gluon but sea+companion pair: go to other.
    if (begCol == 0) {
      iBeg = resolved[iBeg].companion();
      begCol = (hasCol) ? resolved[iBeg].col() : resolved[iBeg].acol();
    }

  // At end of gluon/(sea+companion) list.
  }
 
  // Now begin checks, and also finding junction information.
  // Loop through remnant partons; isolate all colours and anticolours.
  vector<int> colList;
  vector<int> acolList;
  for (int i = 0; i < size(); ++i)
  if ( resolved[i].isFromBeam() )
  if ( resolved[i].col() != resolved[i].acol() ) {
    if (resolved[i].col() > 0) colList.push_back( resolved[i].col() );
    if (resolved[i].acol() > 0) acolList.push_back( resolved[i].acol() );
  }

  // Remove all matching colour-anticolour pairs.
  bool foundPair = true;
  while (foundPair && colList.size() > 0 && acolList.size() > 0) {
    foundPair = false;
    for (int iCol = 0; iCol < int(colList.size()); ++iCol) {
      for (int iAcol = 0; iAcol < int(acolList.size()); ++iAcol) {
        if (acolList[iAcol] == colList[iCol]) {
          colList[iCol] = colList.back();
          colList.pop_back();
          acolList[iAcol] = acolList.back();
          acolList.pop_back();
          foundPair = true;
          break;
        }
      } if (foundPair) break;
    }
  }

  // Usually one unmatched pair left to collapse.
  if (colList.size() == 1 && acolList.size() == 1) {
    int finalFrom = max( colList[0], acolList[0]);
    int finalTo   = min( colList[0], acolList[0]);
    for (int i = 0; i < size(); ++i)
    if ( resolved[i].isFromBeam() ) {
      if (resolved[i].col()  == finalFrom) resolved[i].col(finalTo);
      if (resolved[i].acol() == finalFrom) resolved[i].acol(finalTo);
    }
    colFrom.push_back(finalFrom);
    colTo.push_back(finalTo);

  // Store an (anti)junction when three (anti)coloured daughters.
  } else if (hasJunctionBeam && colList.size() == 3
    && acolList.size() == 0) {
    event.appendJunction( 1, colList[0], colList[1], colList[2]);
    junCol[0] = colList[0];
    junCol[1] = colList[1];
    junCol[2] = colList[2];
  } else if (hasJunctionBeam && acolList.size() == 3
    && colList.size() == 0) {
    event.appendJunction( 2, acolList[0], acolList[1], acolList[2]);
    junCol[0] = acolList[0];
    junCol[1] = acolList[1];
    junCol[2] = acolList[2];

  // Any other nonvanishing values indicate failure.
  } else if (colList.size() > 0 || acolList.size() > 0) {
    infoPtr->errorMsg("Error in BeamParticle::remnantColours: "
      "leftover unmatched colours");
    return false;
  }

  // Store colour assignment of beam particles.
  for (int i = nInit; i < size(); ++i)
    event[resolved[i].iPos()].cols( resolved[i].col(), resolved[i].acol() );

  // Done.
  return true;
}


//--------------------------------------------------------------------------

// Pick unrescaled x values for beam remnant sharing.

double BeamParticle::xRemnant( int i) {

  double x = 0.;

  // Calculation of x of valence quark or diquark, for latter as sum.
  if (resolved[i].isValence()) {

    // Resolve diquark into sum of two quarks.
    int id1 = resolved[i].id();
    int id2 = 0;
    if (abs(id1) > 10) {
      id2 = (id1 > 0) ? (id1/100)%10 : -(((-id1)/100)%10);
      id1 = (id1 > 0) ? id1/1000 : -((-id1)/1000);
    }
 
    // Loop over (up to) two quarks; add their contributions.
    for (int iId = 0; iId < 2; ++iId) {
      int idNow = (iId == 0) ? id1 : id2;
      if (idNow == 0) break;
      double xPart = 0.;

      // Assume form (1-x)^a / sqrt(x).
      double xPow = valencePowerMeson;
      if (isBaryonBeam) {
        if (nValKinds == 3 || nValKinds == 1)
          xPow = (3. * rndmPtr->flat() < 2.)
            ? valencePowerUinP : valencePowerDinP ;
        else if (nValence(idNow) == 2) xPow = valencePowerUinP;
        else xPow = valencePowerDinP;
      }
      do xPart = pow2( rndmPtr->flat() );
      while ( pow(1. - xPart, xPow) < rndmPtr->flat() );

      // End loop over (up to) two quarks. Possibly enhancement for diquarks.
      x += xPart;
    }
   if (id2 != 0) x *= valenceDiqEnhance;
      
  // Calculation of x of sea quark, based on companion association.
  } else if (resolved[i].isCompanion()) {

    // Find rescaled x value of companion.
    double xLeft = 1.;
    for (int iInit = 0; iInit < nInit; ++iInit)
      if (resolved[iInit].isFromBeam()) xLeft -= resolved[iInit].x();
    double xCompanion = resolved[ resolved[i].companion() ].x();
    xCompanion /= (xLeft + xCompanion);

    // Now use ansatz q(x; x_c) < N/(x +x_c) to pick x.
    do x = pow( xCompanion, rndmPtr->flat()) - xCompanion;
    while ( pow( (1. - x - xCompanion) / (1. - xCompanion), companionPower)
      * (pow2(x) + pow2(xCompanion)) / pow2(x + xCompanion)
      < rndmPtr->flat() );

  // Else, rarely, a single gluon remnant, so value does not matter.
  } else x = 1.;
  return x;

}
   
//--------------------------------------------------------------------------

// Print the list of resolved partons in a beam.

void BeamParticle::list(ostream& os) const {

  // Header.
  os << "\n --------  PYTHIA Partons resolved in beam  -----------------"
     << "-------------------------------------------------------------\n"
     << "\n    i  iPos      id       x    comp   xqcomp    pTfact      "
     << "colours      p_x        p_y        p_z         e          m \n";
  
  // Loop over list of removed partons and print it.
  double xSum  = 0.;
  Vec4   pSum;
  for (int i = 0; i < size(); ++i) {
    ResolvedParton res = resolved[i];
    os << fixed << setprecision(6) << setw(5) << i << setw(6) << res.iPos()
       << setw(8) << res.id() << setw(10) << res.x() << setw(6)
       << res.companion() << setw(10) << res.xqCompanion() << setw(10)
       << res.pTfactor() << setprecision(3) << setw(6) << res.col()
       << setw(6) << res.acol() << setw(11) << res.px() << setw(11)
       << res.py() << setw(11) << res.pz() << setw(11) << res.e()
       << setw(11) << res.m() << "\n";

    // Also find sum of x and p values.
    if (res.companion() != -10) {
      xSum  += res.x();
      pSum  += res.p();
    }
  }

  // Print sum and endline.
  os << setprecision(6) << "             x sum:" << setw(10) << xSum
     << setprecision(3) << "                                p sum:"
     << setw(11) << pSum.px() << setw(11) << pSum.py() << setw(11)
     << pSum.pz() << setw(11) << pSum.e()
     << "\n\n --------  End PYTHIA Partons resolved in beam  -----------"
     << "---------------------------------------------------------------"
     << endl;
}
   
//--------------------------------------------------------------------------

// Test whether a lepton is to be considered as unresolved.

bool BeamParticle::isUnresolvedLepton() {
 
  // Require record to consist of lepton with full energy plus a photon.
  if (!isLeptonBeam || resolved.size() > 2 || resolved[1].id() != 22
    || resolved[0].x() < XMINUNRESOLVED) return false;
  return true;
  
}
   
//--------------------------------------------------------------------------

// For a diffractive system, decide whether to kick out gluon or quark.

bool BeamParticle::pickGluon(double mDiff) {
  
  // Relative weight to pick a quark, assumed falling with energy.
  double probPickQuark = pickQuarkNorm / pow( mDiff, pickQuarkPower);
  return  ( (1. + probPickQuark) * rndmPtr->flat() < 1. );
  
}
   
//--------------------------------------------------------------------------

// Pick a valence quark at random. (Used for diffractive systems.)

int BeamParticle::pickValence() {

  // Pick one valence quark at random.
  int nTotVal = (isBaryonBeam) ? 3 : 2;
  double rnVal = rndmPtr->flat() * nTotVal;
  int iVal = (rnVal < 1.) ? 1 : ( (rnVal < 2.) ? 2 : 3 );

  // This valence in slot 1, the rest thereafter.
  idVal1 = 0;
  idVal2 = 0;
  idVal3 = 0;
  int iNow = 0;
  for (int i = 0; i < nValKinds; ++i)
  for (int j = 0; j < nVal[i]; ++j) {
    ++iNow;
    if (iNow == iVal) idVal1 = idVal[i];
    else if ( idVal2 == 0) idVal2 = idVal[i];
    else idVal3 = idVal[i];
  }

  // Construct diquark if baryon.
  if (idVal3 != 0) idVal2 = flavSelPtr->makeDiquark( idVal2, idVal3, idBeam);

  // Done.
  return idVal1;

}
   
//--------------------------------------------------------------------------

// Share lightcone momentum between two remnants in a diffractive system.

double BeamParticle::zShare( double mDiff, double m1, double m2) {

  // Set up as valence in normal beam so can use xRemnant code.
  append(0, idVal1, 0., -3);
  append(0, idVal2, 0., -3);
  double m2Diff = mDiff*mDiff;

  // Begin to generate z and pT until acceptable solution.
  double wtAcc = 0.;
  do {
    double x1 = xRemnant(0);
    double x2 = xRemnant(0);
    zRel = x1 / (x1 + x2);
    pair<double, double> gauss2 = rndmPtr->gauss2();
    pxRel = diffPrimKTwidth * gauss2.first;
    pyRel = diffPrimKTwidth * gauss2.second;

    // Suppress large invariant masses of remnant system.
    double mTS1 = m1*m1 + pxRel*pxRel + pyRel*pyRel;
    double mTS2 = m2*m2 + pxRel*pxRel + pyRel*pyRel;
    double m2Sys = mTS1 / zRel + mTS2 / (1. - zRel);
    wtAcc = (m2Sys < m2Diff)
      ? pow( 1. - m2Sys / m2Diff, diffLargeMassSuppress) : 0.;
  } while (wtAcc < rndmPtr->flat());

  // Done.
  return zRel;

}

//==========================================================================

} // end namespace Pythia8
