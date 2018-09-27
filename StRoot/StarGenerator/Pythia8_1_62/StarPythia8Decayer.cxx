#include "StarPythia8Decayer.h"
#include "TLorentzVector.h"
#include <assert.h>
#include "TParticle.h"
#include "StarGenerator/UTIL/StarRandom.h"
#include "TSystem.h"
#include "StMessMgr.h"

#ifndef Pythia8_version 
#error "Pythia8_version is not defined"
#endif

// ----------------------------------------------------------------------------
// Remap pythia's random number generator to StarRandom
class PyRand : public Pythia8::RndmEngine {
public:
  Double_t flat() { return StarRandom::Instance().flat(); }
};

//..................................................................................................
StarPythia8Decayer::StarPythia8Decayer( Pythia8::Pythia *pythia ) : mPythia(pythia), mOwner(false), mDebug(0)
{

  if ( pythia ) return; // If provided, we expect pythia to be correctly initialized


  // Create the new instance of pythia, specifying the path to the
  // auxilliary files required by pythia. 
  // NOTE:  When adding new versions of Pythia8, we need to specify
  // the version of the code in the source code
  TString path = "StRoot/StarGenerator/"; path+= Pythia8_version; path+="/xmldoc/"; 
  { 
    ifstream in(path); 
    if (!in.good()) { path = "$(STAR)/"+path; }
    path = gSystem->ExpandPathName(path.Data());
  }
  
  Info(GetName(),Form("MC version is %s data at %s",Pythia8_version,path.Data()));
  Info(GetName(),Form("Configuration files at %s",path.Data()));

  mPythia = new Pythia8::Pythia( path.Data() );
  mPythia -> setRndmEnginePtr( new PyRand() );
  mPythia->readString("SoftQCD:elastic = on");
  
  // Initialize for pp @ 510.0 GeV
  mPythia->init( 2212, 2212, 510.0 );
    
}
//..................................................................................................
void StarPythia8Decayer::Init()
{
  
}
//..................................................................................................
void StarPythia8Decayer::Decay( Int_t pdgid, TLorentzVector *_p )
{

  LOG_INFO << "Decay pdgid=" << pdgid << endm;

  // Clear the event from the last run
  ClearEvent();
  // Add the particle to the pythia stack
  AppendParticle( pdgid, _p );
  // Get the particle ID and allow it to decay
  Int_t id = mPythia -> event[0].id();
  mPythia->particleData.mayDecay( id, true );
  mPythia->moreDecays();
  // Print list of pythia lines on debug
  if ( mDebug ) mPythia->event.list();
}
//..................................................................................................
Int_t StarPythia8Decayer::ImportParticles( TClonesArray *_array )
{
  // Save the decay products
  assert(_array);
  TClonesArray &array = *_array;
  array.Clear();

  Int_t nparts = 0;

  for ( Int_t i=0; i<mPythia->event.size();i++ )
    {
      
     // if ( mPythia->event[i].id() == 90 ) continue; //??
      new(array[nparts++]) TParticle ( 		  
          mPythia->event[i].id(),
	  mPythia->event[i].status(),
	  mPythia->event[i].mother1(),
	  mPythia->event[i].mother2(),
	  mPythia->event[i].daughter1(),
	  mPythia->event[i].daughter2(),
	  mPythia->event[i].px(),       // [GeV/c]
	  mPythia->event[i].py(),       // [GeV/c]
	  mPythia->event[i].pz(),       // [GeV/c]
	  mPythia->event[i].e(),        // [GeV]
	  mPythia->event[i].xProd(),    // [mm]
	  mPythia->event[i].yProd(),    // [mm]
	  mPythia->event[i].zProd(),    // [mm]
	  mPythia->event[i].tProd());   // [mm/c]		  
    };
  return nparts;
}
//..................................................................................................
// not implemented.  complain about it.  loudly.
void StarPythia8Decayer::SetForceDecay( Int_t type ){ assert(0); }
void StarPythia8Decayer::ForceDecay(){ assert(0); }
Float_t StarPythia8Decayer::GetPartialBranchingRatio( Int_t ipdg ){ assert(0); }
void StarPythia8Decayer::ReadDecayTable(){ assert(0); }
//..................................................................................................
Float_t StarPythia8Decayer::GetLifetime(Int_t pdg) 
{
  // return lifetime in seconds of teh particle with PDG number pdg
  return (mPythia->particleData.tau0(pdg) * 3.3333e-12) ;
}
//..................................................................................................
void StarPythia8Decayer::AppendParticle(Int_t pdg, TLorentzVector* p)
{
  // Append a particle to the stack
  mPythia->event.append(pdg, 11, 0, 0, p->Px(), p->Py(), p->Pz(), p->E(), p->M());
}
//..................................................................................................
void StarPythia8Decayer::ClearEvent()
{
  mPythia->event.clear();
}
//..................................................................................................
StarPythia8Decayer::~StarPythia8Decayer()
{
  if (mPythia) delete mPythia; // delete 
}
