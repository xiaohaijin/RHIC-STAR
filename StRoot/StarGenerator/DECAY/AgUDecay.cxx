#include "AgUDecay.h"
#include "StMessMgr.h"
#include "TVirtualMCDecayer.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "StarDecayManager.h"
#include "St_geant_Maker/St_geant_Maker.h"
#include "StarGenerator/UTIL/StarParticleData.h"

#define geant3 St_geant_Maker::instance()->Geant3()

StarParticleData &pdb = StarParticleData::instance();

AgUDecay AgUDecay::sInstance;
//
// --------------------------------------------------------------------------------------------------
//
extern "C" {

  void agudcay_() 
  { static AgUDecay &decay = AgUDecay::instance();

    // Decay current particle on G3 stack
    decay();

  };

  void gsking_(int &igk);

}
//
void gsking( int igk ){ gsking_(igk); }
//
// --------------------------------------------------------------------------------------------------
//


Int_t AgUDecay::operator()()
{
   
  Gctrak_t& gctrak = *(geant3->Gctrak()); // kinematics of current track
  Gcking_t& gcking = *(geant3->Gcking()); // kinematics of decay products
  Gckin3_t& gckin3 = *(geant3->Gckin3()); // vertex of decay products

  float x = gctrak.vect[0];
  float y = gctrak.vect[1];
  float z = gctrak.vect[2];

  //  LOG_INFO << Form(">>> decay() called x=%f y=%f z=%f <<<",x,y,z) << endm;
  if (0==mDecayer) return 0; // no decayer registerd

  //  int np = 0;

  mArray   -> Clear();
  mDecayer -> ForceDecay();

  int idGeant3 = geant3->Gckine()->ipart;

  double pmom = double( gctrak.vect[6] );
  double px   = double( gctrak.vect[3] ) * pmom;
  double py   = double( gctrak.vect[4] ) * pmom;
  double pz   = double( gctrak.vect[5] ) * pmom;
  double E    = double( gctrak.getot   ); // Input energy

  mP[0] = px; mP[1] = py; mP[2] = pz; mP[3] = E;

  // Extract PDG ID from idPart
  int idPdg = pdb.GetParticleG3( idGeant3 )->PdgCode();

  // Perform the decay
  mDecayer -> Decay( idPdg, &mP );

  // Retrieve the particles into the clones array
  int np = mDecayer -> ImportParticles( mArray ); if ( np<1 ) return np;

  // 
  // Loop over daughter particles in the decay and stack them for
  // transport.  (Potentially recurses through the entire decay
  // chain until it reaches particles known to the simulator).
  //
  TParticle* mother = (TParticle*)mArray->At(0);
  // Possible that the decay is setup in a special "system" residing at zero, so 
  // start by searching for the PDG id we are decaying
  for ( int i=0;i<np;i++ ) {
      mother = (TParticle*)mArray->At(i);  if ( mother->GetPdgCode() == idPdg ) break; 
  }
  int first = mother->GetFirstDaughter();
  int last  = mother->GetLastDaughter();
  double EnergySum = 0; // Energy conservation...
  for ( int i=first /* first daughter */; i <= last; i++ )
    {

      TParticle    *particle    = (TParticle *)mArray->At(i);
      if ( 0==particle ) {
          LOG_WARN << "Daughter " << i << " not valid for particle " << Form("%s [@0x%p]",mother->GetName(),mother) << endm;
	  mArray->Print();
	  continue;
      }

      EnergySum += StackParticleForTransport( particle );

    }

  double violation;
  if ( violation = TMath::Abs(E - EnergySum)/E > 0.1E-5 ) {
      TParticle    *particle    = (TParticle *)mArray->At(0);
  LOG_WARN << particle->GetName() << " decay violates E conservation by"
  << violation*100 <<"%"
  <<endm;
  mNonConservation++;
  }

  return np;
}
//_____________________________________________________________________________
bool AgUDecay::MayTransport( const TParticle* particle )
{

      int           first       = particle->GetFirstDaughter();
      int           last        = particle->GetLastDaughter();
      int           pdgid       = particle->GetPdgCode();
      int           status      = particle->GetStatusCode();
      TParticlePDG *particlePDG = pdb.GetParticle(pdgid); 
      int           g3id        = particlePDG->TrackingCode();

      if ( 0 == g3id )              // particle not known to G3
      switch( mDiscovery ) {
          case kDecay: return false;
	  case kSpawn: 
	     pdb.AddParticleToG3( particlePDG, mNextG3id++ ); 
	     assert(mNextG3id < 60000);
	     return true;
	     break;
	  default:
	     assert(0);
	     break;
      }
      return true;



}
//_____________________________________________________________________________
double AgUDecay::StackParticleForTransport( const TParticle* particle )
{

  Gctrak_t& gctrak = *(geant3->Gctrak()); // kinematics of current track
  Gcking_t& gcking = *(geant3->Gcking()); // kinematics of decay products
  Gckin3_t& gckin3 = *(geant3->Gckin3()); // vertex of decay products

  double EnergySum = 0.0; // returns total energy of stacked particles

      int           first       = particle->GetFirstDaughter();
      int           last        = particle->GetLastDaughter();
      int           pdgid       = particle->GetPdgCode();
      int           status      = particle->GetStatusCode();
      TParticlePDG *particlePDG = pdb.GetParticle(pdgid); 
      int           g3id        = particlePDG->TrackingCode();

      // Stack the particle for transport by geant if possible,
      // otherwise recursively stack daughters.                       
      if ( 0 == MayTransport(particle) )
      {
      LOG_INFO << Form("%s [@0x%p] decayed in place", particle->GetName(),particle) << endm;
           for ( int kid=first; kid<=last; kid++ )
	   {
	        TParticle* daughter = (TParticle*)mArray->At( kid ); 
      LOG_INFO << Form("    %s [@x%p] visit kid: %s [@0x%p]", particle->GetName(),particle,daughter->GetName(),daughter) << endm;
	        EnergySum += StackParticleForTransport( daughter );
	   }
	   return EnergySum; 
      }

      // Stack this particle for transport 
      LOG_INFO << Form("%s [@0x%p] stacked for transport", particle->GetName(),particle) << endm;

      int &index = gcking.ngkine;

      // Throw particle on the stack
      (gcking.gkin[index][0]) = particle->Px();
      (gcking.gkin[index][1]) = particle->Py();
      (gcking.gkin[index][2]) = particle->Pz();
      (gcking.gkin[index][3]) = particle->Energy(); EnergySum += particle->Energy();

      (gcking.gkin[index][4]) = float(g3id);
      //      particlePDG->Print();

      // Decay vertex
      (gckin3.gpos[index][0]) = gctrak.vect[0];
      (gckin3.gpos[index][1]) = gctrak.vect[1];
      (gckin3.gpos[index][2]) = gctrak.vect[2];

      // time of flight offset (mm)... (huh?)
      (gcking.tofd[index])    = 0.;

      // Set the flag to handle the particle in GSKING.   We currently skip particles
      // which we have decayed with the "continue" statements above.  This is because
      // the STAR stepping routine drops the iflgk flag before the call to gsking which
      // processes it.... 
      //
      // so in order to preserve the particle in the event record, we will need to add
      // a call to gsking in this routine.  (And test test test test...) 
      (gcking.iflgk[index])   = 0;

      // And increase stack counter
      index++;

      // Return the energy sum for validation
      return EnergySum;
}
//_____________________________________________________________________________________
AgUDecay::AgUDecay() : mDecayer( 0 ), 
		       mArray( new TClonesArray("TParticle",1000) ), 
		       mP(),
		       mDiscovery( kDecay ),
		       mNextG3id( 12345 ), // dynamic G3 id
                       mNonConservation(0) 
{

}
//_____________________________________________________________________________________
AgUDecay::~AgUDecay() {
  if ( mNonConservation ) {
    LOG_ERROR << "Energy was not conserved in " << mNonConservation << " decays" << endm; 
  }
}
