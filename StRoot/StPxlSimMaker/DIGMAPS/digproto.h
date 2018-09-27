#ifndef MAPS_DIGPROTO_H
#define MAPS_DIGPROTO_H

#include <TNamed.h>
#include <TList.h>
#include <TGraph.h>
#include "Riostream.h"
#include "vector"

// ROOT classes
#include "TString.h"
#include "TObject.h"
#include "TVector.h"
#include "TFile.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TH1.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
using namespace std;

class DIGInitialize;
class DIGPlane;
class DIGADC;
class DIGBeam;
class DIGParticle;
class DIGAction;
class DIGEvent;
class DIGMAPS;
/*
This class is a dummy class to be used as a framework if we want to create a new class for DIGMAPS.
*/
//==============================================================================
class DIGProto : public TObject {
 public:
  DIGProto();
  DIGProto(Float_t myvar);
  DIGProto(DIGProto& adigproto);
  virtual ~DIGProto();
  void    Clear(const Option_t * /*option*/ = "");
  void PrintInfo();

  Float_t GetMyvar(){return fMyvar;}
  
  void SetMyvar(Float_t Myvar);


 protected:

 Float_t fMyvar;

  ClassDef(DIGProto,1);
};



//==============================================================================

#endif
