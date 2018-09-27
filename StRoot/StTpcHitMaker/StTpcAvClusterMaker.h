#ifndef StTpcAvClusterMaker_H
#define StTpcAvClusterMaker_H
/***************************************************************************
 *
 * $Id: StTpcAvClusterMaker.h,v 1.2 2014/08/06 11:43:50 jeromel Exp $
 ***************************************************************************/
#include "TString.h"
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "THnSparse.h"
class StTpcAvClusterMaker : public StMaker {
 public:
  StTpcAvClusterMaker(const char *name="TpcAvCluster") : StMaker(name),  fAvLaser(0), NoRows(-1), NoInnerPadRows(-1) {}
  virtual ~StTpcAvClusterMaker() {}
  Int_t   InitRun(Int_t runnumber);
  Int_t   Make();
  THnSparseF *CompressTHn(THnSparseF *hist, Double_t compress = 1e3);
  virtual Int_t        Finish();
 private:
  THnSparseF **fAvLaser;
  Int_t    NoRows;
  Int_t    NoInnerPadRows;
 public:
  virtual const char *GetCVS() const    {
    static const char cvs[]="Tag $Name: SL18c $Id: built " __DATE__ " " __TIME__ ; return cvs;
  }
  ClassDef(StTpcAvClusterMaker, 1)    //StTpcAvClusterMaker - class to fille the StEvewnt from DAQ reader
};

#endif
