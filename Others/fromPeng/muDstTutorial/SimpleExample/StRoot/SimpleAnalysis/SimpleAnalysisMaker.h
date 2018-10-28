// Based on the MuDST tools written by Frank Laue.
// Based on the DST Tutorial by Dan Magestro on the STAR Computing/Tutorials page.
// Updated 9/4/2006 by Jim Thomas to include the latest DST format, and Scheduler techniques.

#ifndef SimpleAnalysisMaker_def
#define SimpleAnalysisMaker_def

#include "StMaker.h"
#include "TString.h"

class StMuDstMaker ;
class TFile        ;
class TH1F         ;

#define MaxNumberOfTH1F     10

class SimpleAnalysisMaker : public StMaker
{
  
 private:

  StMuDstMaker* mMuDstMaker ;                      //  Make MuDst pointer available to member functions

  TH1F*         histogram[MaxNumberOfTH1F] ;       //  1D Histograms
  TFile*        histogram_output ;                 //  Histograms outputfile pointer

  UInt_t        mEventsProcessed ;                 //  Number of Events read and processed
  TString       mHistogramOutputFileName ;         //  Name of the histogram output file 


 protected:


 public:

  SimpleAnalysisMaker(StMuDstMaker* maker) ;       //  Constructor
  virtual          ~SimpleAnalysisMaker( ) ;       //  Destructor

  Int_t Init    ( ) ;                              //  Initiliaze the analysis tools ... done once
  Int_t Make    ( ) ;                              //  The main analysis that is done on each event
  Int_t Finish  ( ) ;                              //  Finish the analysis, close files, and clean up.

  void SetOutputFileName(TString name) {mHistogramOutputFileName = name;} // Make name available to member functions
  
  ClassDef(SimpleAnalysisMaker,1)                  //  Macro for CINT compatability
    
};

#endif















