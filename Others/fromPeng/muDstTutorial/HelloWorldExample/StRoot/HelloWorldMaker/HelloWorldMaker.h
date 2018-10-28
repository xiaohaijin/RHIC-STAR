#ifndef HelloWorldMaker_def
#define HelloWorldMaker_def

#include "StMaker.h"

class HelloWorldMaker : public StMaker
{
  
 private:
  
  ULong_t mEventsProcessed ;    

 protected:

 public:

  HelloWorldMaker( ) ;              //  Constructor
  virtual ~HelloWorldMaker( ) ;     //  Destructor

  Int_t Init     ( ) ;              //  Initiliaze the analysis tools ... done once
  Int_t Make     ( ) ;              //  The main analysis that is done on each event
  Int_t Finish   ( ) ;              //  Finish the analysis, close files, and clean up.

  ClassDef(HelloWorldMaker,1)       //  Macro for CINT compatability
    
};

#endif















