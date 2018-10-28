#include "HelloWorldMaker.h"
#include <iostream>

ClassImp(HelloWorldMaker)            // Macro for CINT compatibility


HelloWorldMaker::HelloWorldMaker( ) : StMaker( )
{ 
  // Class constructor ... zero pointers, zero public/private data members, etc. 
}


HelloWorldMaker::~HelloWorldMaker() 
{ 
  // Class destructor ... destroy and/or zero out pointers and public/private data members 
}


Int_t HelloWorldMaker::Init( )
{ 
  // Do once at the start the analysis, open files, create histograms, etc.
  mEventsProcessed = 0 ;
  return kStOK ;
}


Int_t HelloWorldMaker::Make( )
 
{
  // Do once at the start of every event
  cout << "***** Hello World *****" << endl ;
  mEventsProcessed++ ;
  return kStOK ;
}


Int_t HelloWorldMaker::Finish( )

{ // Do once at the end of the analysis, final updates, close files, etc.
  cout << "***** Number of events processed = " << mEventsProcessed << endl ;
  return kStOK ;
}









