// Based on the MuDST tools written by Frank Laue.
// Based on the DST Tutorial by Dan Magestro on the STAR Computing/Tutorials page.
// Updated 9/4/2006 by Jim Thomas to include the latest DST format, and Scheduler techniques.

#include "SimpleAnalysisMaker.h"

#include <iostream>

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"

#include "TH1.h"
#include "TFile.h"
#include "TObjArray.h"

#define NumberOfTH1F      2                     // Number of Histograms

ClassImp(SimpleAnalysisMaker)                   // Macro for CINT compatibility


SimpleAnalysisMaker::SimpleAnalysisMaker( StMuDstMaker* maker ) : StMaker("SimpleAnalysisMaker")

{ // Initialize and/or zero all public/private data members here.

  for ( Int_t i = 0 ; i < NumberOfTH1F ; i++ )  // Zero the histogram pointers
    {
      histogram[i] = NULL ;
    }

  mMuDstMaker      = maker ;                    // Pass MuDst pointer to DstAnlysisMaker Class member functions
  histogram_output = NULL  ;                    // Zero the Pointer to histogram output file
  mEventsProcessed = 0     ;                    // Zero the Number of Events processed by the maker 
  mHistogramOutputFileName = "" ;               // Histogram Output File Name will be set inside the "analysis".C macro

}


SimpleAnalysisMaker::~SimpleAnalysisMaker() 

{ // Destroy and/or zero out all public/private data members here.

}


Int_t SimpleAnalysisMaker::Init( )

{ // Do once at the start of the analysis

  // Create Histogram output file

  histogram_output = new TFile( mHistogramOutputFileName, "recreate" ) ;  // Name was set previously in calling macro 

  // Create Histograms

  const Int_t    nbins    =  100   ;
 
  histogram[0]  = new TH1F( "Vertex", "Event Vertex Z Position", nbins, -25.0, 25.0 ) ; 
  histogram[1]  = new TH1F( "Pt", "Transverse Momentum for all particles", nbins, 0.0, 10.0 ) ;

  return kStOK ; 

}


Int_t SimpleAnalysisMaker::Make( )
 
{ // Do each event

  // Get 'event' data 
  
  StMuEvent* muEvent      =  mMuDstMaker->muDst()->event() ;

  // Cut on the number of vertices in the event.  On old tapes, no-vertex gets reported as VtxPosition=(0,0,0).

  if ( (muEvent->primaryVertexPosition().x() == 0) && (muEvent->primaryVertexPosition().y() == 0) 
       && (muEvent->primaryVertexPosition().z() == 0) ) return kStOK ;  // Skip events that do not have a primary vertex

  // Do 'event' analysis based on event data 

  histogram[0] -> Fill( muEvent->primaryVertexPosition().z() ) ; // Make histogram of the vertex Z distribution

  // Get 'track' data, make cuts on tracks, do physics analysis, histogram results.
 
  TObjArray* tracks = mMuDstMaker->muDst()->primaryTracks() ;    // Create a TObject array containing the primary tracks
  TObjArrayIter GetTracks(tracks) ;                              // Create an iterator to step through the tracks

  StMuTrack* track ;                                             // Pointer to a track
  while ( ( track = (StMuTrack*)GetTracks.Next() ) )             // Main loop for Iterating over tracks
    {
      histogram[1] -> Fill( track->pt() ) ;
    }

  mEventsProcessed++ ;
  return kStOK ;
  
}


Int_t SimpleAnalysisMaker::Finish( )

{ // Do once at the end the analysis

  // Write histograms to disk, output miscellaneous other information

  histogram_output -> Write() ;   // Write all histograms to disk 

  cout << "Total Events Processed in DstMaker " << mEventsProcessed << endl ;

  return kStOk ;  

}









