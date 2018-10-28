// Based on the MuDST tools written by Frank Laue.
// Based on the DST Tutorial by Dan Magestro on the STAR Computing/Tutorials page.
// Updated 9/4/2006 by Jim Thomas to include the latest DST format, and Scheduler techniques.

#include "BiggerAnalysisMaker.h"

#include <iostream>

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"

#include "TH1.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TList.h"

#define NumberOfTH1F      2                     // Number of Histograms

ClassImp(BiggerAnalysisMaker)                   // Macro for CINT compatibility


BiggerAnalysisMaker::BiggerAnalysisMaker( StMuDstMaker* maker ) : StMaker("BiggerAnalysisMaker")

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


BiggerAnalysisMaker::~BiggerAnalysisMaker() 

{ // Destroy and/or zero out all public/private data members here.

}


Int_t BiggerAnalysisMaker::Init( )

{ // Do once at the start of the analysis

  // Create Histogram output file

  histogram_output = new TFile( mHistogramOutputFileName, "recreate" ) ;  // Name was set in "analysis".C macro

  // Create Histograms

  const Int_t    nbins    =  100   ;
 
  histogram[0]  = new TH1F( "Vertex", "Event Vertex Z Position", nbins, -25.0, 25.0 ) ; 
  histogram[1]  = new TH1F( "Pt", "Transverse Momentum for all particles", nbins, 0.0, 10.0 ) ;
  histogram[2]  = new TH1F( "highPt", "High Transverse Momentum particles", nbins, 0.0, 10.0 ) ;

  return kStOK ;

}


Int_t BiggerAnalysisMaker::Make( )
 
{ // Do each event

  // Get 'event' data 
  
  StMuEvent* muEvent      =  mMuDstMaker->muDst()->event() ;

  // Cut on the number of vertices in the event.  On old tapes, no-vertex gets reported as VtxPosition=(0,0,0).
  // On newer tapes, we have a value for the number of vertices ... but this value is reported equal to zero on
  // old tapes that have a good vertex ... so we do a test on both methods of selecting events with only one vertex.

  if ( (muEvent->primaryVertexPosition().x() == 0) && (muEvent->primaryVertexPosition().y() == 0) 
       && (muEvent->primaryVertexPosition().z() == 0) ) return kStOK ;  // Skip events that do not have a primary vertex

  Int_t NumberOfPrimaryVertices = mMuDstMaker->muDst()->numberOfPrimaryVertices() ;
  if (NumberOfPrimaryVertices >= 2 ) return kStOK ;              // Skip events with more than one primary vertex

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

  // Option - alternatively you can store the tracks in a TList and do your physics analysis on the TList.
  // This is a good way to store several sub-sets of data from the same event so you can do pairwise analyses. 
  // Or it might be a way to store several events worth of data to do mixed event background studies.

  TList HighPtTracks ;                                           // Prepare to store a subset of events in a TList

  GetTracks.Reset()  ;                                           // Reset the track iterator so we can do it again
  while ( ( track = (StMuTrack*)GetTracks.Next() ) )
    {
      if ( track->pt() >= 2.0 ) HighPtTracks.Add(track) ;        // Very simple subset of tracks
    }

  TListIter GetHighPtTracks(&HighPtTracks) ;                     // Create the iterator for the events in the TList

  StMuTrack* listedtrack ;
  while ( ( listedtrack = (StMuTrack*)GetHighPtTracks.Next() ))
    {
      histogram[2] -> Fill ( listedtrack->pt() ) ;               // Very simple 'analysis' of the TList'ed data
    }

  mEventsProcessed++ ;
  return kStOK ;
  
}


Int_t BiggerAnalysisMaker::Finish( )

{ // Do once at the end of the analysis

  // Write histograms to disk, output miscellaneous other information

  histogram_output -> Write() ;   // Write all histograms to disk 

  cout << "Total Events Processed in DstMaker " << mEventsProcessed << endl ;

  return kStOK ;

}









