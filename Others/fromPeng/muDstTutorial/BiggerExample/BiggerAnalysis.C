void BiggerAnalysis(Int_t nEvents, Int_t nFiles, TString InputFileList, TString OutputDir, TString JobIdName ) 
{
 
  // Load libraries
  gROOT   -> Macro("loadMuDst.C");
  gSystem -> Load("BiggerAnalysis.so") ;

  // List of member links in the chain
  StChain*                    chain  =  new StChain ;
  StMuDstMaker*          muDstMaker  =  new StMuDstMaker(0,0,"",InputFileList,"MuDst",nFiles) ;
  BiggerAnalysisMaker* AnalysisCode  =  new BiggerAnalysisMaker(muDstMaker) ;

  // Turn off everything but Primary tracks in order to speed up the analysis and eliminate IO
  muDstMaker -> SetStatus("*",0) ;                // Turn off all branches
  muDstMaker -> SetStatus("MuEvent",1) ;          // Turn on the Event data (esp. Event number)
  muDstMaker -> SetStatus("PrimaryTracks",1) ;    // Turn on the primary track data

  // Miscellaneous things we need before starting the chain
  TString Name = JobIdName ;
  Name.Append(".histograms.root") ;
  AnalysisCode -> SetOutputFileName(Name) ;       // Name the output file for histograms
  if ( nEvents == 0 )  nEvents = 10000000 ;       // Take all events in nFiles if nEvents = 0

  // Loop over the links in the chain
  chain -> Init() ;
  chain -> EventLoop(1,nEvents) ;
  chain -> Finish() ;

  // Cleanup
  delete chain ;

}

