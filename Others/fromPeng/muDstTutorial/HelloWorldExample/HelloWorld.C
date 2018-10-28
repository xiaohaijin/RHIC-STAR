void HelloWorld( ) 
{
 
  // Load libraries
  gROOT   -> Macro("loadMuDst.C") ;
  gSystem -> Load("HelloWorldMaker.so") ;

  // List of member links in the chain
  StChain* chain = new StChain ;
  HelloWorldMaker* Hello = new HelloWorldMaker( ) ;

  Int_t nEvents = 10 ;

  // Loop over the links in the chain
  chain -> Init() ;
  chain -> EventLoop(1,nEvents) ;
  chain -> Finish() ;

  // Cleanup
  delete chain ;

}

