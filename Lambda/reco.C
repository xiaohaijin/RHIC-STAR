////////////////////////////////////////////////////////////////////////
/// COPYRIGHT NOTICE
/// COPYRIGHT (c) 2018, 金小海
/// All rights reserved.
///
/// @file    This file is part of the Lambda project.
/// @version 1.0
/// @author  jinxiaohai <xiaohaijin@outlook.com>
/// @date    2018 06 07    13:25:30
/// @brief   topologically reconstruct Lambda and Anti-Lambda
///
/// 修订说明:最初版本
////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/// was added by xiaohai (begin)
///////////////////////////////////////////////////////////////////////
#define XIAOHAI

#ifdef XIAOHAI
#include "COMMON/StMuDstMaker.h"
#include "StChain.h"
#include "StMaker.h"
#include "StRoot/StV0Maker/StV0Maker.h"

#include "TChain.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#endif
///////////////////////////////////////////////////////////////////////
/// was added by xiaohai (end)
///////////////////////////////////////////////////////////////////////

#include "StRoot/StV0Maker/StV0Type.h"

/// give the default parameters ????
void reco(TString InputFileList, Int_t nFiles = 1000, Int_t nEvents = 0,
          TString OutputDir = "output/", TString JobIdName = "testrun");

void reco(TString InputFileList, Int_t nFiles, Int_t nEvents, TString OutputDir,
          TString JobIdName) {
  /// Load files and libraries
  gROOT->Macro("loadMuDst.C");
  gSystem->Load("StV0Maker");

  /// Main base class to control chains for the different STAR "chains"
  /// This class:
  ///     * Initialises the run default parameters
  ///     * Provide API to Set/Get run parameters
  ///     * Creates the support lists (TClonesArrays) for the Event structure
  ///     * Creates the physics objects makers
  /// default construct function is StChain(const char *name="bfcChain", const
  /// Bool_t UseOwnHeader=kFALSE)
  StChain* chain = new StChain;
  /// Class to create and read STAR's common micro dst (StMulDst)
  /// This class is a true maker in the STAR sense. It inherits from "StMaker"
  /// and implements the functions "int Init()", "void Clear()", "int Make()",
  /// and "int Finish()" in order to run as part of an "StChain". Please refer
  /// to the STAR Computing Web pages in case you do not know what "StMaker" and
  /// "StChain" mean.
  /// constructor is: StMuDstMaker(int mode,
  ///                              int nameMode,
  ///                              const char *dirName="./",
  ///                              const char *fileName="",
  ///                              const char *filter=".",
  ///                              int maxfiles=10,
  ///                              const char *name="MuDst")
  StMuDstMaker* muDstMaker =
      new StMuDstMaker(0, 0, "", InputFileList, "MuDst", nFiles);

  /// Turn off everything but Primary tracks in order to speed up the analysis
  /// and eliminate IO
  /// muDstMaker -> SetStatus("MuEvent",1); Turn on the Event data (esp. Event
  /// number)<---
  /// muDstMaker -> SetStatus("GlobalTracks",1); Turn on the
  /// global track data<---
  /// muDstMaker -> SetStatus("BTofHeader",1); Turn on the global track
  /// data<---
  /// Turn off all branches
  muDstMaker->SetStatus("*", 1);

  /// Turn off Debug information
  muDstMaker->SetDebug(0);

  /////////////////////////////////////////////////////////////////////
  /// use StV0Maker to get the v0 candidates (begin)
  /////////////////////////////////////////////////////////////////////
  /// Miscellaneous things we need before starting the chain
  /// Miscellaneous : 混杂的，各式各样的，五花八门的。
  StV0Maker* rawsig = new StV0Maker(muDstMaker, "v0makerfirst_signal");
  /// Name the output file for histograms
  rawsig->setHistoFileName(OutputDir + JobIdName + ".la.histo.root");
  /// V0 candidate tree file for further cuts.
  rawsig->setV0TreeFileName(OutputDir + JobIdName + ".la.picodst.root");
  /// set V0 type. once a time! do not try to mess
  rawsig->setV0Type(kLambda);
  /// things up for the time being.
  rawsig->setRotate(false);
  rawsig->SetDebug(0);

  ///  use StV0Maker to get the anti-v0 candidates
  StV0Maker* rawantisig = new StV0Maker(muDstMaker, "v0makerfirst_antisignal");
  /// Name the output file for histograms
  rawantisig->setHistoFileName(OutputDir + JobIdName + ".ala.histo.root");
  /// V0 candidate tree file for further cuts.
  rawantisig->setV0TreeFileName(OutputDir + JobIdName + ".ala.picodst.root");
  /// set V0 type. once a time! do not try to mess
  rawantisig->setV0Type(kAntiLambda);
  /// things up for the time being.
  rawantisig->setRotate(false);
  rawantisig->SetDebug(0);

  /// use Rotate method to get the background information.
  StV0Maker* rawbg = new StV0Maker(muDstMaker, "v0makerfirst_background");
  rawbg->setHistoFileName(OutputDir + JobIdName + ".labg.histo.root");
  rawbg->setV0TreeFileName(OutputDir + JobIdName + ".labg.picodst.root");
  rawbg->setV0Type(kLambda);
  rawbg->setRotate(true);
  rawbg->SetDebug(0);

  ///  use Rotate method to get the anti-background information.
  StV0Maker* rawantibg =
      new StV0Maker(muDstMaker, "v0makerfirst_antibackground");
  TString Name = JobIdName;
  Name.Append(".alabg.histo.root");
  rawantibg->setHistoFileName(OutputDir + JobIdName + ".alabg.histo.root");
  rawantibg->setV0TreeFileName(OutputDir + JobIdName + ".alabg.picodst.root");
  rawantibg->setV0Type(kAntiLambda);
  rawantibg->setRotate(true);
  rawantibg->SetDebug(0);
  /////////////////////////////////////////////////////////////////////
  /// use StV0Maker to get the v0 candidates (end)
  /////////////////////////////////////////////////////////////////////

  /// Take all events in nFiles if nEvents = 0
  if (nEvents == 0) {
    nEvents = 10000000;
  }

  /// Loop over the links in the chain
  /// Init() function override the Init() in StMaker.
  Int_t iInit = chain->Init();
  if (iInit) {
    chain->Fatal(iInit, "on init");
  }
  Int_t totalE = static_cast<int>(muDstMaker->chain()->GetEntries());

  /// will output lots of useless debugging
  /// chain -> EventLoop(1,nEvents);
  /// info.
  Int_t istat = 0, i = 1;
  cout << " Total entries = " << totalE << endl;
  if (nEvents > totalE) {
    nEvents = totalE;
  }

  while (i <= nEvents && istat != 2) {
    if (i % 1000 == 0) {
      cout << endl << "== Event " << i << " start ==" << endl;
    }

    chain->Clear();
    istat = chain->Make(i);
    /// cout << endl << "== Event " << i << " finish =="<< endl;
    if (istat == 2) {
      cout << "Last  event processed. Status = " << istat << endl;
    }
    if (istat == 3) {
      cout << "Error event processed. Status = " << istat << endl;
    }

    i++;
  }

  if (nEvents > 1) {
    chain->Finish();
  }

  // Cleanup
  delete chain;
}
