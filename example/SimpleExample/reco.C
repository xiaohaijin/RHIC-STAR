///////////////////////////////////////////////////////////////////////
/// was added by xiaohai (begin)
///////////////////////////////////////////////////////////////////////
#include "COMMON/StMuDstMaker.h"
#include "StChain.h"
#include "StRoot/SimpleAnalysisMaker/SimpleAnalysisMaker.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
///////////////////////////////////////////////////////////////////////
/// was added by xiaohai (end)
///////////////////////////////////////////////////////////////////////

void reco(TString InputFileList, Int_t nFiles = 1000, Int_t nEvents = 0,
          TString OutputDir = "output/", TString JobIdName = "testrun");

void reco(TString InputFileList, Int_t nFiles, Int_t nEvents, TString OutputDir,
          TString JobIdName) {
  /// Macro()
  /// Execute a macro in the interpreter.
  /// Equivalent to the command line command ".x filename". if the filename has
  /// "+" or "++" appended the macro will be compiled by ACLiC. The filename
  /// must have the format [path/]macro.C[+|++]. The possible error codes are
  /// defined by TInterpreter::EErrorCode. If padUpdate is true (default) update
  /// the current pad. Returns the macro return value.
  /// 其中loadMuDst.C文件里面是全是gSystem->Load("XXXX")加载动态库的语句，为了模块化
  /// 的简单控制才这样做的。
  gROOT->Macro("loadMuDst.C");
  /// Load()
  /// Load a shared library
  /// Returns 0 on successful loading, 1 in case lib was already loaded, -1 in
  /// case lib does not exist or in case of error and -2 in case of version
  /// mismatch. When entry is specified the loaded lib is searched for this
  /// entry point (return -1 when entry does not exist 0, 0 otherwise). When the
  /// system flag is kTRUE, the library is considered a permanent system library
  /// that should not be unloaded during the course of the session.
  gSystem->Load("SimpleAnalysisMaker.so");

  /// Main base class to control chains for the different STAR "chains"
  /// This class:
  ///     * initialises the run default parameters
  ///     * Provide API to Set/Get run parameters
  ///     * Creates the support lists (TClonesArray) for the Event structure
  ///     * Creates the physics objects makers
  StChain* chain = new StChain;

  ///
  /// \brief muDstMaker
  /// Constructor
  /// StMuDstMaker(int mode, int nameMode, cosnt char *dirName="./", const char
  /// *filter=".", int maxfiles=10, const char *name="MuDst")
  StMuDstMaker* muDstMaker =
      new StMuDstMaker(0, 0, "", InputFileList, "MuDst", nFiles);
  /// Turn off everything but Primary tracks in order to speed up the analysis
  /// and eliminate IO
  /// muDstMaker -> SetStatus("*",0);// Turn off all branches
  /// muDstMaker -> SetStatus("MuEvent",1); // Turn on the Event data (esp.
  /// Event number) muDstMaker -> SetStatus("PrimaryTracks",1);// Turn on the
  /// global track data
  /// Turn off Debug information
  muDstMaker->SetDebug(0);

  /// use StV0Maker to get the v0 candidates
  SimpleAnalysisMaker* test = new SimpleAnalysisMaker(muDstMaker, "testMaker");

  /// Miscellaneous things we need before starting the chain
  TString Name = JobIdName;
  Name.Append(".test.histo.root");
  /// Name the output file for histograms
  test->SetOutputFileName(OutputDir + Name);
  test->SetDebug(0);

  if (nEvents == 0) {
    /// Take all events in nFiles if nEvents = 0
    nEvents = 10000000;
  }

  /// Loop over the links in the chain
  Int_t iInit = chain->Init();
  if (iInit) {
    chain->Fatal(iInit, "on init");
  }
  Int_t totalE = muDstMaker->chain()->GetEntries();

  /// chain -> EventLoop(1,nEvents);//will output lots of useless debugging
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

  /// Cleanup
  delete chain;
}
