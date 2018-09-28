#include "SimpleAnalysisMaker.h"
#include <iostream>

#include "COMMON/StMuDstMaker.h"
#include "COMMON/StMuEvent.h"
#include "COMMON/StMuTrack.h"

#include "TFile.h"
#include "TH1F.h"
#include "TObjArray.h"

ClassImp(SimpleAnalysisMaker)

    SimpleAnalysisMaker::SimpleAnalysisMaker(StMuDstMaker *maker,
                                             const char *name)
    : StMaker("name") {
  mMuDstMaker = maker;

  hVertexZ = NULL;
  hSelectNRefMult = NULL;
  hpt = NULL;
  histogram_output = NULL;

  mEventsProcessed = 0;

  mHistogramOutputFileName = "";
}
SimpleAnalysisMaker::~SimpleAnalysisMaker() {}

Int_t SimpleAnalysisMaker::Init() {
  if (mHistogramOutputFileName != "") {
    histogram_output = new TFile(mHistogramOutputFileName, "recreate");
  } else {
    cout << "Error: Please specify the histogtam output file!!!!!!!!!!!!!!"
         << endl;
    exit(-1);
  }

  hVertexZ = new TH1F("VertexZ", "Event Vertex Z Position", 20000, -100, 100);
  hSelectNRefMult =
      new TH1F("SelectRefMult", "Reference of Multiplicity of selected events",
               1000, 0.0, 1000.0);
  hpt = new TH1F("pt", "the pt of tracks", 200, 0.0, 20.0);
  magn1 = new TH1F("magn1", "mag1", 100, -10.0, 10.0);
  magn2 = new TH1F("magn2", "mag2", 200, -10.0, 10.0);
  test = new TH1F("test", "test", 200, -10.0, 10.0);
  hitsposs = new TH1F("hitsposs", "hitsposs", 100, -100.0, 100.0);
  beta1 = new TH1F("beta1", "beta1", 100, -5.0, 5.0);
  beta2 = new TH1F("beta2", "beta2", 100, -5.0, 5.0);
  return kStOK;
}

Int_t SimpleAnalysisMaker::Make() {
  StMuEvent *muEvent = mMuDstMaker->muDst()->event();

  if (!muEvent) return kStOK;

  //	if(fabs(muEvent->primaryVertexPosition().x())< 1e-5 &&
  // fabs(muEvent->primaryVertexPosition().y())< 1e-5 &&
  // fabs(muEvent->primaryVertexPosition().z())< 1e-5)return kStOK;
  //	if(pow(muEvent->primaryVertexPosition().x(),2)+pow(muEvent->primaryVertexPosition().y(),2)>4)return
  // kStOK;

  //	if(
  //			(!muEvent->triggerIdCollection().nominal().isTrigger(450008))&&
  //			(!muEvent->triggerIdCollection().nominal().isTrigger(450009))&&
  //			(!muEvent->triggerIdCollection().nominal().isTrigger(450010))&&
  //			(!muEvent->triggerIdCollection().nominal().isTrigger(450011))&&
  //			(!muEvent->triggerIdCollection().nominal().isTrigger(450014))&&
  //			(!muEvent->triggerIdCollection().nominal().isTrigger(450015))&&
  //			(!muEvent->triggerIdCollection().nominal().isTrigger(450025))&&
  //			(!muEvent->triggerIdCollection().nominal().isTrigger(450050))&&
  //			(!muEvent->triggerIdCollection().nominal().isTrigger(450060))&&
  //			(!muEvent->triggerIdCollection().nominal().isTrigger(450103))
  //	  ) return kStOK;
  //
  hVertexZ->Fill(muEvent->primaryVertexPosition().z());

  double Magn1 = muEvent->runInfo().magneticField();
  double Magn2 = muEvent->magneticField();

  magn1->Fill(Magn1);
  magn2->Fill(Magn2);
  test->Fill(Magn1 - Magn2);

  //	if(fabs(muEvent->primaryVertexPosition().z())>50.0)return kStOK;

  hSelectNRefMult->Fill(muEvent->refMult());

  TObjArray *tracks = mMuDstMaker->muDst()->primaryTracks();
  TObjArrayIter GetTracks(tracks);

  StMuTrack *track = 0;
  /// deal with each track in an event.
  while ((track = (StMuTrack *)GetTracks.Next())) {
    const StMuBTofPidTraits &tofpid = track->btofPidTraits();

    beta1->Fill(tofpid.beta());
    double tofbeta = -999.0;
    tofbeta = tofpid.beta();
    if (tofbeta > 0) beta2->Fill(tofpid.beta());

    hpt->Fill(track->pt());
    hitsposs->Fill(track->nHitsPoss());
  }
  mEventsProcessed++;
  return kStOK;
}
Int_t SimpleAnalysisMaker::Finish() {
  histogram_output->Write();
  cout << "Total Events Processed in DstMaker  " << mEventsProcessed << endl;
  return kStOK;
}
