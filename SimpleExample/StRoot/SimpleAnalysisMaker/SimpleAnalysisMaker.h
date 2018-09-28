#ifndef SimpleAnalysisMaker_def
#define SimpleAnalysisMaker_def

#include <StMaker.h>
#include <TString.h>

class StMuDstMaker;
class TFile;
class TH1F;

class SimpleAnalysisMaker : public StMaker {
 private:
  ULong_t mEventsProcessed;

  StMuDstMaker *mMuDstMaker;

  TH1F *hVertexZ;
  TH1F *hpt;
  TH1F *hSelectNRefMult;
  TFile *histogram_output;
  TH1F *magn1;
  TH1F *magn2;
  TH1F *test;
  TH1F *hitsposs;
  TH1F *beta1;
  TH1F *beta2;

  TString mHistogramOutputFileName;

 protected:
 public:
  SimpleAnalysisMaker(StMuDstMaker *maker, const char *name);
  virtual ~SimpleAnalysisMaker();

  Int_t Init();
  Int_t Make();
  Int_t Finish();

  void SetOutputFileName(TString name) { mHistogramOutputFileName = name; }

  ClassDef(SimpleAnalysisMaker, 1)
};
#endif
