#ifndef StV0Maker_def
#define StV0Maker_def

#include "StMaker.h"
#include "StPhysicalHelixD.hh"
#include "StThreeVectorD.hh"
#include "TString.h"

class StMuDstMaker;
class StMuTrack;
class TFile;
class TTree;
class TH1F;
class TH2F;

#include "StV0Dst.h"
#include "StV0Type.h"

/// \brief StV0Maker用来重构V0粒子
///
/// 该code来自于陈金辉老师，不过从code内的注释来看，
/// 主要的工作来自Longhui和朱相雷老师。
class StV0Maker : public StMaker {
 private:
  /// Make MuDst pointer available to member functions
  StMuDstMaker* mMuDstMaker;

  /// master switch for the analysis, v0 type
  StV0Type mV0Type;

  /// switch for rotating the coordinates and momenta (in transverse direction)
  /// of one daughter. for background estimation.
  bool mRotate;
  bool mSameSignPlus;
  bool mSameSignMinus;
  bool mDcaAlgoLong;

  /// section of v0type related constants
  Double_t mMassV0;
  Double_t mMass1;
  Double_t mMass2;
  Int_t mChargeV0;
  Int_t mCharge1;
  Int_t mCharge2;

  /// section of parameters for v0 analysis
  double cutAbsVertexZLeEq;
  int cutNHitsGr;
  int cutNHitsDedxGr;
  double cutPtGrEq;

  double cutAbsNSigma2Le;
  double cutDca1GrEq;
  double cutDca2GrEq;

  double cutDca1to2LeEq;
  double cutV0MassWidthLeEq;
  double cutDauPtArmLeEq;
  double cutAbsDausPtShoulderDiffLeEq;
  double cutDau1DecAngGr;
  double cutDau2DecAngGr;
  double cutV0rdotpGr;
  double cutDcaV0Le;
  double cutV0DecLenGrEq;
  double cutDau1Dau2Ang3DLe;
  double cutDau1Dau2DipAngDiffLe;

  double RotDegree;

  /// histograms (mostly for QA purpose)
  TH1F* hNPrimVertex;
  TH1F* hVertexZ;
  TH1F* hVertexZdiff;
  TH1F* hNRefMult;
  TH1F* hSelectNRefMult;
  TH1F* hMagneticField;

  TH2F* hnSigmaProton;
  TH2F* hnSigmaPion;
  TH2F* hnSigmaKaon;
  TH2F* hdEdxP;
  TH2F* hDcaP;
  TH2F* hMassP;
  TH2F* hInvBetaP;
  TH2F* hdau1dEdxP;
  TH2F* hdau2dEdxP;
  TH2F* hdau1ZP;
  TH2F* hdau2ZP;
  TH2F* hdau1MassP;
  TH2F* hdau2MassP;
  TH2F* hdau1InvBetaP;
  TH2F* hdau2InvBetaP;
  TH2F* hdau1DiffInvBetaP;
  TH2F* hdau2DiffInvBetaP;

  TH1F* hInvMass;

  /// dEdx information of tracks.
  double dedx_dau1_th[11901];
  double dedx_dau2_th[11901];

  /// files related
  /// Name of the histogram output file
  TString mHistogramOutputFileName;
  /// Name of the v0 tree output file
  TString mV0TreeOutputFileName;
  /// Histograms outputfile pointer
  TFile* histogram_output;
  /// V0 Tree outputfile pointer
  TFile* v0tree_output;
  /// V0 Tree outputfile pointer
  TTree* mV0Tree;
  /// V0 event (picoDst), to fill the tree
  StV0Dst mV0Dst;

  /// statistic information
  ///  Number of Events read and processed
  UInt_t mEventsProcessed;

  /// some diagnosing variables
  Double_t mTestVZ;
  UInt_t mTestNTrack;
  bool mPassEventCut;

  /// temporary variables (put here for performance consideration)
  vector<StMuTrack*> mDauVec1;
  vector<StMuTrack*> mDauVec2;
  vector<double> mDauDcaVec1;
  vector<double> mDauDcaVec2;
  vector<double> mDauZVec1;
  vector<double> mDauZVec2;
  vector<double> mDauMassVec1;
  vector<double> mDauMassVec2;
  vector<double> mDauBetaVec1;
  vector<double> mDauBetaVec2;
  vector<double> mDauDiffInvBetaVec1;
  vector<double> mDauDiffInvBetaVec2;

  /// private member functions
  /// initialize parameters
  void initParam();
  /// initialize constants for v0type
  void initConst();
  /// book histograms
  void initHisto();
  /// book tree
  void initTree();

 protected:
  /// I do not expect some class inherits this maker!

 public:
  /// Constructor
  StV0Maker(StMuDstMaker* maker, const char* name);
  /// Destructor
  virtual ~StV0Maker();

  /// Initiliaze the analysis tools ... done once
  Int_t Init();
  /// The main analysis that is done on each event
  Int_t Make();
  /// Finish the analysis, close files, and clean up.
  Int_t Finish();

  /// set the name of output file for histograms
  void setHistoFileName(TString name) { mHistogramOutputFileName = name; }
  /// set the name of output file for StV0Dst
  void setV0TreeFileName(TString name) { mV0TreeOutputFileName = name; }
  /// set the v0 type: kLambda,kAntiLambda,kKs
  void setV0Type(StV0Type v0type) { mV0Type = v0type; }
  /// set rotation option
  void setRotate(bool brot) { mRotate = brot; }
  /// set same sign plus
  void setSameSignPlus(bool brot) { mSameSignPlus = brot; }
  /// set same sign minus
  void setSameSignMinus(bool brot) { mSameSignMinus = brot; }

  const StV0Dst& getV0Dst() const { return mV0Dst; }
  bool passEventCut() const { return mPassEventCut; }
  /// Rotate dau1 or dau2 tracks by certain degrees in the azimuthal plane with
  /// respect to the primary vertex.
  StPhysicalHelixD RotHelix(StPhysicalHelixD Helix, StThreeVectorD Pv,
                            double degree, double magn, int charge);

  int index_for_p(double);

  /// Macro for CINT compatability
  ClassDef(StV0Maker, 1)
};

#endif
