#include <limits>

#include "TMath.h"
#include "St_base/StMessMgr.h"
#include "StPicoEvent/StPicoBTowHit.h"

ClassImp(StPicoBTowHit)

//_________________
StPicoBTowHit::StPicoBTowHit(): mId(0), mAdc(0), mE(0) {
  /* empty */
}

//_________________
StPicoBTowHit::StPicoBTowHit(int id, int adc, float e): StPicoBTowHit() {

  if (id  < 0 || adc < 0) return;
  mId   = (id  > std::numeric_limits<unsigned short>::max()) ? std::numeric_limits<unsigned short>::max() : (UShort_t)id;
  mAdc  = (adc > std::numeric_limits<unsigned short>::max()) ? std::numeric_limits<unsigned short>::max() : (UShort_t)adc;
  mE    = (e * 1000. > std::numeric_limits<short>::max()) ? std::numeric_limits<short>::max() : (Short_t)(TMath::Nint(e * 1000.));
}

//_________________
StPicoBTowHit::StPicoBTowHit(const StPicoBTowHit &hit) {
  mId = hit.mId;
  mAdc = hit.mAdc;
  mE = hit.mE;
}

//_________________
StPicoBTowHit::~StPicoBTowHit() {
  /* empty */
}

//_________________
void StPicoBTowHit::Print(const Char_t* option) const {
  LOG_INFO << " Id = " << id() << " Adc = " << adc() << " Energy = " << energy() << endm;
}
