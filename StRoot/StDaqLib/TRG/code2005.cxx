// This is the analog of duplicated.code for the version of trgStructures.h
// which I (Herb) label "2005".

#define PREPOST 11 // This is also defined in TRG_Reader.cxx.

#include <assert.h>
#include "trgStructures2005.h"
#include "TRG_Reader.hh"

typedef struct {
  EvtDescData2005 EvtDesc;  /* L1 Event Descriptor Data */
  TrgSumData2005  TrgSum;   /* summary data */
  RawTrgDet2005   RAW[PREPOST];      
} MarilynMonroe;
MarilynMonroe *gs2005;
using namespace OLDEVP;

int Bank_TRGD::HerbSwap2005(char *ptr) {
  int numToSwap,returnValue,i;

  gs2005=(MarilynMonroe*)ptr;

  assert(header.ByteOrder==0x01020304||header.ByteOrder==0x04030201);
  if(header.ByteOrder==0x04030201) return 0;
  returnValue=header.swap();
  assert(header.ByteOrder==0x04030201);

  swapHerb2bytes(&(gs2005->EvtDesc.TCUdataBytes),1);
  swapHerb4bytes(&(gs2005->EvtDesc.bunchXing_hi),1);
  swapHerb4bytes(&(gs2005->EvtDesc.bunchXing_lo),1);
  swapHerb2bytes(&(gs2005->EvtDesc.actionWdDetectorBitMask),1);
  swapHerb2bytes(&(gs2005->EvtDesc.TrgToken),1);
  swapHerb2bytes(&(gs2005->EvtDesc.addBits),1);
  swapHerb2bytes(&(gs2005->EvtDesc.DSMInput),1);
  swapHerb2bytes(&(gs2005->EvtDesc.externalBusy),1);
  swapHerb2bytes(&(gs2005->EvtDesc.modifiedBusyStatus),1);
  swapHerb2bytes(&(gs2005->EvtDesc.physicsWord),1);
  swapHerb2bytes(&(gs2005->EvtDesc.TriggerWord),1);
  swapHerb2bytes(&(gs2005->EvtDesc.DSMAddress),1);
  swapHerb2bytes(&(gs2005->EvtDesc.contaminationBusyStatus),1);
  swapHerb2bytes(&(gs2005->EvtDesc.npre),1);
  swapHerb2bytes(&(gs2005->EvtDesc.npost),1);
  swapHerb2bytes(&(gs2005->EvtDesc.dummy),1);

  swapHerb2bytes(&(gs2005->TrgSum.TrgSumBytes),1);
  swapHerb2bytes(&(gs2005->TrgSum.TrgSumHeader),1);
  swapHerb4bytes(&(gs2005->TrgSum.L1Sum[0]),2);
  swapHerb4bytes(&(gs2005->TrgSum.L2Sum[0]),2);
  swapHerb2bytes(&(gs2005->TrgSum.L0SumBytes),1);
  swapHerb2bytes(&(gs2005->TrgSum.L0SumHeader),1);
  swapHerb2bytes(&(gs2005->TrgSum.DSMdata.CPA[0]),32);
  swapHerb2bytes(&(gs2005->TrgSum.DSMdata.quadDSM[0]),8);
  swapHerb2bytes(&(gs2005->TrgSum.DSMdata.lastDSM[0]),8);
  swapHerb2bytes(&(gs2005->TrgSum.DSMdata.VTX[0]),8);
  swapHerb2bytes(&(gs2005->TrgSum.DSMdata.EMC[0]),8);
  swapHerb2bytes(&(gs2005->TrgSum.DSMdata.BCdata[0]),16);
  swapHerb2bytes(&(gs2005->TrgSum.DSMdata.specialTriggers[0]),8);
  swapHerb2bytes(&(gs2005->TrgSum.DSMdata.FPD[0]),8);
  swapHerb2bytes(&(gs2005->TrgSum.L1SumBytes),1);
  swapHerb2bytes(&(gs2005->TrgSum.L1SumHeader),1);
  swapHerb4bytes(&(gs2005->TrgSum.L1Result[0]),32);
  swapHerb2bytes(&(gs2005->TrgSum.L2SumBytes),1);
  swapHerb2bytes(&(gs2005->TrgSum.L2SumHeader),1);
  swapHerb4bytes(&(gs2005->TrgSum.L2Result[0]),32);
  
  numToSwap=1+gs2005->EvtDesc.npost+gs2005->EvtDesc.npre; assert(numToSwap<50&&numToSwap>0);
  assert(numToSwap>=0&&numToSwap<=PREPOST);
  int maxToSwap = 4*header.BankLength-(sizeof(EvtDescData2005)+sizeof(TrgSumData2005));
  maxToSwap/=sizeof(RawTrgDet2005);
  if (numToSwap>maxToSwap)  return -1;

    for(i=0;i<numToSwap;i++) { // loop over NPRE, NPOST as well
    swapHerb2bytes(&(gs2005->RAW[i].RawDetBytes),1);
    swapHerb2bytes(&(gs2005->RAW[i].CTBdataBytes),1);
    swapHerb2bytes(&(gs2005->RAW[i].MWCdataBytes),1);
    swapHerb4bytes(&(gs2005->RAW[i].MWCfiller),1);
    swapHerb2bytes(&(gs2005->RAW[i].BEMCdataBytes),1);
    swapHerb4bytes(&(gs2005->RAW[i].BEMCfiller),1);
    swapHerb2bytes(&(gs2005->RAW[i].BEMClayer1[0]),48);
    swapHerb2bytes(&(gs2005->RAW[i].EEMCdataBytes),1);
    swapHerb4bytes(&(gs2005->RAW[i].EEMCfiller),1);
    swapHerb2bytes(&(gs2005->RAW[i].EEMClayer1[0]),16);
    swapHerb2bytes(&(gs2005->RAW[i].FPDdataBytes),1);
    swapHerb4bytes(&(gs2005->RAW[i].FPDfiller),1);
    swapHerb2bytes(&(gs2005->RAW[i].FPDEastNSLayer1[0]),8);
    swapHerb2bytes(&(gs2005->RAW[i].FPDEastTBLayer1[0]),8);
    swapHerb2bytes(&(gs2005->RAW[i].FPDWestNSLayer1[0]),8);
    swapHerb2bytes(&(gs2005->RAW[i].FPDWestTBLayer1[0]),8);
    swapHerb2bytes(&(gs2005->RAW[i].BBCdataBytes),1);
    swapHerb4bytes(&(gs2005->RAW[i].BBCfiller),1);
    swapHerb2bytes(&(gs2005->RAW[i].BBClayer1[0]),16);
    swapHerb2bytes(&(gs2005->RAW[i].ZDClayer1[0]),8);
  }

  return returnValue;
}

/*!
  Helps ensure we actually have the trigger data.

  If one of these asserts()s fails, it probably means that the
  trigger group wrote the .daq bank with a new version of trgStructures.h .
  If so, then you will have to modify all the offline code to switch
  between versions of trgStructures.h .

*/
void TRG_Reader::SanityCheck2005(char *ptr, int check_s=1) { 

  gs2005=(MarilynMonroe*)ptr;
  unsigned short x;

  x=gs2005->TrgSum.L1SumBytes; assert(x==0x0084||x==0x8400);
  x=gs2005->TrgSum.L2SumBytes; assert(x==0x0084||x==0x8400);

  if (check_s){
    assert( gs2005->RAW[0].RawDetHeader[0]  =='R');   
    assert( gs2005->RAW[0].RawDetHeader[1]  =='D');   
    assert( gs2005->RAW[0].CTBdataHeader[0] =='C');   
    assert( gs2005->RAW[0].CTBdataHeader[1] =='T');   
    assert( gs2005->RAW[0].MWCdataHeader[0] =='M');   
    assert( gs2005->RAW[0].MWCdataHeader[1] =='W');   
    assert( gs2005->RAW[0].BEMCdataHeader[0]=='E'); 
    assert( gs2005->RAW[0].BEMCdataHeader[1]=='M');  
  } else {
    cout << "TRG_Reader::SanityCheck2005 : Data position sanity check is disabled" << endl;
  }
}
