//
// Pibero Djawotho <pibero@comp.tamu.edu>
// Texas A&M University Cyclotron Institute
// 7 Jan 2009
//

#include "DSM.hh"
#include "sumTriggerPatchChannels.hh"
#include "DSMAlgo_BW003_2009.hh"

void DSMAlgo_BW003_2009::operator()(DSM& dsm)
{
  // INPUT:

  // 10 12-bit BEMC channels
  // (5-0)  high-tower
  // (11-6) trigger-patch

  // J1 - ch0/2/4/6/8 (even channels) -> (16-31) 
  // J3 - ch1/3/5/7/9 (odd  channels) -> ( 0-15)  

  // REGISTERS:

  // R0: BEMC-High-Tower-Th0 (6)
  // R1: BEMC-High-Tower-Th1 (6)
  // R2: BEMC-High-Tower-Th2 (6)
  // R3: BEMC-High-Tower-Th3 (6)

  // ACTION:

  // JP1 - odd channels - to upper DSM

  int highTowerBitsJP1 = 0;
  int  lowEtaSumJP1 = 0;
  int highEtaSumJP1 = 0;

  // Args: dsm, chMin, chMax, step, targetPedestal, sum, highTowerBits

  sumTriggerPatchChannels(dsm, 1, 3, 2, 1,  lowEtaSumJP1, highTowerBitsJP1);
  sumTriggerPatchChannels(dsm, 5, 9, 2, 1, highEtaSumJP1, highTowerBitsJP1);

  // JP6 - even channels - to lower DSM

  int highTowerBitsJP6 = 0;
  int lowEtaSumJP6 = 0;
  int highEtaSumJP6 = 0;

  // Args: dsm, chMin, chMax, step, targetPedestal, sum, highTowerBits

  sumTriggerPatchChannels(dsm, 0, 2, 2, 1,  lowEtaSumJP6, highTowerBitsJP6);
  sumTriggerPatchChannels(dsm, 4, 8, 2, 1, highEtaSumJP6, highTowerBitsJP6);

  // OUTPUT (32):

  // JP1 (0-15) to upper DSM

  // (0-5) TP sum for low-eta group (6)
  // (6-11) TP sum for high-eta group (6)
  // (12-15) HT bits (4)

  // JP6 (16-31) to lower DSM

  // (16-21) TP sum for low-eta group (6)
  // (22-27) TP sum for high-eta group (6)
  // (28-31) HT bits (4)

  // JP1 (0-15)

  int out = 0;

  out |=  lowEtaSumJP1;
  out |= highEtaSumJP1    <<  6;
  out |= highTowerBitsJP1 << 12;

  // JP6 (16-31)

  out |=  lowEtaSumJP6    << 16;
  out |= highEtaSumJP6    << 22;
  out |= highTowerBitsJP6 << 28;

  dsm.output = out;
}
