//
// Pibero Djawotho <pibero@comp.tamu.edu>
// Texas A&M University Cyclotron Institute
// 7 Jan 2009
//

#include "DSM.hh"
#include "sumTriggerPatchChannels.hh"
#include "DSMAlgo_BE001_2009.hh"

void DSMAlgo_BE001_2009::operator()(DSM& dsm)
{
  // INPUT:

  // 10 x 12-bit BEMC channels
  // (5-0) high-tower
  // (11-6) trigger-patch

  // REGISTERS:

  // R0: BEMC-High-Tower-Th0 (6)
  // R1: BEMC-High-Tower-Th1 (6)
  // R2: BEMC-High-Tower-Th2 (6)
  // R3: BEMC-High-Tower-Th3 (6)

  int highTowerBits = 0;
  int  lowEtaSum = 0;
  int highEtaSum = 0;

  // Args: dsm, chMin, chMax, step, targetPedestal, sum, highTowerBits

  sumTriggerPatchChannels(dsm, 0, 5, 1, 1,  lowEtaSum, highTowerBits);
  sumTriggerPatchChannels(dsm, 6, 9, 1, 1, highEtaSum, highTowerBits);

  // OUTPUT (16):

  // (0-5) TP sum for low-eta group (6)
  // (6-11) TP sum for high-eta group (6)
  // (12-15) HT bits (4)

  int out = 0;

  out |=  lowEtaSum;
  out |= highEtaSum    <<  6;
  out |= highTowerBits << 12;

  dsm.output = out;
}
