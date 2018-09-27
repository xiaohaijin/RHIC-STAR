//
// Pibero Djawotho <pibero@comp.tamu.edu>
// Texas A&M University Cyclotron Institute
// 7 Jan 2009
//

#ifndef SUM_TRIGGER_PATCH_CHANNELS_HH
#define SUM_TRIGGER_PATCH_CHANNELS_HH

#include "DSM.hh"

inline void sumTriggerPatchChannels(const DSM& dsm, int chMin, int chMax, int step, int targetPedestal, int& sum, int& highTowerBits)
{
  // Loop over channels, make ADC sum and compare high towers to thresholds

  int nChannels = 0;

  for (int ch = chMin; ch <= chMax; ch += step) {
    int highTower = dsm.channels[ch] & 0x3f;

    for (int reg = 0; reg < 4; ++reg)
      highTowerBits |= (highTower > dsm.registers[reg]) << reg;

    int triggerPatch = dsm.channels[ch] >> 6 & 0x3f;

    sum += triggerPatch;

    ++nChannels;
  }

  // Reset pedestal to 1

  if (sum < nChannels)
    sum = 0;
  else
    sum -= (nChannels - targetPedestal);

  // If overflow, set to max

  if (sum > 63) sum = 63;
}

#endif	// SUM_TRIGGER_PATCH_CHANNELS_HH
