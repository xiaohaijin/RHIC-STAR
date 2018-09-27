#ifndef _ST_DAQ_SIM_ADC_TB_H
#define _ST_DAQ_SIM_ADC_TB_H
#ifndef __CINT__
#include "DAQ_READER/daq_dta_structs.h"
#else
class daq_sim_adc_tb;
#endif
#include "TTable.h"
#include "Ttypes.h"
class St_daq_sim_adc_tb : public TTable {
 public:
  ClassDefTable(St_daq_sim_adc_tb,daq_sim_adc_tb)
  ClassDef(St_daq_sim_adc_tb,1) // 
};
#endif
