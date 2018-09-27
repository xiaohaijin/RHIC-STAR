#ifndef DSM_ALGO_EM201_2013_HH
#define DSM_ALGO_EM201_2013_HH

#include "../DSMAlgo.hh"

struct DSMAlgo_EM201_2013 : public DSMAlgo {
  int ajpBarrel(DSM& dsm, int offset) const;
//  int ajpEndcap(const DSM& dsm) const;
  int Dijet(DSM& dsm, int jpTH) const;
  int EndDijet(DSM& dsm, int jpTH) const;

  void operator()(DSM& dsm);
};

#endif	// DSM_ALGO_EM201_2013_HH
