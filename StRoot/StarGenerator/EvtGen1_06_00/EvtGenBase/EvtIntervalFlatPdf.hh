/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtIntervalFlatPdf.hh,v 1.1 2016/09/23 18:37:30 jwebb Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

/*
 * Uniform PDF defined on a 1D interval.
 */

#ifndef EVT_INTERVAL_FLAT_PDF_HH
#define EVT_INTERVAL_FLAT_PDF_HH

#include <assert.h>
#include "EvtGenBase/EvtPdf.hh"
#include "EvtGenBase/EvtPoint1D.hh"

class EvtIntervalFlatPdf : public EvtPdf<EvtPoint1D> {
public:
  
  EvtIntervalFlatPdf(double min, double max);
  EvtIntervalFlatPdf(const EvtIntervalFlatPdf& other);
  virtual ~EvtIntervalFlatPdf();
  virtual EvtPdf<EvtPoint1D>* clone() const;
  
  virtual EvtValError compute_integral() const;
  virtual EvtPoint1D randomPoint();
  
protected:

  virtual double pdf(const EvtPoint1D&) const;

  double _min;
  double _max;
};

#endif
