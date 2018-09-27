//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: EvtRanf.cpp,v 1.1 2016/10/28 16:52:48 jwebb Exp $
//
// Description:
//	subroutine evtranf_.
//      Provides FORTRAN calable interface to EvtRandom::Flat()
//      Can be used as EVTRANF instead of RANF in FORTRAN programs
//      or as evtranf_ instead of ranf_ in C/C++ programs.
//      No header file is provided, as C++ programs should use EvtRandom
//
//
//------------------------------------------------------------------------
#include "EvtGenBase/EvtRandom.hh"
#ifdef WIN32
extern "C" {
  double __stdcall EVTRANF( ) 
  {
    return EvtRandom::Flat() ;
  }
}
#else
extern "C" {
  double evtranf_( ) 
  {
    return EvtRandom::Flat() ;
  }
}
#endif
