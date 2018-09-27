// @(#)root/star:$Name:  $:$Id: TFileSet.h,v 1.1.1.2 2000/09/05 21:08:48 fisyak Exp $
// Author: Valery Fine(fine@mail.cern.ch)   03/07/98

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFileSet
#define ROOT_TFileSet

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFileSet                                                             //
//                                                                      //
// TFileSet class is a class to convert the                             //
//      "native file system structure"                                  //
// into an instance of the TDataSet class                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TDataSet.h"
#include "TString.h"

class TFileSet : public TDataSet
{
 public:
    TFileSet();
    TFileSet(const TString &dirname, const Char_t *filename=0,Bool_t expand=kTRUE);
    virtual ~TFileSet();
    virtual Long_t HasData() const;
    virtual Bool_t IsEmpty() const;
    virtual Bool_t IsFolder() const;
    ClassDef(TFileSet,1)
};

#endif
