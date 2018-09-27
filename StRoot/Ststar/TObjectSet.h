// @(#)root/star:$Name:  $:$Id: TObjectSet.h,v 1.1.1.1 2000/05/19 12:46:09 fisyak Exp $
// Author: Valery Fine(fine@bnl.gov)   25/12/98

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifndef ROOT_TObjectSet
#define ROOT_TObjectSet

#include "TDataSet.h"

//////////////////////////////////////////////////////////////////////////////////////
//                                                                                  //
//  TObjectSet  - is a container TDataSet                                           //
//                  This means this object has an extra pointer to an embedded      //
//                  TObject.                                                        //
//  Terminology:    This TOvjectSet may be an OWNER of the embeded TObject          //
//                  If the container is the owner it can delete the embeded object  //
//                  otherwsie it leaves that object "as is"                         //
//                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////

class TObjectSet : public TDataSet {
protected:
  enum EOwnerBits { kIsOwner         = BIT(23) };
  TObject *fObj;                              // TObject to be inserted
  virtual Bool_t DoOwner(Bool_t done=kTRUE);

public:
  TObjectSet(const Char_t *name, TObject *obj=0,Bool_t makeOwner=kTRUE);
  TObjectSet(TObject *obj=0,Bool_t makeOwner=kTRUE);
  virtual ~TObjectSet();
  virtual void     Browse(TBrowser *b);
  virtual void     Delete(Option_t *opt="");
  virtual TObject *GetObject() const {return fObj;};
  virtual void     SetObject(TObject *obj) { SetObject(obj,kTRUE);}
  virtual TObject *SetObject(TObject *obj,Bool_t makeOwner);
  virtual TObject *AddObject(TObject *obj,Bool_t makeOwner=kTRUE);
  virtual Long_t   HasData() const {return fObj ? 1 : 0;}
  virtual Bool_t   IsOwner() const {return TestBit(kIsOwner);}

  ClassDef(TObjectSet,1)
};

#endif

