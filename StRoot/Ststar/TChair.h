// @(#)root/star:$Name: v2-25-03vp $:$Id: TChair.h,v 1.7 2000/09/22 22:25:13 fine Exp $
// Author: Valery Fine(fine@bnl.gov)   13/03/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TChair
#define ROOT_TChair

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TChair                                                              //
//                                                                      //
//  It is a base class to create a custom interface for TTable objects  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TTable.h"

class TChair : public TDataSet {

private:
   TTable *fTable;

protected:
   ULong_t  fLastIndx;  // index pof the last used  table row;
   void    *fLastRow;   // pointer to the last used table row; fLastRow = table[fLastIndx]

   TTable    *GetThisTable() const {return fTable; }
   static void  *GetOffset(const void *base,ULong_t offset) { return (void  *)((Char_t *)base + offset);}
   TChair(){ fTable = 0; }

public:

   TChair(TTable *table){ fTable = table; }
   TChair(const TChair &chair){assert(0);}
//   TChair     &operator=(const TChair &rhs){ assert(0); return rhs;}
   virtual    ~TChair(){;}

   virtual     void       Adopt(Int_t n, void *array){GetThisTable()->Adopt(n,array);}
   virtual     void       AddAt(const void *c, Int_t i){GetThisTable()->AddAt(c,i);}
              const void *At(Int_t i) const {return GetThisTable()->At(i);}
   virtual     void       Browse(TBrowser *b){GetThisTable()->Browse(b);}
   virtual     void       CopySet(TChair &chair){GetThisTable()->CopySet(*chair.GetThisTable());}
               Int_t      CopyRows(const TChair *srcChair, Int_t srcRow=0, Int_t dstRow=0, Int_t nRows=0, Bool_t expand=kFALSE)
                          {return GetThisTable()->CopyRows(srcChair->GetThisTable(),srcRow,dstRow,nRows,expand);}
   virtual     void       Draw(Option_t *opt){GetThisTable()->Draw(opt);}
   virtual     TH1       *Draw(TCut varexp, TCut selection, Option_t *option="",
                          Int_t nentries=1000000000, Int_t firstentry=0)
                          {return GetThisTable()->Draw(varexp,selection,option,nentries,firstentry);}
   virtual     TH1       *Draw(const Text_t *varexp, const Text_t *selection, Option_t *option="",
                               Int_t nentries=1000000000, Int_t firstentry=0) {
                           return GetThisTable()->Draw(varexp,selection,option,nentries,firstentry);}
   virtual     Char_t    *GetArray() const    {return (Char_t *)GetThisTable()->GetArray();}
   virtual     TClass    *GetRowClass() const {return GetThisTable()->GetRowClass();}
   virtual     Long_t     GetNRows() const    {return GetThisTable()->GetNRows();}
   virtual     Long_t     GetRowSize() const  {return GetThisTable()->GetRowSize();}
   virtual     Long_t     GetTableSize() const{return GetThisTable()->GetTableSize();}
               const TTable  *Table() const {return fTable; }
   virtual     TTableDescriptor *GetRowDescriptors()   const {return GetThisTable()->GetRowDescriptors();}
   virtual     const Char_t       *GetType()             const {return GetThisTable()->GetType();}
   virtual     void       Fit(const Text_t *formula ,const Text_t *varexp, const Text_t *selection="",Option_t *option="",Option_t *goption="",
                              Int_t nentries=1000000000, Int_t firstentry=0) {
                           GetThisTable()->Fit(formula,varexp,selection,option,goption,nentries,firstentry);}
   virtual     Long_t     HasData() const  { return GetThisTable()->HasData();}
   virtual     Bool_t     IsFolder() const { return GetThisTable()->IsFolder();}
   virtual     void       ls(Option_t *option=""){GetThisTable()->ls(option);}
   virtual     void       ls(Int_t deep)  {GetThisTable()->ls(deep);}
               Int_t      NaN()           {return GetThisTable()->NaN();}
   virtual     Char_t    *MakeExpression(const Char_t *expressions[],Int_t nExpressions)
                         {return GetThisTable()->MakeExpression(expressions,nExpressions);}
   virtual     Char_t    *Print(Char_t *buf,Int_t n) const { return GetThisTable()->Print(buf, n);}
   virtual     void       Print(Option_t *opt="")          {GetThisTable()->Print(opt);}
   virtual  const Char_t *Print(Int_t row, Int_t rownumber=10,
                                const Char_t *colfirst="",const Char_t *collast="") const {
                           return GetThisTable()->Print(row,rownumber,colfirst,collast); }
   virtual  const Char_t *PrintHeader() const {return GetThisTable()->PrintHeader();}
   virtual  Int_t         Purge(Option_t *opt="")    {return GetThisTable()->Purge(opt);}

               void      *ReAllocate(Int_t newsize) { return GetThisTable()->ReAllocate(newsize); }
               void      *ReAllocate()              { return GetThisTable()->ReAllocate(); }
   virtual     void       SavePrimitive(ofstream &out, Option_t *option) {GetThisTable()->SavePrimitive(out,option);}

   virtual     void       Set(Int_t n)                                   {GetThisTable()->Set(n);}
   virtual     void       Set(Int_t n, Char_t *array)                    {GetThisTable()->Set(n,array);}
   virtual     void       SetNRows(Int_t n)                              {GetThisTable()->SetNRows(n);}
   virtual     void       Reset(Int_t c=0)                               {GetThisTable()->Reset(c) ;}
   virtual     void       Update()                                       {GetThisTable()->Update();}
   virtual     void       Update(TDataSet *set, UInt_t opt=0)            {GetThisTable()->Update(set,opt);}
               void      *operator[](Int_t i);
              const void *operator[](Int_t i) const;

   ClassDef(TChair,0)  // A base class to provide a user custom interface to TTable class objects
};

inline void *TChair::operator[](Int_t i)
{

//   if (!GetThisTable()->BoundsOk("TChair::operator[]", i))
//      i = 0;
    return (void *)((char *)GetArray()+i*GetRowSize());
}

inline const void *TChair::operator[](Int_t i) const
{
//   if (!GetThisTable()->BoundsOk("TChair::operator[]", i))
//      i = 0;
    return (const void *)((char *)GetArray()+i*GetRowSize());
}
#endif
