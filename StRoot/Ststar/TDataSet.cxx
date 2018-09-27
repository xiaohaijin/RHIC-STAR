// @(#)root/star:$Name: v2-25-03vp $:$Id: TDataSet.cxx,v 1.4 2000/09/29 18:20:42 fine Exp $
// Author: Valery Fine(fine@mail.cern.ch)   03/07/98
const char *gCoPyRiGhT[] = {
     "STAR dataset C++ base class library:",
     "Copyright(c) 1997~2000  [BNL] Brookhaven National Laboratory, Valeri Fine  (fine@bnl.gov). All right reserved",
     "************************************************************************",
     "This program is distributed in the hope that it will be useful,",
     "but WITHOUT ANY WARRANTY; without even the implied warranty of",
     "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.",
     "",
     "Permission to use, copy, modify and distribute this software and its",
     "documentation for any purpose is hereby granted without fee,",
     "provided that the above copyright notice appear in all copies and",
     "that both that copyright notice and this permission notice appear",
     "in supporting documentation.  Brookhaven National Laboratory makes no",
     "representations about the suitability of this software for any",
     "purpose.  It is provided \"as is\" without express or implied warranty.",
     "************************************************************************"
};

const char *Id = {
    "$Id: TDataSet.cxx,v 1.4 2000/09/29 18:20:42 fine Exp $"
};
#include <iostream.h>
#include "TSystem.h"
#include "TDataSetIter.h"
#include "TDataSet.h"

#include "TROOT.h"
#include "TBrowser.h"

#include "TSystem.h"

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TDataSet                                                             //
//                                                                      //
// TDataSet class is to create a special compound object-container:     //
//                                                                      //
// ==================================================================== //
//    TDataSet object ::= the "named" list of TDataSet objects          //
// ==================================================================== //
// where the "list" (the pointer to TList object) may contain no object //
//                                                                      //
//  TDataSet object has a back pointer to its "parent" TDataSet         //
//  object, the "character" *name* and "character" *title*              //
//                                                                      //
//  The service this class does provide is to help the user to build    //
//  and manage the hierarchy of his/her data but the data itself.       //
//                                                                      //
//  So it is not "Container" itself rather the basement (base class)    //
//  to built the containers.                                            //
//                                                                      //
//  One may derive the custom container classes from TDataSet.          //
//  See for example TObjectSet, TTable, TVolume, TFileSet               //
//  These classes  derived from TDataSet:                               //
//                                                                      //
//   Class Name                                                         //
//   ----------                                                         //
//  TObjectSet::public TDataSet - is a container for TObject            //
//  TTable::    public TDataSet - is a container for the array          //
//                                    of any "plain" C-structure        //
//  TNode::     public TDataSet - is a container for 3D objects         //
//  TMaker::     public TDataSet - is a container for STAR "control"    //
//                                    objects                           //
//   etc etc                                                            //
//                                                                      //
//  TDataSet class is a base class to implement the directory-like      //
//  data structures and maintain it via TDataSetIter class iterator     //
//                                                                      //
// TDataSet can be iterated using an iterator object (see TDataSetIter) //
//            or by TDataSet::Pass method (see below)                   //
//                                                                      //
//  Terms:    Dataset       - any object from the list above            //
//  =====     Member          is called "DataSet Member"                //
//                                                                      //
//          Structural      - the "Dataset Member" is its               //
//            member          "Structural member" if its "back pointer" //
//                            points to this object                     //
//                                                                      //
//           Dataset        - we will say this TDataSet object "OWNs"   //
//            Owner           (or is an OWNER / PARENT of ) another     //
//          (parent)          TDataSet object if the last one is its    //
//                            "Structural Member"                       //
//                                                                      //
//          Associated      - If some object is not "Structural member" //
//            member          of this object we will say it is an       //
//                            "Associated Member" of this dataset       //
//                                                                      //
//           Orphan         - If some dataset is a member of NO other   //
//           dataset          TDataSet object it is called an "orphan"  //
//                            dataset object                            //
//                                                                      //
// - Any TDataSet object may be "Owned" by one and only one another     //
//   TDataSet object if any.                                            //
//                                                                      //
// - Any TDataSet object can be the "Structural Member" of one and      //
//   only one another TDataSet                                          //
//                                                                      //
// - Any TDataSet object may be an "Associated Member" for any number   //
//   of other TDataSet objects if any                                   //
//                                                                      //
// - NAME issue:                                                        //
//   Each "dataset member" is in possession of some "alpha-numerical"   //
//   NAME as defined by TNamed class.                                   //
//   The NAME may contain any "printable" symbols but "SLASH" - "/"     //
//   The symbol "RIGHT SLASH" - "/" can not be used as any part of the  //
//   "DataSet Member" NAME                                              //
//    Any DataSet  can be found by its NAME with TDataSetIter object    //
//                                                                      //
// - TITLE issue:                                                       //
//   Each "dataset member" is in possession of the "alpha-numerical"    //
//   TITLE as defined by TNamed class. The meaning of the TITLE is      //
//   reserved for the derived classes to hold there some indetification //
//   that is special for that derived class.                            //
//                                                                      //
//   This means the user must be careful about  the "TDataSet           //
//   NAME and TITLE since this may cause some "side effects" of the     //
//   particular class functions                                         //
//                                                                      //
// - It is NOT required those all "DataSet Members" are in possession   //
//   of the unique names, i.e. any number of "DataSet Members"          //
//   may bear one and the same name                                     //
//                                                                      //
//   Actions:                                                           //
//   ========                                                           //
//   Create  DataSet is born either as "Orphan" or                      //
//                                  as "Structural Member"              //
//           of another TDataSet object                                 //
//                                                                      //
//   Add     One dataset can be included into another dataset.          //
//           Upon adding:                                               //
//           -  the "Orphan dataset" becomes "Structural Member"        //
//           - "Structural Members" of another dataset becomes the      //
//             "Associated Member" of this datatset                     //
//                                                                      //
//   Delete  - Upon deleting the "Structural Member":                   //
//             - "REMOVES" itself  from the "Parent DataSet".           //
//             - Its "Associated memberships" is not changed though     //
//                                                                      //
//              The last means the DataSet with the "Associated Members"//
//              may contain a DIED pointers to unexisting "Associated"  //
//              objects !!!                                             //
//                                                                      //
//  Further information is provided my the particular method            //
//  descriptions.                                                       //
//                                                                      //
//  The TDataSet class has several methods to control object('s)        //
//  memberships                                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

TDataSet mainSet("DSMAIN");
TDataSet *TDataSet::fgMainSet = &mainSet;

ClassImp(TDataSet)

//______________________________________________________________________________
TDataSet::TDataSet(const Char_t *name, TDataSet *parent, Bool_t arrayFlag)
           : TNamed(name,"TDataSet")
{
  //  cout << "ctor for " << GetName() << " - " << GetTitle() << endl;
   if (strchr(name,'/')) {
      Error("TDataSet::TDataSet","dataset name cannot contain a slash", name);
      return;
   }
   fList =0; fParent=0;
   if (arrayFlag) SetBit(kArray);
   if (parent) parent->Add(this);
//   else AddMain(this);
}

//______________________________________________________________________________
TDataSet *TDataSet::GetRealParent(){
  TDataSet *p = GetParent();
  if (fgMainSet && p == fgMainSet) p = 0;
  return p;
}
//______________________________________________________________________________
TDataSet::TDataSet(const TDataSet &pattern,EDataSetPass iopt)
{
  //
  // Creates TDataSet (clone) with a topology similar with TDataSet *pattern
  //
  //  Parameters:
  //  -----------
  //  pattern        - the pattern dataset
  //  iopt = kStruct - clone only my structural links
  //         kAll    - clone all links
  //         kRefs   - clone only refs
  //         kMarked - clone marked (not implemented yet) only
  //
  //   All new-created sets become the structural ones anyway.
  //
  //  cout << "ctor for " << GetName() << " - " << GetTitle() << endl;

  SetName(pattern.GetName());
  SetTitle(pattern.GetTitle());

  TDataSet *set = 0;
  TDataSetIter next((TDataSet *)&pattern);
  Bool_t optsel = (iopt == kStruct);
  Bool_t optall = (iopt == kAll);
  while ((set = next())) {
//              define the parent of the next set
     TDataSet *parent = set->GetParent();
     if ( optall || (optsel && parent == this) )
       Add((TDataSet *)(set->Clone()));
  }
}
//______________________________________________________________________________
TDataSet::~TDataSet()
{
//              cout << "Default destructor for " << GetName() << " - " << GetTitle() << endl;
   Shunt(0); Delete();
}

//______________________________________________________________________________
void  TDataSet::MakeCollection()
{
  if (!fList)
    fList = TestBit(kArray) ? (TSeqCollection *)new TObjArray : (TSeqCollection *) new TList;
}

//______________________________________________________________________________
void TDataSet::AddAt(TDataSet *dataset,Int_t idx)
{
//
// Add TDataSet object at the "idx" position in ds
// or at the end of the dataset
// The final result is defined by either TList::AddAt or TObjArray::AddAt
// methods
//
  if (!dataset) return;

  MakeCollection();

  // Check whether this new child has got any parent yet
  if (!dataset->GetRealParent()) dataset->SetParent(this);
  fList->AddAt(dataset,idx);
}

//______________________________________________________________________________
void TDataSet::AddAtAndExpand(TDataSet *dataset, Int_t idx)
{
//   !!!! Under construction !!!!!
// Add TDataSet object at the "idx" position in ds
// or at the end of the dataset
// The final result is defined by either TList::AddAt or TObjArray::AddAt
// methods
//
  if (!dataset) return;

  MakeCollection();

  // Check whether this new child has got any parent yet
  if (!dataset->GetRealParent()) dataset->SetParent(this);
  if (TestBit(kArray)) ((TObjArray *) fList)->AddAtAndExpand(dataset,idx);
  else                  fList->AddAt(dataset,idx);
}

//______________________________________________________________________________
void TDataSet::AddLast(TDataSet *dataset)
{
// Add TDataSet object at the end of the dataset list of this dataset
  if (!dataset) return;

  MakeCollection();

  // Check whether this new child has got any parent yet
  if (!dataset->GetRealParent()) dataset->SetParent(this);
  fList->AddLast(dataset);
}

//______________________________________________________________________________
void TDataSet::AddFirst(TDataSet *dataset)
{
 // Add TDataSet object at the beginning of the dataset list of this dataset
  if (!dataset) return;

  MakeCollection();

  // Check whether this new child has got any partent yet
  if (!dataset->GetRealParent()) dataset->SetParent(this);
  fList->AddFirst(dataset);
}

//______________________________________________________________________________
void TDataSet::Browse(TBrowser *b)
{
  // Browse this dataset (called by TBrowser).
   TDataSetIter next(this);
   TDataSet *obj;
   if (b)
       while ((obj = next())) b->Add(obj,obj->GetName());
}

//______________________________________________________________________________
TObject *TDataSet::Clone() {
   return new TDataSet(*this);
}

//______________________________________________________________________________
void TDataSet::Delete(Option_t *opt)
{
//
// Delete - deletes the list of the TDataSet objects and all "Structural Members"
//          as well
//          This method doesn't affect the "Associated Members"
//
  if(opt){/*unused*/}

//      Delete list of the TDataSet
  TSeqCollection     *thisList = GetCollection();
  if (!thisList) return;
  fList = 0;
  TIter next(thisList);
  TDataSet *son = 0;
  //  Delete the "Structural Members" of this TDataSet only
  while ((son = (TDataSet *)next())) {
    if ( (!son->IsOnHeap()) || (this != son->GetParent()) ) continue;
    // mark the object is deleted from the TDataSet dtor or Delete method
    son->SetBit(kCanDelete);
    delete son;
  }
  //  Cleare list
  thisList->Clear("nodelete");
  delete thisList;
}

//______________________________________________________________________________
TDataSet  *TDataSet::FindByPath(const Char_t *path) const
{
  // Aliase for TDataSet::Find(const Char_t *path) method
   return Find(path);
}

//______________________________________________________________________________
TDataSet *TDataSet::Find(const Char_t *path) const
{
  //
  // Full description see: TDataSetIter::Find
  //
  // Note. This method is quite expansive.
  // ----- It is done to simplify the user's code when one wants to find ONLY object.
  //       If you need to find more then 1 object in this dataset,
  //       regard using TDataSetIter class yourself.
  //
  TDataSetIter next((TDataSet*)this);
  return next.Find(path);
}

//______________________________________________________________________________
TObject  *TDataSet::FindObject(const Char_t *name) const
{
  // This method is not recommended.
  // It is done to back TObject::FindObject method only.
  // One is recommnened to use FindByName method instead.

  return FindByName(name);
}

//______________________________________________________________________________
TDataSet *TDataSet::FindByName(const Char_t *name,const Char_t *path,Option_t *opt) const
{
  //
  // Full description see: TDataSetIter::FindDataSet
  //
  // Note. This is method is quite expansive.
  // ----- It is done to simplify the user's code when one wants to find ONLY object.
  //       If you need to find more then 1 object in this dataset,
  //       regard using TDataSetIter class yourself.
  //

  TDataSetIter next((TDataSet*)this);
  return next.FindByName(name,path,opt);
}

//______________________________________________________________________________
TDataSet *TDataSet::First() const
{
 //  Return the first object in the list. Returns 0 when list is empty.
 if (fList) return (TDataSet *)(fList->First());
 return 0;
}

//______________________________________________________________________________
TDataSet *TDataSet::Last() const
{
 // Return the last object in the list. Returns 0 when list is empty.
 if (fList) return (TDataSet *)(fList->Last());
 return 0;
}

//______________________________________________________________________________
void  TDataSet::ls(Option_t *option) const
{
 /////////////////////////////////////////////////////////////////////
 //                                                                 //
 //  ls(Option_t *option)                                           //
 //                                                                 //
 //    option       - defines the path to be listed                 //
 //           = "*" -  means print all levels                       //
 //                                                                 //
 /////////////////////////////////////////////////////////////////////

  if (option && !strcmp(option,"*")) ls(Int_t(0));
  else {
    TDataSet *set = 0;
    if (option && strlen(option) > 0) {
      TDataSetIter local((TDataSet*)this);
      set = local(option);
    }
    else
      set = (TDataSet*)this;
    if (set) set->ls(Int_t(1));
    else
      if (option) Warning("ls","Dataset <%s> not found",option);
  }
}

//______________________________________________________________________________
void TDataSet::ls(Int_t depth) const
{
 /////////////////////////////////////////////////////////////////////
 //                                                                 //
 //  ls(Int_t depth)                                                //
 //                                                                 //
 //  Prints the list of the this TDataSet.                          //
 //                                                                 //
 //  Parameter:                                                     //
 //  =========                                                      //
 //    Int_t depth >0 the number of levels to be printed            //
 //               =0 all levels will be printed                     //
 //            No par - ls() prints only level out                  //
 //                                                                 //
 /////////////////////////////////////////////////////////////////////

  printf("%3d - %s\t%s\n",TROOT::GetDirLevel(),(const char*)Path(),(char*)GetTitle());
  if (!fList || depth == 1 ) return;
  if (!depth) depth = 99999;

  TIter next(fList);
  TDataSet *d=0;
  while ((d = (TDataSet *)next())) {
    TROOT::IncreaseDirLevel();
    d->ls(depth-1);
    TROOT::DecreaseDirLevel();
  }
}
//______________________________________________________________________________
Bool_t TDataSet::IsThisDir(const Char_t *dirname,int len,int ignorecase) const
{
  if (!ignorecase) {
    if (len<0) {return !strcmp (GetName(),dirname);
    } else     {return !strncmp(GetName(),dirname,len);}
  } else {
    const char *name = GetName();
    if (len==-1) len = strlen(dirname);
    for (int i=0;i<len;i++) { if ( tolower(name[i])!=tolower(dirname[i])) return 0;}
    return 1;
  }
}

//______________________________________________________________________________
void TDataSet::MarkAll()
{
  // Mark all members of this dataset
  Mark();
  TDataSetIter nextMark(this,0);
  TDataSet *set = 0;
  while ( (set = nextMark()) ) set->Mark();
}

//______________________________________________________________________________
void TDataSet::UnMarkAll()
{
  // UnMark all members of this dataset
  Mark(kMark,kReset);
  TDataSetIter nextMark(this,0);
  TDataSet *set = 0;
  while ( (set = nextMark()) ) set->Mark(kMark,kReset);
}

//______________________________________________________________________________
void TDataSet::InvertAllMarks()
{
 // Invert mark bit for all members of this dataset
  if (IsMarked()) Mark(kMark,kReset);
  else Mark();
  TDataSetIter nextMark(this,0);
  TDataSet *set = 0;
  while (( set = nextMark()) ) {
   if (set->IsMarked()) set->Mark(kMark,kReset);
   else set->Mark();
 }
}

//______________________________________________________________________________
Bool_t TDataSet::IsEmpty() const
{
   // return kTRUE if the "internal" collection has no member
   return First() ? kFALSE : kTRUE ;
}

//______________________________________________________________________________
TString TDataSet::Path() const
{
 // return the full path of this data set
   TString str;
   TDataSet *parent = GetParent();
   if (parent) {
       str = parent->Path();
       str += "/";
   }
   str +=  GetName();
   return str;
}

//______________________________________________________________________________
void TDataSet::Remove(TDataSet *set)
{
  if (fList && set) {
    if (set->GetParent() == this) set->SetParent(0);
    fList->Remove(set);
  }

}

//______________________________________________________________________________
TDataSet  *TDataSet::RemoveAt(Int_t idx)
{
  //
  // Remove object from the "idx" cell of this set and return
  // the pointer to the removed object if any
  //
  TDataSet *set = 0;
  if (fList) {
      set = (TDataSet *)fList->At(idx);
      fList->RemoveAt(idx);
      if (set && (set->GetParent() == this) ) set->SetParent(0);
  }
  return set;
}

//______________________________________________________________________________
EDataSetPass TDataSet::Pass(EDataSetPass ( *callback)(TDataSet *),Int_t depth)
{
 /////////////////////////////////////////////////////////////////////
 //                                                                 //
 // Pass (callback,depth)                                           //
 //                                                                 //
 // Calls callback(this) for all datasets those recursively         //
 //                                                                 //
 //  Parameter:                                                     //
 //  =========                                                      //
 //    Int_t depth >0 the number of levels to be passed             //
 //                =0 all levels will be passed                     //
 //                                                                 //
 //  Return (this value mast be returned by the user's callback):   //
 //  ======                                                         //
 //  kContinue - continue passing                                   //
 //  kPrune    - stop passing the current branch, go to the next one//
 //  kUp       - stop passing, leave the current branch,            //
 //              return to previous level and continue              //
 //  kStop     - stop passing, leave all braches                    //
 //                                                                 //
 /////////////////////////////////////////////////////////////////////

  if (!callback) return kStop;

  EDataSetPass condition = callback(this);

  if (condition == kContinue){
    if (fList && depth != 1 ) {
      TIter next(fList);
      TDataSet *d=0;
      while ( (d = (TDataSet *)next()) ) {
         condition = d->Pass(callback, depth == 0 ? 0 : --depth);
         if (condition == kStop || condition == kUp) break;
      }
    }
  }
  return condition==kUp ? kContinue:condition;
}

//______________________________________________________________________________
EDataSetPass TDataSet::Pass(EDataSetPass ( *callback)(TDataSet *,void*),void *user,Int_t depth)
{
 /////////////////////////////////////////////////////////////////////
 //                                                                 //
 // Pass (callback,user,depth)                                      //
 //                                                                 //
 // Calls callback(this,user) for all datasets those recursively    //
 //                                                                 //
 //  Parameter:                                                     //
 //  =========                                                      //
 //    Int_t depth >0 the number of levels to be passed             //
 //                =0 all levels will be passed                     //
 //                                                                 //
 //  Return (this value mast be returned by the user's callback):   //
 //  ======                                                         //
 //  kContinue - continue passing                                   //
 //  kPrune    - stop passing the current branch, go to the next one//
 //  kUp       - stop passing, leave the current branch,            //
 //              return to previous level and continue              //
 //  kStop     - stop passing, leave all braches                    //
 //                                                                 //
 /////////////////////////////////////////////////////////////////////

  if (!callback) return kStop;

  EDataSetPass condition = callback(this,user);

  if (condition == kContinue){
    if (fList && depth != 1 ) {
      TIter next(fList);
      TDataSet *d=0;
      while ((d = (TDataSet *)next())) {
        condition = d->Pass(callback, user, depth == 0 ? 0 : --depth);
        if (condition == kStop) break;
        if (condition == kUp  ) break;
      }
    }
  }
  return (condition==kUp) ? kContinue:condition;
}

//______________________________________________________________________________
Int_t TDataSet::Purge(Option_t *)
{
//
// Purge  - deletes all "dummy" "Structural Members" those are not ended
//          up with some dataset with data inside (those return HasData() = 0)
//
// Purge does affect only the "Structural Members" and doesn't "Associated" ones
//

 if (!fList) return 0;
 TIter next(fList);
 TDataSet *son = 0;
   // Purge "Structural Members" only
   TList garbage;
   while ((son = (TDataSet *)next())) {
     if (this == son->GetParent()) continue;
     //   mark the object is deleted from the TDataSet dtor
     son->Purge();
     if (son->HasData() || son->GetListSize()) continue;
     delete son;
 }
 return 0;
}
//______________________________________________________________________________
void  TDataSet::SetParent(TDataSet *parent)
{
//
//  Break the "parent" relationship with the current object parent if present
//  parent != 0   Makes this object the "Structural Member"
//                of the "parent" dataset
//          = 0   Makes this object the "pure Associator", i.e it makes this
//                object the "Structural Member" of NO other TDataSet
//
   fParent = parent;
}

//______________________________________________________________________________
void TDataSet::SetWrite()
{
 // One should not use this method but TDataSet::Write instead
 // This method os left here for the sake of the backward compatibility
 // To Write object first we should temporary break the
 // the backward fParent pointer (otherwise ROOT follows this links
 // and will pull fParent out too.
 //
 Write();
}

//______________________________________________________________________________
void TDataSet::Shunt(TDataSet *newParent)
{
  //
  //  Remove the object from the original and add it to dataset
  //  TDataSet dataset   != 0  -  Make this object the "Structural Member"
  //                                of "dataset"
  //                        = 0  -  Make this object "Orphan"
  //
  if (fParent)   fParent->Remove(this);
  if (newParent) newParent->Add(this);
}

//______________________________________________________________________________
void TDataSet::Update(TDataSet* set,UInt_t opt)
{
//
// Update this TDataSet with "set"
//
// ATTENTION !!!
// ---------
// This method changes the parent relationships of the input "set"
//
  if(opt){/*unused*/}
  if(!set) return;

  SetTitle(set->GetTitle());
  TDataSetIter nextnew(set);
  TDataSet *newset = 0;
  while((newset = nextnew())) {
    Bool_t found = kFALSE;
//              Check whether this has the list of the sons
    if (fList) {
      TIter nextold(fList);
      const Char_t *newname = newset->GetName();
      TDataSet *oldset = 0;
      while ( ((oldset = (TDataSet *)nextold())!=0) && !found) {
        // if the "new" set does contain the dataset
        // with the same name as ours update it too
        if (oldset->IsThisDir(newname)) {
           oldset->Update(newset);
           found = kTRUE;
        }
      }
    }
    // If the new "set" contains some new dataset with brand-new name
    // move it into the our dataset and remove it from its old location
   if (!found) newset->Shunt(this);
  }
}

//______________________________________________________________________________
void TDataSet::Update()
{
 //
 //  Update()
 //
 //  Recursively updates all tables for all nested datasets
 //  in inverse order
 //

  TDataSetIter next(this);
  TDataSet *set = 0;
  while(( set = next())) set->Update();
}

//______________________________________________________________________________
void TDataSet::Sort()
{
  TDataSetIter next(this,0);
  TDataSet *ds;
  TList *list;
  while ((ds=next())) {
    list = ds->GetList();
    if (!list) continue;
    list->Sort(); ds->Sort();
  }
}

//______________________________________________________________________________
Int_t TDataSet::Write(const Text_t *name, Int_t option, Int_t bufsize)
{
 //
 // To Write object first we should temporary break the
 // the backward fParent pointer (otherwise ROOT follows this links
 // and will pull fParent out too.
 //
  TDataSet *saveParent = fParent; // GetParent();
  fParent = 0;
  Int_t nbytes = TObject::Write(name,option, bufsize);
  fParent = saveParent;
  return nbytes;
}


class Copyright {public: Copyright(const Char_t *id=0){
  if (!gSystem) return;
  if (gSystem->Getenv("No_Copyright") ||
     ( gSystem->Getenv("STAR") && !gSystem->Getenv("Copyright")) ) return;
  Int_t i = 0;
  cout <<  " " << gCoPyRiGhT[i++] << endl
       <<  " " << id << endl
       <<  " " << gCoPyRiGhT[i++] << endl;
  if (gSystem->Getenv("Copyright") && id) {
    Int_t lCopy = sizeof(gCoPyRiGhT)/4;
    for (Int_t l =i; l < lCopy; l++) cout << " " << gCoPyRiGhT[l] << endl;
  }}};
static Copyright TDataSetLibraryCopyright(Id);
