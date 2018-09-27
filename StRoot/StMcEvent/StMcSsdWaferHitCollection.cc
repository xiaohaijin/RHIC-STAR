/***************************************************************************
 *
 * $Id: StMcSsdWaferHitCollection.cc,v 2.2 2012/03/01 16:48:29 perev Exp $
 *
 * Author: Manuel Calderon de la Barca Sanchez, Oct 1999
 ***************************************************************************
 *
 * Description: Monte Carlo Ssd Wafer Hit Collection class
 *
 ***************************************************************************
 *
 * $Log: StMcSsdWaferHitCollection.cc,v $
 * Revision 2.2  2012/03/01 16:48:29  perev
 * method Browse() added
 *
 * Revision 2.1  2005/11/22 21:44:52  fisyak
 * Add compress Print for McEvent, add Ssd collections
 *
 * Revision 2.4  2005/01/27 23:40:48  calderon
 * Adding persistency to StMcEvent as a step for Virtual MonteCarlo.
 *
 * Revision 2.3  2003/02/19 03:16:05  calderon
 * Introduction of Ctb Hit Class and Ctb Hit Collection class, modified
 * StMcTrack, and StMcEvent accordingly.  Clearing of hits in StMcSsdWaferHitCollection.
 *
 * Revision 2.2  2000/03/06 18:05:23  calderon
 * 1) Modified SSD Hits storage scheme from layer-ladder-wafer to
 * barrel-ladder-wafer.
 * 2) Added Rich Hit class and collection, and links to them in other
 * classes.
 *
 * Revision 2.1  1999/11/19 19:06:33  calderon
 * Recommit after redoing the files.
 *
 * Revision 2.0  1999/11/17 02:01:00  calderon
 * Completely revised for new StEvent
 *
 *
 **************************************************************************/
#include "TBrowser.h"
#include "StMcSsdWaferHitCollection.hh"
#include "StMcSsdHit.hh"

static const char rcsid[] = "$Id: StMcSsdWaferHitCollection.cc,v 2.2 2012/03/01 16:48:29 perev Exp $";
ClassImp(StMcSsdWaferHitCollection);
//_____________________________________________________________________________
StMcSsdWaferHitCollection::StMcSsdWaferHitCollection() { /* noop */ }

//_____________________________________________________________________________
StMcSsdWaferHitCollection::~StMcSsdWaferHitCollection()
{
  // StMcSsdHit provides its own new/delete operator, and
  // mHits is a polymorphic container, so we need to do this.
  Clear();
}
//_____________________________________________________________________________
const StSPtrVecMcSsdHit& StMcSsdWaferHitCollection::hits() const { return mHits; }

//_____________________________________________________________________________
StSPtrVecMcSsdHit& StMcSsdWaferHitCollection::hits() { return mHits; }
//_____________________________________________________________________________
void StMcSsdWaferHitCollection::Clear(const char*)
{
  for (int i=0; i<(int)mHits.size(); i++) 
  {
    delete mHits[i]; mHits[i] = 0;
  }
  mHits.clear();
}
//_____________________________________________________________________________
void StMcSsdWaferHitCollection::Browse(TBrowser *b)
{
  // Browse this event (called by TBrowser).
   for (int i=0; i<(int)mHits.size(); i++) {
     TObject *obj = mHits[i]; if (!obj) continue;
     TString ts(obj->GetName()); ts+="#"; ts+=i;
     b->Add(obj,ts.Data());
   }
}
