/***************************************************************************
 *
 * $Id: StTrackTopologyMap.cxx,v 2.20 2016/02/24 18:51:09 ullrich Exp $
 *
 * Author: Thomas Ullrich, Aug 1999
 ***************************************************************************
 *
 * Description:
 *
 * The topo map is stored in two 32 bit words.
 * Cases: Depending on the context the bits have slightly
 *        different meaning. For FTPC tracks bit 1-20 describe
 *        the FTPC. There's no confusion with the TPC since
 *        there are only tracks in either one. Starting 2007
 *        the 6 layers of the SVT are re-used for the HFT.
 *      
 * 
 * bit             Case 1               Case 2            Case 3
 *--------------------------------------------------------------------
 *  0            primary-vertex-used
 *  1            SVT layer=1      FTPC-West row=1     PXL layer=1
 *  2            SVT layer=2      FTPC-West row=2     PXL layer=2
 *  3            SVT layer=3      FTPC-West row=3     PXL layer=3
 *  4            SVT layer=4      FTPC-West row=4     IST layer=1
 *  5            SVT layer=5      FTPC-West row=5     IST layer=2
 *  6            SVT layer=6      FTPC-West row=6     SSD layer=1
 *  7            SSD              FTPC-West row=7     SSD layer=2
 *  8            TPC row=1        FTPC-West row=8
 *  9            TPC row=2        FTPC-West row=9
 *  10           TPC row=3        FTPC-West row=10
 *  11           TPC row=4        FTPC-East row=1
 *  12           TPC row=5        FTPC-East row=2
 *  13           TPC row=6        FTPC-East row=3
 *  14           TPC row=7        FTPC-East row=4
 *  15           TPC row=8        FTPC-East row=5
 *  16           TPC row=9        FTPC-East row=6
 *  17           TPC row=10       FTPC-East row=7
 *  18           TPC row=11       FTPC-East row=8
 *  19           TPC row=12       FTPC-East row=9
 *  20           TPC row=13       FTPC-East row=10
 *  21           TPC row=14
 *  22           TPC row=15
 *  23           TPC row=16
 *  24           TPC row=17
 *  25           TPC row=18
 *  26           TPC row=19
 *  27           TPC row=20
 *  28           TPC row=21
 *  29           TPC row=22
 *  30           TPC row=23
 *  31           TPC row=24
 *  -------------------------- word boundary
 *  32     0     TPC row=25
 *  33     1     TPC row=26
 *  34     2     TPC row=27
 *  35     3     TPC row=28
 *  36     4     TPC row=29
 *  37     5     TPC row=30
 *  38     6     TPC row=31
 *  39     7     TPC row=32
 *  40     8     TPC row=33
 *  41     9     TPC row=34
 *  42     10    TPC row=35
 *  43     11    TPC row=36
 *  44     12    TPC row=37
 *  45     13    TPC row=38
 *  46     14    TPC row=39
 *  47     15    TPC row=40
 *  48     16    TPC row=41
 *  49     17    TPC row=42
 *  50     18    TPC row=43
 *  51     19    TPC row=44
 *  52     20    TPC row=45
 *  53     21    Mwpc
 *  54     22    CTB
 *  55     23    ToF
 *  56     24    RICH
 *  57     25    Barrel EMC/SMD
 *  58     26    Endcap EMC/SMD
 *  59     27
 *  60     28
 *  61     29    HFT Format (case 3) - TPC tracks
 *  62     30    turn around flag  (flags that track spirals back)
 *  63     31    FTPC Format (flags TOC or FTPC)
 *
 ***************************************************************************
 *
 * $Log: StTrackTopologyMap.cxx,v $
 * Revision 2.20  2016/02/24 18:51:09  ullrich
 * Made kSsd and kSst equivalent in numberOfHits()
 *
 * Revision 2.19  2014/03/18 13:23:04  fisyak
 * Xin\'s fix for Ssd numberOfHits  with hftFormat
 *
 * Revision 2.18  2014/03/16 16:06:24  fisyak
 * Xin\'s fix for HFT
 *
 * Revision 2.17  2007/11/07 00:54:54  ullrich
 * Added PXL and IST.
 *
 * Revision 2.16  2005/06/23 19:04:24  ullrich
 * Added overloaded version of hasHitInDetector() taking up to 6 args.
 *
 * Revision 2.15  2004/08/11 05:21:22  genevb
 * tpcOnly excludes SSD
 *
 * Revision 2.14  2002/03/18 17:11:31  jeromel
 * Typo found by Jamie while testing P02gd. Corrected.
 *
 * Revision 2.13  2001/05/10 19:12:15  genevb
 * Switch FTPC definitions
 *
 * Revision 2.12  2001/04/24 21:32:07  genevb
 * Additional helper functions
 *
 * Revision 2.11  2001/04/05 04:00:58  ullrich
 * Replaced all (U)Long_t by (U)Int_t and all redundant ROOT typedefs.
 *
 * Revision 2.10  2000/12/08 20:21:08  genevb
 * Changed kTofPatchId -> kTofId
 *
 * Revision 2.9  2000/07/28 19:49:28  akio
 * Change in Detector Id for Endcap SMD
 *
 * Revision 2.8  2000/05/17 17:21:34  ullrich
 * New method largestGap() and new output operator.
 *
 * Revision 2.7  2000/04/12 19:44:03  genevb
 * Reimplement mMap data members as individual unsigned ints
 *
 * Revision 2.6  2000/04/10 19:59:33  genevb
 * StRoot/StEvent/doc/tex/
 *
 * Revision 2.5  2000/03/29 00:16:30  ullrich
 * Fixed off-by-one error.
 *
 * Revision 2.4  2000/01/07 18:20:34  ullrich
 * Buf fixed. Wrong indices for mMap in bit().
 *
 * Revision 2.3  1999/12/13 20:16:36  ullrich
 * Changed numbering scheme for hw_position unpack methods (STAR conventions).
 *
 * Revision 2.2  1999/10/28 22:27:55  ullrich
 * Adapted new StArray version. First version to compile on Linux and Sun.
 *
 * Revision 2.1  1999/10/13 19:45:46  ullrich
 * Initial Revision
 *
 **************************************************************************/
#include "StTrackTopologyMap.h"
#include <vector>
#include <algorithm>
#include <numeric>
#if !defined(ST_NO_NAMESPACES)
using std::vector;
using std::adjacent_difference;
using std::max_element;
#endif

static const char rcsid[] = "$Id: StTrackTopologyMap.cxx,v 2.20 2016/02/24 18:51:09 ullrich Exp $";

ClassImp(StTrackTopologyMap)

StTrackTopologyMap::StTrackTopologyMap()
{
    mMap0 = mMap1 = 0;
}

StTrackTopologyMap::StTrackTopologyMap(unsigned int m1, unsigned int m2) : mMap0(m1), mMap1(m2) { /* noop */ }

StTrackTopologyMap::StTrackTopologyMap(const unsigned int* m) : mMap0(m[0]), mMap1(m[1]) { /* noop */ }

StTrackTopologyMap::StTrackTopologyMap(const unsigned long* m) : mMap0(m[0]), mMap1(m[1]) { /* noop */ }

StTrackTopologyMap::~StTrackTopologyMap() { /* noop */ }

bool
StTrackTopologyMap::bit(int i) const
{
    return i>31 ? (mMap1>>(i-32) & 1U) : (mMap0>>i & 1U);
}

bool
StTrackTopologyMap::ftpcFormat() const
{
    return bit(63);
}

bool
StTrackTopologyMap::hftFormat() const
{
    return bit(61);
}

unsigned int
StTrackTopologyMap::data(unsigned int i) const
{
    return static_cast<unsigned int>(i<2 ? (i<1 ? mMap0 : mMap1) : 0);
}

bool
StTrackTopologyMap::primaryVertexUsed() const { return bit(0); }

bool
StTrackTopologyMap::turnAroundFlag() const { return bit(62); }

bool
StTrackTopologyMap::hasHitInDetector(StDetectorId id) const
{
    return ((numberOfHits(id)) ? 1U : 0U);
}


bool
StTrackTopologyMap::hasHitInDetector(StDetectorId d1, StDetectorId d2,
                               StDetectorId d3, StDetectorId d4,
                               StDetectorId d5, StDetectorId d6) const
{
    //
    //  Note d3 - d6 are optional, if not given they will have
    //  the value kUnknownId and shouldn't be used.
    //
    return (hasHitInDetector(d1) && hasHitInDetector(d2) &&
          (d3 == kUnknownId ? true : hasHitInDetector(d3)) &&
          (d4 == kUnknownId ? true : hasHitInDetector(d4)) &&
          (d5 == kUnknownId ? true : hasHitInDetector(d5)) &&
          (d6 == kUnknownId ? true : hasHitInDetector(d6)));
}

bool
StTrackTopologyMap::hasHitInSvtLayer(unsigned int layer) const
{
    if (ftpcFormat())
        return false;
    else
        return bit(layer);
}

bool
StTrackTopologyMap::hasHitInPxlLayer(unsigned int layer) const
{
    if(ftpcFormat())
        return false;
    else
        return bit(layer);
}

bool
StTrackTopologyMap::hasHitInIstLayer(unsigned int layer) const
{
    if(ftpcFormat())
        return false;
    else
//        return bit(layer+2);
        return bit(layer+3);
}

bool
StTrackTopologyMap::hasHitInSsdLayer(unsigned int layer) const
{
    if(ftpcFormat())
        return false;
    else {
      if(hftFormat())
        return bit(layer+5);
      else
        return bit(layer+6);
    }
}

bool
StTrackTopologyMap::hasHitInRow(StDetectorId id, unsigned int row) const
{
    switch (id) {
    case kTpcId:
        return !ftpcFormat() && bit(row+7);
        break;
    case kFtpcWestId:
        return ftpcFormat() && bit(row);
        break;
    case kFtpcEastId:
        return ftpcFormat() && bit(row+10);
        break;
    default:
        return false;
        break;
    }
}

unsigned int
StTrackTopologyMap::numberOfHits(StDetectorId id) const
{
    if (ftpcFormat() &&
        !(id == kFtpcWestId || id == kFtpcEastId))
        return 0;
    
    int i;
    int n = 0;

    switch (id) {
    case kSvtId:
        for (i=1; i<7; i++)
            if (hasHitInSvtLayer(i)) n++;
        break;
    case kSstId:
    case kSsdId:
        if(! hftFormat()) {
          if (bit(7)) n++;        
        } else {
          for(int i=1;i<=2;i++) {
            if (hasHitInSsdLayer(i)) n++;
          }
        }
        break;
    case kPxlId:
        for(i=1;i<=3;i++)
	  if(hasHitInPxlLayer(i)) n++;
        break;
    case kIstId:
//        for(i=1;i<4;i++)
        for(i=1;i<=2;i++)
	  if(hasHitInIstLayer(i)) n++;
        break;
    case kFtpcWestId:
    case kFtpcEastId:
        for (i=1; i<11; i++)
            if (hasHitInRow(id, i)) n++;
        break;
    case kTpcId:
        for (i=1; i<46; i++)
            if (hasHitInRow(id, i)) n++;
        break;
    case kMwpcWestId:
    case kMwpcEastId:
        if (bit(53)) n++;
        break;
    case kCtbId:
        if (bit(54)) n++;
        break;
    case kTofId:
        if (bit(55)) n++;
        break;
    case kRichId:
        if (bit(56)) n++;
        break;
    case kBarrelEmcTowerId:
    case kBarrelEmcPreShowerId:
    case kBarrelSmdEtaStripId:
    case kBarrelSmdPhiStripId:
        if (bit(57)) n++;
        break;
    case kEndcapEmcTowerId:
    case kEndcapEmcPreShowerId:
    case kEndcapSmdUStripId:
    case kEndcapSmdVStripId:
        if (bit(58)) n++;
        break;
    default:
        n = 0;
        break;
    }
    return n;
}

bool
StTrackTopologyMap::trackTpcOnly() const
{
    if(hftFormat())
      return ((hasHitInDetector(kTpcId)) & 
             ~((hasHitInDetector(kPxlId)) | (hasHitInDetector(kIstId)) | (hasHitInDetector(kSsdId))));
    else
      return ((hasHitInDetector(kTpcId)) &
             ~((hasHitInDetector(kSvtId)) | (hasHitInDetector(kSsdId))));
}

bool
StTrackTopologyMap::trackSvtOnly() const
{
    return ((hasHitInDetector(kSvtId)) & ~(hasHitInDetector(kTpcId)));
}

bool
StTrackTopologyMap::trackTpcSvt() const
{
    return ((hasHitInDetector(kTpcId)) & (hasHitInDetector(kSvtId)));
}

bool
StTrackTopologyMap::trackFtpcEast() const
{
    return (hasHitInDetector(kFtpcEastId));
}

bool
StTrackTopologyMap::trackFtpcWest() const
{
    return (hasHitInDetector(kFtpcWestId));
}

bool
StTrackTopologyMap::trackFtpc() const
{
    return ((hasHitInDetector(kFtpcWestId)) | (hasHitInDetector(kFtpcEastId)));
}

int
StTrackTopologyMap::largestGap(StDetectorId id) const
{
    if (ftpcFormat() && !(id == kFtpcWestId || id == kFtpcEastId))
        return -1;

    vector<int> rows;
    int i;

    switch (id) {
    case kSvtId:
        for (i=1; i<7; i++)
            if (hasHitInSvtLayer(i)) rows.push_back(i);
        break;
    case kFtpcWestId:
    case kFtpcEastId:
        for (i=1; i<11; i++)
            if (hasHitInRow(id, i)) rows.push_back(i);
        break;
    case kTpcId:
        for (i=1; i<46; i++)
            if (hasHitInRow(id, i)) rows.push_back(i);
        break;
    default:
        return -1;
    }

    if (rows.size() < 2) return -1;

    vector<int> diffs(rows.size());
    adjacent_difference(rows.begin(), rows.end(), diffs.begin());
    return *max_element(diffs.begin()+1, diffs.end()) - 1;  // skip first
}

ostream& operator<< (ostream& os, const StTrackTopologyMap& m)
{
    for (int i=0; i<64; i++) {
        if (i>31)
            os << ((m.data(1)>>(i-32) & 1U) ? 1 : 0);
        else
            os << ((m.data(0)>>i & 1U) ? 1 : 0);
    }
    return os;
}

