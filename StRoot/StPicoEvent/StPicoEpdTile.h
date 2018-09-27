#ifndef StPicoEpdTile_h
#define StPicoEpdTile_h

#include <cstdlib>
#include "TObject.h"

#include "StPicoEvent/StPicoCommon.h"

using namespace StarPicoDst;


/**
 * EPD tile class for STAR picoDst
 * Total size of this object is 6 bytes
 *
 * If all EPD tiles are saved, this makes 2*12*31*6 = 4.4 kB
 *
 * 1/8 install in 2017:  if all EPD tiles are saved, this makes 1*3*31*6 = 558 bytes
 * in 2017, we will save all EPD tiles.  After 2017, we will not.
 *
 * - Mike Lisa 20 May 2017
 */
//_________________
class StPicoEpdTile : public TObject {

 public:
  //Default constructor
  StPicoEpdTile();
  //Constructor that takes values
  StPicoEpdTile(int positionId, int tileId, DetectorSide EW, int ADC, int TAC, int TDC, bool hasTAC, bool statusIsGood = true);
  //Copy constructor
  StPicoEpdTile(const StPicoEpdTile &tile);
  //Destructor
  virtual ~StPicoEpdTile();
  //Print EPD tile information
  virtual void Print(const Char_t *option = "") const;

  bool hasTac() const;
  int  adc() const;
  int  tac() const;
  int  tdc() const;
  DetectorSide side() const;

  int id() const;
  int position() const;     // 1...12
  int tile() const;         // 1...31

  /// false if tile is bad or missing
  bool isGood() const;

 protected:
  
  /// Packed channel Id: 100*positionId + tileId
  /// sign: +/- = West/East
  /// positionId and tileId are phototube indices start at 1, [1, 12] and [1, 31] respectively
  Short_t mId;

  /// Packed channel data: bits  0-11 are ADC; bits 12-23 are TAC;
  ///                      bits 24-28 are TDC; bit 29 is noTAC flag
  ///                      bit 30 is the good/bad (1/0) status flag
  ULong_t mQTdata;

  ClassDef(StPicoEpdTile, 1)
};

inline DetectorSide StPicoEpdTile::side() const { return mId < 0 ? DetectorSide::East : DetectorSide::West;}
inline int  StPicoEpdTile::id() const { return mId; }
inline int  StPicoEpdTile::position() const { return std::abs(mId / 100); }
inline int  StPicoEpdTile::tile() const { return std::abs(mId % 100); }
inline int  StPicoEpdTile::adc() const { return mQTdata & 0x0FFF; }
inline int  StPicoEpdTile::tac() const { return (mQTdata >> 12) & 0x0FFF; }
inline int  StPicoEpdTile::tdc() const { return (mQTdata >> 24) & 0x001F; }
inline bool StPicoEpdTile::hasTac() const { return (mQTdata >> 29) & 0x1; }
inline bool StPicoEpdTile::isGood() const { return (mQTdata >> 30) & 0x1; }
#endif
