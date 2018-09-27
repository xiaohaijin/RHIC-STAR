///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author: jwebb $: author of last commit
// $Date: 2012/11/27 22:27:31 $: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef EVENTCHANNEL_H
#define EVENTCHANNEL_H

#include <vector>

#include "starlightconstants.h"
#include "readinluminosity.h"
#include "beambeamsystem.h"
#include "randomgenerator.h"
#include "upcevent.h"



class eventChannel : public readLuminosity
{

public:

	eventChannel(const inputParameters& input,
	             beamBeamSystem&        bbsystem);
	virtual ~eventChannel();

	unsigned long nmbAttempts() const {return _nmbAttempts;}  ///< returns number of attempted events
	unsigned long nmbAccepted() const {return _nmbAccepted;}  ///< returns number of accepted events

	virtual starlightConstants::event produceEvent(int &ievent) = 0;

	virtual upcEvent produceEvent() = 0;
 
	static void transform(const double betax,
	                      const double betay,
	                      const double betaz,
	                      double&      E,
	                      double&      px,
	                      double&      py,
	                      double&      pz,
	                      int&         iFbadevent); ///< Lorentz-transforms given 4-vector

	randomGenerator _randy;
	beamBeamSystem  _bbs;

	static double pseudoRapidity(const double px,
	                             const double py,
	                             const double pz);  ///< calculates pseudorapidity for given 3-momentum

protected:

	unsigned long _nmbAttempts;  ///< number of attempted events
	unsigned long _nmbAccepted;  ///< number of accepted events

	bool   _ptCutEnabled;   ///< en/disables cut in pt
	double _ptCutMin;       ///< minimum pt, if cut is enabled
	double _ptCutMax;       ///< maximum pt, if cut is enabled
	bool   _etaCutEnabled;  ///< en/disables cut in eta
	double _etaCutMin;      ///< minimum eta, if cut is enabled
	double _etaCutMax;      ///< maximum eta, if cut is enabled
		
};


#endif  // EVENTCHANNEL_H
