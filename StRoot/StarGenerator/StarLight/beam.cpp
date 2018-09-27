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


#include <iostream>
#include <fstream>
#include <cmath>

#include "starlightconstants.h"
#include "inputParameters.h"
#include "reportingUtils.h"
#include "bessel.h"
#include "beam.h"


using namespace std;
using namespace starlightConstants;


//______________________________________________________________________________
beam::beam(const int              Z,
           const int              A,
           const double           bdeuteron,
           const bool             dAuCoherentProduction,
           const inputParameters& input)
	: nucleus(Z,
	          A,
	          bdeuteron,
	          dAuCoherentProduction)
{
  // setting needed inputparameters to protected variables
	_beamLorentzGamma = input.beamLorentzGamma();
}


//______________________________________________________________________________
beam::~beam()
{ }


//______________________________________________________________________________
double beam::photonFlux(const double impactparameter, 
                        const double photonEnergy) const
{
  // function for the calculation of the "photon density".
  // photonFlux = number of photons / (energy * area)
  // assume beta = 1 and gamma >> 1, i.e. neglect the (1 / gamma^2) * K_0(x) term
  
  const double X
	  = (impactparameter * photonEnergy) / (_beamLorentzGamma * hbarc);
  if (X <= 0) 
	  printWarn << "X = " << X << endl;
  
  const double factor1 = (double(Z() * Z()) * alpha) / (pi * pi);  
  const double factor2 = 1. / (photonEnergy * impactparameter * impactparameter);
  const double bessel  = bessel::dbesk1(X);
  const double factor3 = X * X * bessel * bessel;

  return factor1 * factor2 * factor3;
}
