//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2011      University of Warwick, UK
//
// Module: EvtExternalGenFactory
//
// Description: A factory type method to create engines for external physics
// generators like Pythia.
//
// Modification history:
//
//    John Back       April 2011            Module created
//
//------------------------------------------------------------------------------
//

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenExternal/EvtExternalGenFactory.hh"

#ifdef EVTGEN_PYTHIA
#include "EvtGenExternal/EvtPythiaEngine.hh"
#endif

#ifdef EVTGEN_PHOTOS
#include "EvtGenExternal/EvtPhotosEngine.hh"
#endif

#ifdef EVTGEN_TAUOLA
#include "EvtGenExternal/EvtTauolaEngine.hh"
#endif

#include <iostream>
using std::endl;

EvtExternalGenFactory::EvtExternalGenFactory() {

  _extGenMap.clear();

}

EvtExternalGenFactory::~EvtExternalGenFactory() {

  ExtGenMap::iterator iter;
  for (iter = _extGenMap.begin(); iter != _extGenMap.end(); ++iter) {

    EvtAbsExternalGen* theGenerator = iter->second;
    delete theGenerator;

  }
  
  _extGenMap.clear();

}

EvtExternalGenFactory* EvtExternalGenFactory::getInstance() {

  static EvtExternalGenFactory* theFactory = 0;
  
  if (theFactory == 0) {
    theFactory = new EvtExternalGenFactory();
  }

  return theFactory;

}

void EvtExternalGenFactory::definePythiaGenerator(std::string xmlDir, 
						  bool convertPhysCodes,
						  bool useEvtGenRandom) {

  // Only define the generator if we have the external ifdef variable set
#ifdef EVTGEN_PYTHIA

  int genId = EvtExternalGenFactory::PythiaGenId;

  EvtGenReport(EVTGEN_INFO,"EvtGen")<<"Defining EvtPythiaEngine: data tables defined in "
		       <<xmlDir<<endl;

  if (convertPhysCodes == true) {
    EvtGenReport(EVTGEN_INFO,"EvtGen")<<"Pythia 6 codes in decay files will be converted to Pythia 8 codes"<<endl;
  } else {
    EvtGenReport(EVTGEN_INFO,"EvtGen")<<"Pythia 8 codes need to be used in decay files"<<endl;
  }

  if (useEvtGenRandom == true) {
    EvtGenReport(EVTGEN_INFO,"EvtGen")<<"Using EvtGen random engine for Pythia 8 as well"<<endl;
  }

  EvtAbsExternalGen* pythiaGenerator = new EvtPythiaEngine(xmlDir, convertPhysCodes, useEvtGenRandom);
  _extGenMap[genId] = pythiaGenerator;

#endif

}

void EvtExternalGenFactory::definePhotosGenerator(std::string photonType, bool useEvtGenRandom) {

#ifdef EVTGEN_PHOTOS

  int genId = EvtExternalGenFactory::PhotosGenId;
  EvtGenReport(EVTGEN_INFO,"EvtGen")<<"Defining EvtPhotosEngine using photonType = "<<photonType<<endl;
  EvtAbsExternalGen* photosGenerator = new EvtPhotosEngine(photonType, useEvtGenRandom);
  _extGenMap[genId] = photosGenerator;

#endif

}

void EvtExternalGenFactory::defineTauolaGenerator(bool useEvtGenRandom) {

#ifdef EVTGEN_TAUOLA

  int genId = EvtExternalGenFactory::TauolaGenId;
  EvtGenReport(EVTGEN_INFO,"EvtGen")<<"Defining EvtTauolaEngine."<<endl;
  EvtAbsExternalGen* tauolaGenerator = new EvtTauolaEngine(useEvtGenRandom);
  _extGenMap[genId] = tauolaGenerator;

#endif

}

EvtAbsExternalGen* EvtExternalGenFactory::getGenerator(int genId) {

  EvtAbsExternalGen* theGenerator(0);

  ExtGenMap::iterator iter;

  if ((iter = _extGenMap.find(genId)) != _extGenMap.end()) {

    // Retrieve the external generator engine
    theGenerator = iter->second;

  }

  return theGenerator;

}

void EvtExternalGenFactory::initialiseAllGenerators() {

  ExtGenMap::iterator iter;
  for (iter = _extGenMap.begin(); iter != _extGenMap.end(); ++iter) {

    EvtAbsExternalGen* theGenerator = iter->second;
    if (theGenerator != 0) {
      theGenerator->initialise();
    }

  }
  
}
