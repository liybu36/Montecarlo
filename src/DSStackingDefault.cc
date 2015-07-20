#include "DSStorage.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSEventHandler.hh"
#include "DSStackingDefault.hh"
#include "DSLogger.hh"
#include "G4StackManager.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4ParticleTable.hh"
#include "G4NeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4OpticalPhoton.hh" 
#include "G4GenericIon.hh" 
#include "globals.hh"
#include "CLHEP/Random/RandExponential.h"
#include "DSStorage.hh"

using namespace std;

DSStackingDefault::DSStackingDefault() {

 
}


DSStackingDefault::~DSStackingDefault(){;}

G4ClassificationOfNewTrack DSStackingDefault::DSClassifyNewTrack (const G4Track* aTrack) {
  if(aTrack->GetGlobalTime() > 10*ms) return fKill ;
  
  return fUrgent;
}

void DSStackingDefault::DSNewStage() { 


 
  
}

void DSStackingDefault::DSPrepareNewEvent() {  
}

/*
 * $Log: DSStackingDefault.cc,v $
 * Revision 1.1  2014/05/07 12:21:05  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.1  2013/07/07 09:42:17  dfranco
 * Added a default stacking for killing particles from long living nucleus decays
 *
 *
 */
