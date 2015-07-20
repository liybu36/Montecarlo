#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "DSGeneratorMultiEventMessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "DSLogger.hh"
#include <stdio.h>

using namespace std;

DSGeneratorMultiEventMessenger::DSGeneratorMultiEventMessenger(DSGeneratorMultiEvent* gen){
  generator = gen;
  fDirectory = new G4UIdirectory("/ds/generator/multi/");
  fDirectory->SetGuidance("Control of DSMultiEvent generator");
  
  fPDGCmd  = new G4UIcmdWithAString("/ds/generator/multi/event",this);
}


DSGeneratorMultiEventMessenger::~DSGeneratorMultiEventMessenger() {

  delete fDirectory;
  delete fPDGCmd;
}


void DSGeneratorMultiEventMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
  if (cmd == fPDGCmd){
    generator->SetListOfParticles(newValue);
  }
}


/*
 * Revision 1.1  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 *
 */
