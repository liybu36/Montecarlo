#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "DSGeneratorEnergyDepositMessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "DSLogger.hh"
#include <stdio.h>

using namespace std;

DSGeneratorEnergyDepositMessenger::DSGeneratorEnergyDepositMessenger(DSGeneratorEnergyDeposit* gen){
  generator = gen;
  fDirectory = new G4UIdirectory("/ds/generator/endep/");
  fDirectory->SetGuidance("Control of DSEnergyDeposit generator");
  
  fFileCmd  = new G4UIcmdWithAString("/ds/generator/endep/file",this);
}


DSGeneratorEnergyDepositMessenger::~DSGeneratorEnergyDepositMessenger() {

  delete fDirectory;
  delete fFileCmd;
}


void DSGeneratorEnergyDepositMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
  if (cmd == fFileCmd){
    DSIO::Get()->SetG4DSFile(newValue);
    DSIO::Get()->SetIsG4DS(true);
  }
}


/*
 * Revision 1.1  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 *
 */
