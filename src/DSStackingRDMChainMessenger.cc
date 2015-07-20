//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// *********--***********************************************************
//

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "DSStackingRDMChainMessenger.hh"
#include "DSStackingRDMChain.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "DSStorage.hh"
#include "DSLogger.hh"
using namespace std;

DSStackingRDMChainMessenger::DSStackingRDMChainMessenger(DSStackingRDMChain* stack){
  stacking = stack;
  fDirectory = new G4UIdirectory("/ds/stack/rdmchain/");
  
  
  fLifeTimeCmd = new G4UIcmdWithADoubleAndUnit("/ds/stack/rdmchain/maxlifetime",this);
  fLifeTimeCmd->SetGuidance("Set the max life time acceptable for a decay");
  fLifeTimeCmd->SetDefaultUnit("s");
  fLifeTimeCmd->SetUnitCandidates("ps ns mus ms s d y");
  


}


DSStackingRDMChainMessenger::~DSStackingRDMChainMessenger() {
  delete fDirectory;
  delete fLifeTimeCmd;
}


void DSStackingRDMChainMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
   if (cmd == fLifeTimeCmd){
       stacking->SetMaxLifeTime(fLifeTimeCmd->GetNewDoubleValue(newValue));
       DSLog(routine) << "Max Decay Mean Life < " << newValue <<  endlog ;
   }
}
