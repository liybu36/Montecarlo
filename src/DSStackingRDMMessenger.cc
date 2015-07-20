#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "DSStackingRDMMessenger.hh"
#include "DSStackingRDM.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "DSEventHandler.hh"
#include "DSLogger.hh"
#include "G4Tokenizer.hh"
//using namespace std;

DSStackingRDMMessenger::DSStackingRDMMessenger(DSStackingRDM* stack){
  stacking = stack;
  
  fDirectory = new G4UIdirectory("/ds/stack/rdm/");
  


  fKillParticleCmd = new G4UIcmdWithAnInteger("/ds/stack/rdm/kill",this);

  fKillLEParticleCmd = new G4UIcommand("/ds/stack/rdm/killLE",this);
  fKillLEParticleCmd->SetGuidance("Set properties of ion to be generated.");
  fKillLEParticleCmd->SetGuidance("[usage]/ds/stack/rdm/killLE  pdg E");
  fKillLEParticleCmd->SetGuidance("	 pdg:(int) particle code");
  fKillLEParticleCmd->SetGuidance("	 E:(double) kinetic energy (in keV)");

  G4UIparameter* param;
  param = new G4UIparameter("pdg",'i',false);
  param->SetDefaultValue("-1000");
  fKillLEParticleCmd->SetParameter(param);
  param = new G4UIparameter("E",'d',true);
  param->SetDefaultValue("0.0");
  fKillLEParticleCmd->SetParameter(param);
 
}


DSStackingRDMMessenger::~DSStackingRDMMessenger() {

  delete fDirectory;
  delete fKillParticleCmd;
  delete fKillLEParticleCmd;
 
}


void DSStackingRDMMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
  if (cmd == fKillParticleCmd){
      stacking->KillParticles(fKillParticleCmd->ConvertToInt(newValue));
      DSLog(routine) << "Kill particle with PDG code: " << newValue   << endlog ;
   }  else if( cmd == fKillLEParticleCmd) {
      G4Tokenizer next( newValue );
      G4int PDG   = StoI(next());
      G4double ene = StoD(next())*keV ;
      stacking->KillLEParticles(PDG,ene); 
      DSLog(routine) << "Kill particle with PDG code: " <<  PDG 
       << " and kinetic energy below " << ene/keV << " keV"   << endlog ;
    }
}
