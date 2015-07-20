#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "DSGeneratorAmBeSourceMessenger.hh"
#include "DSGeneratorAmBeSource.hh"
#include "DSEventHandler.hh"
#include "DSLogger.hh"

using namespace std;

DSGeneratorAmBeSourceMessenger::DSGeneratorAmBeSourceMessenger(DSGeneratorAmBeSource* gen){
  
  generator = gen;
  fDirectory = new G4UIdirectory("/ds/generator/AmBe/");
  fDirectory->SetGuidance("Control of DSAmBeSource event generator");

  fNeutrinoCmd = new G4UIcmdWithAString("/ds/generator/AmBe/source",this);
  G4String candidates = "all neutron0G neutron1G neutron2G";
  fNeutrinoCmd->SetCandidates(candidates);

  fDisableParticleType = new G4UIcmdWithAString("/ds/generator/AmBe/disable",this);
  fDisableParticleType->SetCandidates("n gamma");

}

DSGeneratorAmBeSourceMessenger::~DSGeneratorAmBeSourceMessenger() {

  delete fDirectory;
  delete fNeutrinoCmd;
  delete fDisableParticleType;
}

void DSGeneratorAmBeSourceMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) {

  if (cmd == fNeutrinoCmd) {
    DSLog(routine) << "AmBeSource source: " << newValue << endlog;
    if (newValue == "all") {
     generator->SetNeutrinoType(DSGeneratorAmBeSource::AllPN);
    }
    else if (newValue == "neutron0G") {
     generator->SetNeutrinoType(DSGeneratorAmBeSource::NeutronNoG);
    }
    else if (newValue == "neutron1G") {
     generator->SetNeutrinoType(DSGeneratorAmBeSource::Neutron1G);
    }
    else if (newValue == "neutron2G") {
     generator->SetNeutrinoType(DSGeneratorAmBeSource::Neutron2G);
    }
    else {
      DSLog(error) << "invalid value: " << newValue << ", should be one of 'all', 'neutron0G', 'neutron1G' or 'neutron2G'" << endlog;
    }
  }

  if (cmd == fDisableParticleType) {
    DSLog(routine) << "AmBeSource disable particle type: " << newValue << endlog;
    if (newValue == "n") {
     generator->SetDisableParticleType(1);
    }
    else if (newValue == "gamma") {
     generator->SetDisableParticleType(2);
    }
    else {
      DSLog(error) << "invalid value: " << newValue << ", should be either 'n' or 'gamma'" << endlog;
    }
  }
    
}
