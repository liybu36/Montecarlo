#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "DSGeneratorRDMDecayGunMessenger.hh"
#include "DSGeneratorRDMDecayGun.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

using namespace std;

DSGeneratorRDMDecayGunMessenger::DSGeneratorRDMDecayGunMessenger(DSGeneratorRDMDecayGun* gen){
  generator = gen;

  fDirectory = new G4UIdirectory("/ds/generator/rdm/");
  fDirectory->SetGuidance("Control of DSRDMDecayGun event generator");

  fIonCmd = new DSGeneratorRDMUIcmdWithNucleusAndUnit("/ds/generator/rdm/ion",this);
  fIonCmd->SetGuidance("define the primary ion (a,z,e)");
  fIonCmd->SetParameterName("A","Z","E",true);
  fIonCmd->SetDefaultUnit("keV");
  fIonCmd->SetUnitCandidates("keV MeV");




}


DSGeneratorRDMDecayGunMessenger::~DSGeneratorRDMDecayGunMessenger() {

  delete fDirectory;
  delete fIonCmd;
}


void DSGeneratorRDMDecayGunMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
   if (cmd == fIonCmd) {
     generator->SetNucleus(fIonCmd->GetNewNucleusValue(newValue));
   }
 }
