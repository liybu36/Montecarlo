#include "DSGeneratorSCSMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"



DSGeneratorSCSMessenger::DSGeneratorSCSMessenger( DSGeneratorSCS* gen ){

  fGenerator = gen;
  fDirectory = new G4UIdirectory("/ds/generator/scs/");
  fDirectory->SetGuidance("Control of Special Cross Sections generator");

  fIsotopeCmd = new G4UIcmdWithAString("/ds/generator/scs/isotope",this);
  G4String ioncandidates = "Ar39";
  fIsotopeCmd->SetCandidates( ioncandidates );

}


DSGeneratorSCSMessenger::~DSGeneratorSCSMessenger(){
  delete fDirectory;
  delete fIsotopeCmd;
}


void DSGeneratorSCSMessenger::SetNewValue( G4UIcommand* cmd, G4String newValue ){

  if(cmd == fIsotopeCmd){
    fGenerator->SetIsotope( newValue );
  }

}
