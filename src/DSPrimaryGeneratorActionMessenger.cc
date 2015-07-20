#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "DSGeneratorG4Gun.hh"
#include "DSGeneratorCosmicRayMuons.hh"
#include "DSGeneratorNeutronsAtGS.hh"
#include "DSGeneratorMultiEvent.hh"
#include "DSGeneratorEnergyDeposit.hh"
#include "DSGeneratorRDMDecayGun.hh"
#include "DSGeneratorSCS.hh"
#include "DSGeneratorAmBeSource.hh"
#include "DSPrimaryGeneratorActionMessenger.hh"
#include "DSPrimaryGeneratorAction.hh"
#include "DSLogger.hh"


using namespace std;

DSPrimaryGeneratorActionMessenger::DSPrimaryGeneratorActionMessenger(DSPrimaryGeneratorAction *act){
  fGeneratorPrimary = act;

  fDirectory = new G4UIdirectory("/ds/generator/");
  fDirectory->SetGuidance("Control commands for generators:");
  fDirectory->SetGuidance("/ds/generator/select: Select generator.");

  // /DS/generator/select command
  fSelectCmd = new G4UIcmdWithAString("/ds/generator/select", this);
  fSelectCmd->SetGuidance("Selects generator for events.");
  fSelectCmd->SetGuidance("Options are:");
  fSelectCmd->SetGuidance("G4gun: Standard G4 gun.");
  fSelectCmd->SetGuidance("CosmicRayMuons: cosmic-ray muons generator"); 
  fSelectCmd->SetGuidance("NeutronsAtGS: neutrons generator");
  fSelectCmd->SetGuidance("AmBeSource: AmBe generator");
 

  G4String candidates 
  = "G4Gun CosmicRayMuons NeutronsAtGS MultiEvent RDM SCS EnergyDeposit AmBeSource";
  fSelectCmd->SetCandidates(candidates);
 }


DSPrimaryGeneratorActionMessenger::~DSPrimaryGeneratorActionMessenger() {

  delete fDirectory;
  delete fGeneratorPrimary;

}


void DSPrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 
   if(command == fSelectCmd) {   
   
     DSLog(routine) << "Generator: " << newValue << endlog;
    
     if(newValue == "G4Gun") {
       fGeneratorPrimary->SetDSGenerator(new DSGeneratorG4Gun);        
     }  else if (newValue == "CosmicRayMuons") {     
       fGeneratorPrimary->SetDSGenerator(new DSGeneratorCosmicRayMuons);  
     }  else if (newValue == "NeutronsAtGS") {     
       fGeneratorPrimary->SetDSGenerator(new DSGeneratorNeutronsAtGS);  
     }  else if (newValue == "MultiEvent") {     
       fGeneratorPrimary->SetDSGenerator(new DSGeneratorMultiEvent); 
     }  else if (newValue == "EnergyDeposit") {     
       fGeneratorPrimary->SetDSGenerator(new DSGeneratorEnergyDeposit);  
     }  else if (newValue == "RDM") {     
       fGeneratorPrimary->SetDSGenerator(new DSGeneratorRDMDecayGun);  
     }  else if (newValue == "SCS") {
       fGeneratorPrimary->SetDSGenerator(new DSGeneratorSCS);
     }  else if (newValue == "AmBeSource") {
       fGeneratorPrimary->SetDSGenerator(new DSGeneratorAmBeSource);
     }
   }

}


/*
 * Revision 1.3  2013/03/22 13:50:09  dfranco
 * Added radioactive decay and radioactive decay chain generators. They require the correspondent stacking actions. Two mac files are included as examples
 *
 */
