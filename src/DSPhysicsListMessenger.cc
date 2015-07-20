#include <G4UnitsTable.hh>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSPhysicsListMessenger.hh"
#include "DSPhysicsList.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "DSIO.hh"

class DSPhysicsListMessenger;

using namespace std;


DSPhysicsListMessenger::DSPhysicsListMessenger(DSPhysicsList *phys){
  fPhysicsList = phys;
  fDirectory = new G4UIdirectory("/ds/physics/");
  fDirectory->SetGuidance("Control of physical processes");

  fHadCmd = new G4UIcmdWithAString("/ds/physics/hadronic_list",this);
  G4String hadcandidates = "none HP QGSP_BERT_HP QGSP_BIC_HP FTF_BIC_HP FTFP_BERT_HP Shielding";  
  fHadCmd->SetCandidates(hadcandidates);
  
  fEMCmd = new G4UIcmdWithAString("/ds/physics/em_list",this);
  G4String emcandidates = "standard livermore";  
  fEMCmd->SetCandidates(emcandidates);

  fHPRangeCutsCmd     = new G4UIcmdWithABool("/ds/physics/HPRangeCuts",this);

  fDepositCutsCmd     = new G4UIcmdWithABool("/ds/physics/DepositCuts",this);
  
  fOpticsCmd          = new G4UIcmdWithAnInteger("/ds/physics/optics",this);
  
  fKillS1S2Cmd        = new G4UIcmdWithABool("/ds/physics/killS1S2",this);
  
  fKillS2Cmd          = new G4UIcmdWithABool("/ds/physics/killS2",this);

  fKillS1Cmd          = new G4UIcmdWithABool("/ds/physics/killS1",this);
   
  fScaleS2Cmd         = new G4UIcmdWithADouble("/ds/physics/scaleS2",this);
  
  fTunedS1At200VCmd   = new G4UIcmdWithABool("/ds/physics/tuned200V",this);
  
  fDriftFieldCmd      = new G4UIcmdWithADoubleAndUnit("/ds/physics/DriftField",this);
  new G4UnitDefinition("V/cm",  "V/cm",  "Electric field", 100 * volt/m );
  new G4UnitDefinition("kV/cm", "kV/cm", "Electric field", 1e5 * volt/m );
  fDriftFieldCmd->SetUnitCategory("Electric field");
  fDriftFieldCmd->SetDefaultUnit("V/cm");
  fDriftFieldCmd->SetUnitCandidates("V/cm kV/cm");

  fExtractionFieldCmd = new G4UIcmdWithADoubleAndUnit("/ds/physics/ExtractionField",this);
  fExtractionFieldCmd->SetUnitCategory("Electric field");
  fExtractionFieldCmd->SetDefaultUnit("kV/cm");
  fExtractionFieldCmd->SetUnitCandidates("V/cm kV/cm");

  fThomasImelNFCmd    = new G4UIcmdWithADouble("/ds/physics/ThomasImelNF",this);

  fThomasImelEp0Cmd   = new G4UIcmdWithADouble("/ds/physics/ThomasImelEp0",this);

  fThomasImelEp1Cmd   = new G4UIcmdWithADouble("/ds/physics/ThomasImelEp1",this);

  fDokeBirksNFp1Cmd   = new G4UIcmdWithADouble("/ds/physics/DokeBirksNFp1",this);

  fDokeBirksNFp3Cmd   = new G4UIcmdWithADouble("/ds/physics/DokeBirksNFp3",this);

  fDokeBirksEp1Cmd    = new G4UIcmdWithADouble("/ds/physics/DokeBirksEp1",this);

  fDokeBirksEp2Cmd    = new G4UIcmdWithADouble("/ds/physics/DokeBirksEp2",this);

  fDokeBirksEp3Cmd    = new G4UIcmdWithADouble("/ds/physics/DokeBirksEp3",this);

}


DSPhysicsListMessenger::~DSPhysicsListMessenger() {
  delete fDirectory;
  delete fHadCmd; 
  delete fEMCmd; 
  delete fS1S2Cmd;
  delete fHPRangeCutsCmd;
  delete fOpticsCmd;
  delete fKillS2Cmd;
  delete fKillS1S2Cmd;
  delete fScaleS2Cmd;
  delete fDriftFieldCmd;
  delete fExtractionFieldCmd;
  delete fThomasImelNFCmd;
  delete fThomasImelEp0Cmd;
  delete fThomasImelEp1Cmd;
  delete fDokeBirksNFp1Cmd;
  delete fDokeBirksNFp3Cmd;
  delete fDokeBirksEp1Cmd;
  delete fDokeBirksEp2Cmd;
  delete fDokeBirksEp3Cmd;
  delete fTunedS1At200VCmd;
}


void DSPhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 

  if(command == fS1S2Cmd)    {
    fPhysicsList->SetS1S2(fS1S2Cmd->ConvertToBool(newValue));
    DSLog(routine) <<" S1S2 Physics: "<< newValue << endlog;
  } else if(command == fHadCmd)    {
    if(newValue == "none")         fPhysicsList->SetHadronicList(-1);
    if(newValue == "HP")           fPhysicsList->SetHadronicList(0);
    if(newValue == "QGSP_BERT_HP") fPhysicsList->SetHadronicList(1);
    if(newValue == "QGSP_BIC_HP")  fPhysicsList->SetHadronicList(2);
    if(newValue == "FTF_BIC_HP")   fPhysicsList->SetHadronicList(3);
    if(newValue == "FTFP_BERT_HP") fPhysicsList->SetHadronicList(4);
    if(newValue == "Shielding")    fPhysicsList->SetHadronicList(5);
    DSLog(routine) <<" Hadronic List: "<< newValue << endlog;
  } else if(command == fEMCmd)    {
    if(newValue == "standard")       fPhysicsList->SetEMList(1);
    if(newValue == "Livermore")      fPhysicsList->SetEMList(2);
    fPhysicsList->SetS1S2(fS1S2Cmd->ConvertToBool(newValue));
    DSLog(routine) <<" EM List: "<< newValue << endlog;
  } else if(command == fHPRangeCutsCmd) {
    fPhysicsList->SetHPRangeCuts( fHPRangeCutsCmd->ConvertToBool(newValue) );
    DSLog(routine) <<" High Precision Range Cuts in LAr: " << newValue << endlog;
  } else if(command == fDepositCutsCmd) {
    fPhysicsList->SetDepositCuts( fDepositCutsCmd->ConvertToBool(newValue) );
    DSLog(routine) <<" Deposit Physics Cuts:" << newValue << endlog;
  } else if(command == fOpticsCmd)    {
    fPhysicsList->SetOptics(fOpticsCmd->ConvertToInt(newValue));
    DSLog(routine) <<"Optical Processes: "<< newValue << endlog;
  } else if(command == fKillS1S2Cmd)    {
    DSStorage::Get()->SetKillS1S2(fKillS1S2Cmd->ConvertToBool(newValue));
    DSLog(routine) <<"Kill S1 and S2 Light: "<< newValue << endlog;
  } else if(command == fKillS2Cmd)    {
    DSStorage::Get()->SetKillS2(fKillS2Cmd->ConvertToBool(newValue));
    DSLog(routine) <<"Kill S2 Light: "<< newValue << endlog;
  } else if(command == fKillS1Cmd)    {
    DSStorage::Get()->SetKillS1(fKillS1Cmd->ConvertToBool(newValue));
    DSLog(routine) <<"Kill S1 Light: "<< newValue << endlog;
  } else if(command == fScaleS2Cmd)    {
    DSStorage::Get()->SetScaleS2(fScaleS2Cmd->ConvertToDouble(newValue));
    DSLog(routine) <<"Scale S2 Light by: "<< newValue << endlog;
  } else if(command == fDriftFieldCmd)    {
    DSStorage::Get()->SetDriftField(fDriftFieldCmd->ConvertToDimensionedDouble(newValue));
    DSStorage::Get()->SetTunedS1At200V(false);
    DSLog(routine)<<"Drift Field: " << newValue << endlog;
  } else if(command == fExtractionFieldCmd)    {
    DSStorage::Get()->SetExtractionField(fExtractionFieldCmd->ConvertToDimensionedDouble(newValue));
    DSLog(routine)<<"Extraction Field: " << newValue << endlog;
  } else if(command == fThomasImelNFCmd)   {
    DSStorage::Get()->SetThomasImelNullField( fThomasImelNFCmd->ConvertToDouble(newValue));
    DSLog(routine) <<"ThomasImel parameter at null field set to " << newValue << endlog;
  } else if(command == fThomasImelEp0Cmd)   {
    DSStorage::Get()->SetThomasImelEp0( fThomasImelEp0Cmd->ConvertToDouble(newValue));
    DSLog(routine) <<"ThomasImel parameter0 at non-null field set to " << newValue << endlog;
  } else if(command == fThomasImelEp1Cmd)   {
    DSStorage::Get()->SetThomasImelEp1( fThomasImelEp1Cmd->ConvertToDouble(newValue));
    DSLog(routine) <<"ThomasImel parameter1 at non-null field set to " << newValue << endlog;
  } else if(command == fDokeBirksNFp1Cmd)   {
    DSStorage::Get()->SetDokeBirksNFp1( fDokeBirksNFp1Cmd->ConvertToDouble(newValue));
    DSLog(routine) <<"DokeBirks parameter1 at null field set to " << newValue << endlog;
  } else if(command == fDokeBirksNFp3Cmd)   {
    DSStorage::Get()->SetDokeBirksNFp3( fDokeBirksNFp3Cmd->ConvertToDouble(newValue));
    DSLog(routine) <<"DokeBirks parameter3 at null field set to " << newValue << endlog;
  } else if(command == fDokeBirksEp1Cmd)   {
    DSStorage::Get()->SetDokeBirksEp1( fDokeBirksEp1Cmd->ConvertToDouble(newValue));
    DSLog(routine) <<"DokeBirks parameter1 at non-null field set to " << newValue << endlog;
  } else if(command == fDokeBirksEp2Cmd)   {
    DSStorage::Get()->SetDokeBirksEp2( fDokeBirksEp2Cmd->ConvertToDouble(newValue));
    DSLog(routine) <<"DokeBirks parameter2 at non-null field set to " << newValue << endlog;
  } else if(command == fDokeBirksEp3Cmd)   {
    DSStorage::Get()->SetDokeBirksEp3( fDokeBirksEp3Cmd->ConvertToDouble(newValue));
    DSLog(routine) <<"DokeBirks parameter3 at non-null field set to " << newValue << endlog;
  } else if(command == fTunedS1At200VCmd)   {
    DSStorage::Get()->SetTunedS1At200V(fTunedS1At200VCmd->ConvertToBool(newValue));
    DSLog(routine) <<"Tuned S1 and S2 at 200 V/cm: " << newValue << endlog;
  }
}
 
/*
 * $Log: DSPhysicsListMessenger.cc,v $
 * Revision 1.4  2015/04/29 14:49:13  dfranco
 * Fixed a bug in the DSGeneratorEnergyDeposit
 *
 * Revision 1.3  2014/12/22 14:40:43  dfranco
 * added the option to activate the recombination probability at 200 V/cm (/ds/physics/tunedS1); this option is by default true; selecting a specific drift field automatically switch off the tunedS1 option
 *
 * Revision 1.2  2014/07/23 14:52:41  pagnes
 * write thermal e- and kill S1 commands added
 *
 * Revision 1.1  2014/05/07 12:21:04  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.13  2014/01/29 13:13:41  perassos
 * Update of the electric field handling and of the Nuclear Recoils generator
 *
 * Revision 1.12  2014/01/21 10:50:49  perassos
 * Range cuts set to 1um in LAr and to 1 mm elsewhere
 *
 * Revision 1.11  2014/01/07 14:10:36  perassos
 * Added the commands to set in the macfile the electric field and the Thomas-Imel parameters
 *
 * Revision 1.10  2013/07/25 09:55:55  dfranco
 * Added a second version of NEST, still in development, and actually not working. The default version is DSLight
 *
 * Revision 1.9  2013/07/24 09:49:02  dfranco
 * Added S1 and S2 equivalent energies, the material index for each deposit, the command killS1S2 to kill photons and electrons generated by DSLight (after storing the equivalent energies)
 *
 * Revision 1.8  2013/07/16 08:52:22  dfranco
 * added high precision (HP) neutron physics list for low energies
 *
 * Revision 1.7  2013/06/10 14:15:39  dfranco
 * Added two commands: /ds/physics/killS2 and /ds/physics/scaleS2 to kill or scale the S2 light
 *
 * Revision 1.6  2013/06/03 08:35:16  dfranco
 * added command to kill the optical processes
 *
 * Revision 1.5  2013/04/18 13:28:17  dfranco
 * added process to the physics list
 *
 * Revision 1.4  2013/04/18 12:55:33  dfranco
 * Added hadronic and em physics lists
 *
 * Revision 1.3  2013/04/10 21:06:06  meregaglia
 * added S1S2 from NEST and possibility to switch on or off this physics. Work in progress on the physics itsel
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
