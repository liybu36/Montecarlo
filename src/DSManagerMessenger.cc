#include "G4UIcmdWithAString.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIdirectory.hh"
#include "DSLogger.hh"
#include "DSManagerMessenger.hh"
#include "DSManager.hh"
#include "DSStorage.hh"
#include "DSIO.hh"


DSManagerMessenger::DSManagerMessenger(DSManager *Manager){
  fDirectory = new G4UIdirectory("/ds/detector/");
  
  fDSLogCmd = new G4UIcmdWithAString("/ds/manager/log", this);
  fDSLogCmd->SetGuidance("Set severity of logs to report to stdout.");
  fDSLogCmd->SetGuidance("Options, in ascending order of severity, are:");
  fDSLogCmd->SetGuidance("debugging: Displays all logs ");
  fDSLogCmd->SetGuidance("trace: All logs, except debugging(default)");
  fDSLogCmd->SetGuidance("routine: All logs, except debugging and trace");
  fDSLogCmd->SetGuidance("warning: All logs, except trace, debugging and routine:");
  fDSLogCmd->SetGuidance("error: Only error and fatal logs.");
  fDSLogCmd->SetGuidance("fatal: Only fatal logs.");


  fDSOpticsCmd    = new G4UIcmdWithAString("/ds/manager/optics_tuning", this);
  fTMBfractionCmd = new G4UIcmdWithADouble("/ds/manager/TMB_fraction", this);
  fTimeCutCmd     = new G4UIcmdWithADoubleAndUnit("/ds/manager/time_cut", this);

  fOverLapCmd = new G4UIcmdWithABool("/ds/manager/checkoverlap", this);

  fGDMLCmd = new G4UIcmdWithABool("/ds/manager/GDML",this);
  fGDMLCmd->SetGuidance("Switch ON or OFF GDML export geometry");

  fWritePhotonsCmd = new G4UIcmdWithABool("/ds/manager/writephotons",this);
  fWriteDaughtersCmd = new G4UIcmdWithABool("/ds/manager/writedaughters",this);
  fWriteDepositsCmd = new G4UIcmdWithABool("/ds/manager/writedeposits",this);
  fWriteThermalElectronsCmd = new G4UIcmdWithABool("/ds/manager/writethermalelectrons",this);
  
  fEventCounterCmd = new G4UIcmdWithAnInteger("/ds/manager/eventcounter",this);
  fVerbosityCmd = new G4UIcmdWithAnInteger("/ds/manager/verbosity",this);

  fFastSimulationCmd = new G4UIcmdWithAnInteger("/ds/manager/fast_simulation",this);
  
  fDaughterDepthCmd= new G4UIcmdWithAnInteger("/ds/manager/daughterdepth",this);
}

DSManagerMessenger::~DSManagerMessenger() {

  delete fDSLogCmd;
  delete fDSOpticsCmd;
  delete fOverLapCmd;
  delete fGDMLCmd;
  delete fWritePhotonsCmd;
  delete fWriteDepositsCmd;
  delete fWriteThermalElectronsCmd;
  delete fWriteDaughtersCmd;
  delete fEventCounterCmd;
  delete fVerbosityCmd;
  delete fDaughterDepthCmd;
  delete fTMBfractionCmd;
  delete fTimeCutCmd;}


void DSManagerMessenger::SetNewValue
(G4UIcommand* command, G4String newValue) { 
  if(command == fDSLogCmd) {
    if(newValue == "debugging")
      DSLogger::SetSeverity(DSLogger::debugging);
    if(newValue == "development")
      DSLogger::SetSeverity(DSLogger::development);
    else if(newValue == "trace")
      DSLogger::SetSeverity(DSLogger::trace);
    else if(newValue == "routine")
      DSLogger::SetSeverity(DSLogger::routine);
    else if(newValue == "warning")
      DSLogger::SetSeverity(DSLogger::warning);
    else if(newValue == "error")
      DSLogger::SetSeverity(DSLogger::error);
    else if(newValue == "fatal")
      DSLogger::SetSeverity(DSLogger::fatal);
    else
      DSLog(error) << "Unknown option." << endlog;

    G4cout << "     ####################################" <<  endl ;
    G4cout << "     ######## Logger Severity: " << DSLogger::GetSeverity() << " ########" << endl ;
    G4cout << "     ####################################" <<  endl ;

  } else if(command == fOverLapCmd) {     
    DSStorage::Get()->SetCheckOverlap(fOverLapCmd->ConvertToBool(newValue));
    DSLog(routine) << " Check Volume Overlaps: " << newValue << endlog;
  } else if(command == fGDMLCmd) {
    DSStorage::Get()->SetExportGDML(fGDMLCmd->ConvertToBool(newValue));
    DSLog(routine) << " GDML export geometry: " << newValue << endlog;
  } else if(command == fDSOpticsCmd) {
    DSIO::Get()->SetDSOpticsFileName(newValue);
    DSLog(routine) << " Optics filename: " << newValue << endlog;
  } else if(command == fWritePhotonsCmd) {
    DSStorage::Get()->SetWritePhotons(fWritePhotonsCmd->ConvertToBool(newValue));
    DSLog(routine) << " Write photons: " << newValue << endlog;
  } else if(command == fWriteDepositsCmd) {
    DSStorage::Get()->SetWriteDeposits(fWriteDepositsCmd->ConvertToBool(newValue));
    DSLog(routine) << " Write deposits: " << newValue << endlog;
  } else if(command == fWriteThermalElectronsCmd) {
    DSStorage::Get()->SetWriteThermalElectrons(fWriteThermalElectronsCmd->ConvertToBool(newValue));
    DSLog(routine) << " Write Thermal Electrons: " << newValue << endlog;
  } else if(command == fWriteDaughtersCmd) {
    DSStorage::Get()->SetWriteDaughters(fWriteDaughtersCmd->ConvertToBool(newValue));
    DSLog(routine) << " Write daughters: " << newValue << endlog;
  } else if(command == fEventCounterCmd) {
    DSStorage::Get()->SetEventCounter(fEventCounterCmd->ConvertToInt(newValue));
    DSLog(routine) << " Event Counter: " << newValue << endlog;
  } else if(command == fFastSimulationCmd) {
    DSStorage::Get()->SetFastSimulation(fFastSimulationCmd->ConvertToInt(newValue));
    DSLog(routine) << " Fast simulation: " << newValue << endlog;
  } else if(command == fVerbosityCmd) {
    DSStorage::Get()->SetVerbosity(fVerbosityCmd->ConvertToInt(newValue));
    DSLog(routine) << " Trace Verbosity: " << newValue << endlog;
  } else if(command == fDaughterDepthCmd) {
    DSStorage::Get()->SetNDaughters(fDaughterDepthCmd->ConvertToInt(newValue));
    DSLog(routine) << " Trace Verbosity: " << newValue << endlog;
  } 
    else if(command == fTMBfractionCmd) {
    DSStorage::Get()->SetTMBfraction(fTMBfractionCmd->ConvertToDouble(newValue));
  } 
    else if(command == fTimeCutCmd) {
    DSStorage::Get()->SetTimeCut(fTimeCutCmd->ConvertToDimensionedDouble(newValue));
  } 


 }
/*
 * $Log: DSManagerMessenger.cc,v $
 * Revision 1.3  2014/07/23 14:52:41  pagnes
 * write thermal e- and kill S1 commands added
 *
 * Revision 1.2  2014/07/16 08:23:04  pagnes
 * QE scaling to 1.0 added (/ds/manager/fast_simulation xxx)
 *
 * Revision 1.1  2014/05/07 12:21:04  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.5  2014/04/11 12:33:30  pagnes
 * command to set TMB/PC ratio inside veto added
 *
 * Revision 1.4  2014/03/19 16:37:30  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.3  2013/03/22 16:24:40  dfranco
 * added a command to set the "daughter depth" level to store in the output file
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
