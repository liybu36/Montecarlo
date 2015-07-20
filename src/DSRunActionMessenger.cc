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
// ********************************************************************
//
// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the DS mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// * DSRunActionMessenger class                                     *
// *                                                                  *
// ********************************************************************
//
// HISTORY
// 22/02/2004: migrated from DS-V04
//
// ********************************************************************

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "CLHEP/Random/Random.h"

#include "DSRunActionMessenger.hh"
#include "DSRunAction.hh"
#include "DSLogger.hh"
#include "DSIO.hh"
#include "DSEventAction.hh"

//class  DSIO;
//class  DSRoot;
namespace CLHEP {} 
using namespace CLHEP;

DSRunActionMessenger::DSRunActionMessenger
(DSRunAction* runAct) : runAction(runAct){
  fDirectory = new G4UIdirectory("/run/");

  fSetAutoSeedCmd = new G4UIcmdWithABool("/run/autoSeed",this);
  fSetAutoSeedCmd->SetGuidance("Switch on/off time-based random seeds");
  fSetAutoSeedCmd->SetGuidance("true: run seeds determined by system time");
  fSetAutoSeedCmd->SetGuidance("false: use command 'random/resetEngineFrom'");
  fSetAutoSeedCmd->SetGuidance("Default = true");
  //fSetAutoSeedCmd->SetParameterName("autoSeed", false);
  //fSetAutoSeedCmd->AvailableForStates(G4State_Idle);

  fHEPRandomSeedCmd = new G4UIcmdWithAnInteger("/run/heprandomseed",this);
  fHEPRandomSeedCmd->SetGuidance("Sets random number generator seed.");


  // /run/rootfilename
  fSetFileNameCmd = new G4UIcmdWithAString("/run/filename", this);
  fSetFileNameCmd->SetGuidance("Name for output files");
 
  fSetIsBinaryCmd = new G4UIcmdWithABool("/run/setbinary", this);
  fSetIsBinaryCmd->SetGuidance("Set Binary");


  fSetG4DSNameCmd = new G4UIcmdWithAString("/run/g4ds", this);
  fSetG4DSNameCmd->SetGuidance("Name of the  output files");

  fSetFileNameCmd->SetGuidance("Decay nuclide");

  fRunCmd = new G4UIcmdWithAnInteger("/run/id",this);
  
  fRateCmd = new G4UIcmdWithADoubleAndUnit("/run/rate",this);
  fRateCmd->SetUnitCandidates("Hz kHz MHz");
  
  //fLiveTimeCmd = new G4UIcmdWithADoubleAndUnit("/run/livetime",this);
  //fLiveTimeCmd->SetUnitCandidates("ms s");

  fWriteDepositsCmd = new G4UIcmdWithABool("/run/writedeposits", this);

  fWriteIsotopesCmd = new G4UIcmdWithABool("/run/writeisotopes", this);

  fWriteEBCmd = new G4UIcmdWithABool("/run/writeeb", this);

  fSetCommentCmd = new G4UIcmdWithAString("/run/comment", this);

}


DSRunActionMessenger::~DSRunActionMessenger() {

  delete fSetAutoSeedCmd;  
  delete fHEPRandomSeedCmd;
  delete fSetFileNameCmd;
  delete fSetIsBinaryCmd;
  delete fDirectory;
  delete fRunCmd;
  delete fRateCmd;
  delete fWriteDepositsCmd;
  delete fSetCommentCmd;
  delete fSetG4DSNameCmd;
  delete fWriteIsotopesCmd;
  delete fWriteEBCmd;
}


void DSRunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 

  if(command == fSetAutoSeedCmd) {
    runAction->SetAutoSeed(fSetAutoSeedCmd->ConvertToBool(newValue));    
  } else if (command == fHEPRandomSeedCmd) {
    runAction->SetAutoSeed(false);
    HepRandom::setTheSeed(fHEPRandomSeedCmd->GetNewIntValue(newValue));
    DSLog(trace) << "HepRandom seed set to: "<< fHEPRandomSeedCmd->GetNewIntValue(newValue) << endlog;
  } else if(command == fSetFileNameCmd) {      
    ; // Already set in g4ds.cc
  } else if(command == fSetIsBinaryCmd) {
    DSIO::Get()->SetIsBinary(fSetIsBinaryCmd->ConvertToBool(newValue));
    //if(!fSetIsBinaryCmd->ConvertToBool(newValue)) 
    // DSLog(routine) <<"The binary file is disabled"  << endlog;  
//  } else if(command == fSetG4DSNameCmd) {      
//      DSIO::Get()->SetG4DSFileName(newValue);
//      DSIO::Get()->SetIsG4DS(true);  
  } else if(command == fRunCmd) { 
      //DSOutputVertex::Get()->SetRun(fRunCmd->ConvertToInt(newValue));
  } else if(command == fSetCommentCmd) { 
      //DSOutputVertex::Get()->SetComment(newValue);
  } else if(command == fRateCmd) { 
      //DSOutputVertex::Get()->SetRate(fRateCmd->ConvertToDouble(newValue));
  //} else if(command == fLiveTimeCmd) { 
  //    DSOutputVertex::Get()->SetLiveTime(fLiveTimeCmd->ConvertToDouble(newValue));
  } else if(command == fWriteDepositsCmd) {
      DSLog(routine) << "Write Deposits: "<< newValue << endlog;
      //DSOutputVertex::Get()->SetWriteDeposits(fWriteDepositsCmd->ConvertToBool(newValue));    
  } else if(command == fWriteIsotopesCmd) {
      DSLog(routine) << "Write Isotopes: "<< newValue << endlog;
      //DSOutputVertex::Get()->SetWriteIsotope(fWriteIsotopesCmd->ConvertToBool(newValue));    
  } else if(command == fWriteEBCmd) {
      DSLog(routine) << "Write External Background: "<< newValue << endlog;
      //DSOutputVertex::Get()->SetWriteEB(fWriteEBCmd->ConvertToBool(newValue));    
  }
}
/*
 * $Log: DSRunActionMessenger.cc,v $
 * Revision 1.1  2014/05/07 12:21:05  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
