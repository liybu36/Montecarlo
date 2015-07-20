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

#include "DSDetectorConstructionMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSDetectorConstruction.hh"
#include "DSMaterial.hh"
#include "DSStorage.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIdirectory.hh"
#include "DSLogger.hh"
#include "DSIO.hh"

class DSDetectorConstructionMessenger;

using namespace std;


DSDetectorConstructionMessenger::DSDetectorConstructionMessenger(DSDetectorConstruction *det){
  fDetectorConstruction = det;
  fDirectory = new G4UIdirectory("/ds/detector/");
  fDirectory->SetGuidance("Control detector gemetry and materials");

  fConfigurationCmd = new G4UIcmdWithAnInteger("/ds/detector/configuration",this);
 
  fScintillatorCmd = new G4UIcmdWithAString("/ds/detector/scintillator",this);
  G4String candidates = "BoronScintillator GdScintillator";  
  fScintillatorCmd->SetCandidates(candidates);
  
  fSourceCmd  = new G4UIcmdWith3VectorAndUnit("/ds/detector/source_position",this );
  fSourceCmd->SetUnitCandidates("mm cm m");

  fExtLArScintillatingCmd = new G4UIcmdWithABool("/ds/detector/ExtLarScintillating", this);

  fVetoYieldFactorCmd = new G4UIcmdWithADouble("/ds/detector/vetoyieldfactor", this);

  fHolderRadiusCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/holderRadius",this);
  fHolderRadiusCmd->SetUnitCandidates("mm cm m");

  fHolderZCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/holderZ",this);
  fHolderZCmd->SetUnitCandidates("mm cm m");

  fHolderPhiCmd = new G4UIcmdWithADoubleAndUnit("/ds/detector/holderPhi",this);
  fHolderPhiCmd->SetUnitCandidates("degree deg rad");

}


DSDetectorConstructionMessenger::~DSDetectorConstructionMessenger() {
  delete fDirectory;
  delete fConfigurationCmd;
  delete fScintillatorCmd;
  delete fSourceCmd;
  delete fExtLArScintillatingCmd ;
  delete fVetoYieldFactorCmd;
  delete fHolderRadiusCmd;
  delete fHolderZCmd;
  delete fHolderPhiCmd;
}


void DSDetectorConstructionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 
  if(command == fConfigurationCmd) {
    DSLog(routine) << "Detector Configuration: " << newValue << endlog ;
    fDetectorConstruction->SetDetectorConfiguration(fConfigurationCmd->ConvertToInt(newValue));
  } else if(command == fScintillatorCmd) {
    DSLog(routine) << "Choosen Scintillator: " << newValue << endlog ;
    if(newValue == "BoronScintillator") DSStorage::Get()->SetScintillator(0);
    if(newValue == "GdScintillator")    DSStorage::Get()->SetScintillator(1);    
  } else  if (command == fSourceCmd) {   
    fDetectorConstruction->SetIsSource(true);    
    DSStorage::Get()->SetSourcePosition(fSourceCmd->ConvertToDimensioned3Vector(newValue));
    DSLog(routine) << "Source position: " << newValue << endlog ;
  } else  if (command == fExtLArScintillatingCmd) {           
    DSStorage::Get()->SetIsExternalLArScintillating(fExtLArScintillatingCmd->ConvertToBool(newValue));
    DSLog(routine) << "External Liquid Argon Scintillation: " << newValue << endlog ;
  } else  if (command == fVetoYieldFactorCmd) {     
    DSStorage::Get()->SetVetoYieldFactor(fVetoYieldFactorCmd->ConvertToDouble(newValue));
    DSLog(routine) << "Veto Yield Factor: " << newValue << endlog ;
  } else  if (command == fHolderRadiusCmd) {     
    DSStorage::Get()->SetHolderRadius(fHolderRadiusCmd->ConvertToDimensionedDouble(newValue));
    DSStorage::Get()->SetHolderSource(true);
    DSLog(routine) << "Holder Radius: " << newValue << endlog ;
  } else  if (command == fHolderZCmd) {     
    DSStorage::Get()->SetHolderZ(fHolderZCmd->ConvertToDimensionedDouble(newValue));
    DSStorage::Get()->SetHolderSource(true);
    DSLog(routine) << "Holder Z: " << newValue << endlog ;
  } else  if (command == fHolderPhiCmd) {     
    DSStorage::Get()->SetHolderPhi(fHolderPhiCmd->ConvertToDimensionedDouble(newValue));
    DSStorage::Get()->SetHolderSource(true);
    DSLog(routine) << "Holder Phi: " << newValue << endlog ;
  }
}

 

/*
 * $Log: DSDetectorConstructionMessenger.cc,v $
 * Revision 1.6  2015/01/14 16:58:36  dfranco
 * added source vial and correspondent commands by Laura and Bianca; manual updated
 *
 * Revision 1.5  2014/11/21 10:19:00  dfranco
 * added a command to scale the veto scintillation yield factor and fixed the visible energy variable in the veto
 *
 * Revision 1.4  2014/11/20 15:32:06  dfranco
 * added a command to remove scintillation process from liquid argon between TPC and cryostat
 *
 * Revision 1.3  2014/11/06 17:39:44  dfranco
 * added source and fixed a bug: removed gaseous nitrogen from the neutron veto
 *
 * Revision 1.2  2014/05/07 14:27:31  dfranco
 * fixed some bugs and added GdScintillator
 *
 * Revision 1.1  2014/05/07 12:21:01  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
