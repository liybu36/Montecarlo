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
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DSEventAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSEventActionMessenger.hh"
#include "DSEventHandler.hh"
#include "DSStorage.hh"
#include "DSLogger.hh"
#include "DSIO.hh"
#include "Randomize.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4String.hh"
#include "G4Timer.hh"
#include <iostream>
#include <fstream>
#include <string>

#include <vector>
#include <algorithm>

using namespace std;



DSEventAction::DSEventAction(){
  fTotNPE    = 0;  
  fTimer     = new G4Timer;
  fTimer->Start();
  fMessenger = new DSEventActionMessenger(this);
}


DSEventAction::~DSEventAction(){
  delete fMessenger;
}


void DSEventAction::BeginOfEventAction(const G4Event* evt){
 ;
}


void DSEventAction::EndOfEventAction(const G4Event* evt){

  // Save event infos  
  DSEventHandler::Get()->SetEventID(evt->GetEventID ());
  DSEventHandler::Get()->SetNPE(int(DSEventHandler::Get()->GetVPhotoElectrons().size()));
  DSEventHandler::Get()->SetMuNPE(int(DSEventHandler::Get()->GetVMuPhotoElectrons().size()));
  DSEventHandler::Get()->SetVetoNPE(int(DSEventHandler::Get()->GetVVetoPhotoElectrons().size()));
  DSEventHandler::Get()->SetNPH(int(DSEventHandler::Get()->GetVPhotons().size()));
  DSEventHandler::Get()->SetNDeposits(int(DSEventHandler::Get()->GetVDeposits().size()));
  DSEventHandler::Get()->SetNDaughters(int(DSEventHandler::Get()->GetVDaughters().size()));
  DSEventHandler::Get()->SetNUsers(int(DSEventHandler::Get()->GetVUsers().size()));


  // Write info on standard output and log file  
  int EventCounter =  DSStorage::Get()->GetEventCounter()  ;
  fTotNPE += DSEventHandler::Get()->GetNPE();
  fTotNPE += DSEventHandler::Get()->GetVetoNPE();
  fTotNPE += DSEventHandler::Get()->GetMuNPE();
  fTotNPH += DSEventHandler::Get()->GetNPH() ; 
  int OverAllNPE = DSEventHandler::Get()->GetNPE() +  DSEventHandler::Get()->GetVetoNPE() + DSEventHandler::Get()->GetMuNPE()+fTotNPH;
  

  if (!evt->GetEventID ())  { fTimer->Stop();    fTimer->Start(); }
  if ((evt->GetEventID()%EventCounter) == 0  ) {  
    fTimer->Stop();
    if(evt->GetEventID() == 0 ) {
      DSLog(routine) << ">>> Event " << evt->GetEventID()  << ";  NPE = " <<  OverAllNPE 
                     << ";  NPE/event = " <<  G4float(fTotNPE)  
                     << ";  CPUTime/event = " <<  fTimer->GetRealElapsed() << " s" << endlog ;
    } else {
      DSLog(routine) << ">>> Event " << evt->GetEventID()  << ";  NPE = " <<  OverAllNPE 
                     << ";  NPE/event = " <<  G4float(fTotNPE)/G4float(EventCounter)  
                     << ";  CPUTime/event = " <<  fTimer->GetRealElapsed()/G4float(EventCounter) << " s" << endlog ;   
    }
    DSLog(trace)   << "    Starting Position: " <<  DSEventHandler::Get()->GetPosition()/cm << " cm" << endlog;
    DSLog(trace)   << "    Energy           : " <<  DSEventHandler::Get()->GetEnergy()/MeV << " MeV" << endlog;
    DSLog(trace)   << "    PDG              : " <<  DSEventHandler::Get()->GetPDG() <<  endlog;
    if(DSStorage::Get()->GetVerbosity() > 0 && DSStorage::Get()->GetNumberOfHits() )
      DSLog(trace)   << "    Fraction of true positions:  " << float( evt->GetEventID() ) / DSStorage::Get()->GetNumberOfHits() << endlog; 
    DSLog(trace)   << endlog;
   
    fTimer->Start();
    fTotNPE = 0;
    fTotNPH =0 ; 
    // Extra info on standard output
    if(DSLogger::GetSeverity() <= DSLogger::trace) {
      if(DSStorage::Get()->GetVerbosity() > 0) DSEventHandler::Get()->DumpEvent();
      if(DSStorage::Get()->GetVerbosity() > 1) DSEventHandler::Get()->DumpDaughter();
      if(DSStorage::Get()->GetVerbosity() > 2) DSEventHandler::Get()->DumpDeposit();
      if(DSStorage::Get()->GetVerbosity() > 2) DSEventHandler::Get()->DumpUser();
      if(DSStorage::Get()->GetVerbosity() > 3) DSEventHandler::Get()->DumpPhotoElectron();
      if(DSStorage::Get()->GetVerbosity() > 3) DSEventHandler::Get()->DumpMuPhotoElectron();
      if(DSStorage::Get()->GetVerbosity() > 3) DSEventHandler::Get()->DumpVetoPhotoElectron();
      if(DSStorage::Get()->GetVerbosity() > 4) DSEventHandler::Get()->DumpPhoton();
    }  
  }



  // Write output on binary file
  int SIZE =                                          sizeof(EventStructureDiskFormat) 
             + DSEventHandler::Get()->GetNDaughters()*sizeof(DaughterStructure)
             + DSEventHandler::Get()->GetNDeposits() *sizeof(DepositStructure)
             + DSEventHandler::Get()->GetNUsers()    *sizeof(UserStructure)
             + DSEventHandler::Get()->GetNPH()       *sizeof(PhotonStructure)
             + DSEventHandler::Get()->GetNPE()       *sizeof(PhotoElectronStructure)
             + DSEventHandler::Get()->GetVetoNPE()   *sizeof(PhotoElectronStructure)
             + DSEventHandler::Get()->GetMuNPE()     *sizeof(PhotoElectronStructure);


  DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&SIZE)   ,sizeof( int  ));
  DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetEvent())                   ,sizeof(EventStructureDiskFormat));
  for(int i = 0; i <DSEventHandler::Get()->GetNDaughters(); i++ ) 
     DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVDaughters()[i])        ,sizeof(DaughterStructure));
  for(int i = 0; i <DSEventHandler::Get()->GetNDeposits(); i++ ) 
     DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVDeposits()[i])         ,sizeof(DepositStructure));       
  for(int i = 0; i <DSEventHandler::Get()->GetNUsers(); i++ ) 
    DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVUsers()[i])             ,sizeof(UserStructure));
  for(int i = 0; i <DSEventHandler::Get()->GetNPH(); i++ ) 
    DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVPhotons()[i])           ,sizeof(PhotonStructure));
  for(int i = 0; i <DSEventHandler::Get()->GetNPE(); i++ ) 
    DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVPhotoElectrons()[i])    ,sizeof(PhotoElectronStructure));    
  for(int i = 0; i <DSEventHandler::Get()->GetMuNPE(); i++ ) 
    DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVMuPhotoElectrons()[i])  ,sizeof(PhotoElectronStructure));    
  for(int i = 0; i <DSEventHandler::Get()->GetVetoNPE(); i++ ) 
    DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetVVetoPhotoElectrons()[i]),sizeof(PhotoElectronStructure));    
  DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&SIZE) ,sizeof( int    ));


  DSEventHandler::Get()->ClearAll();
  
}

/*
 * $Log: DSEventAction.cc,v $
 * Revision 1.5  2015/04/28 10:15:52  pagnes
 * SiPM photoelectrons stored in DS20k (nph structure)
 *
 * Revision 1.4  2014/11/13 16:47:04  dfranco
 * removed variables which were creating conflicts with the previous version of g4ds10
 *
 * Revision 1.3  2014/10/13 18:43:54  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.2  2014/07/17 13:24:20  dfranco
 * Added gaussian smearing in the photon generation; fixed a bug in the std output of the eventaction
 *
 * Revision 1.1  2014/05/07 12:21:02  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.7  2014/04/11 10:20:40  perassos
 * Added generation in materials
 *
 * Revision 1.6  2013/10/22 18:11:12  dfranco
 * (hopefully) fixed bug in simulating radioactive chains
 *
 * Revision 1.5  2013/08/20 03:25:52  swesterd
 * added G2 TPC geoemtry (not complete) and added monoenergetic energy distribution to generator
 *
 * Revision 1.4  2013/04/03 10:14:25  dfranco
 * Fixed bugs with RDM and RDMChain staking actions. The logic of Geant4 is changed. Different excited states of a nucleus correspond to new particles (trackID). Code adapted.
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
