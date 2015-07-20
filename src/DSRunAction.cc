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

#include "G4Timer.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSEventStructure.hh"
#include "DSEventHandler.hh"
#include "DSManager.hh"
#include "DSG4DSReader.hh"
#include "DSRunAction.hh"
#include "DSRunActionMessenger.hh"
#include "DSLogger.hh"
#include "DSIO.hh"
#include "G4ios.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include <time.h>
#include <iostream>
#include "G4String.hh"
namespace CLHEP {} 
using namespace CLHEP;
using namespace std ;

DSRunAction::DSRunAction() {
  timer = new G4Timer;
  autoSeed = true;
  DSIO::Get()->SetIsBinary(false);
  DSIO::Get()->SetIsG4DS(false);
  runMessenger = new DSRunActionMessenger(this);
}

DSRunAction::~DSRunAction() {
  delete timer;
  delete runMessenger;
}

void DSRunAction::BeginOfRunAction(const G4Run* aRun) {
  G4UImanager* UI = G4UImanager::GetUIpointer();

  
  
  timer->Start();

  

  if(autoSeed) {
    DSLog(routine) << "*******************" << endlog;
    DSLog(routine) << "*** AUTOSEED ON ***" << endlog;
    DSLog(routine) << "*******************" << endlog;
    long seeds[2];
    time_t systime = time(NULL);
    seeds[0] = (long) systime;
    seeds[1] = (long) (systime*G4UniformRand());
    DSLog(routine) << "Seed: " << seeds[1] << endlog;
    HepRandom::setTheSeed(seeds[1]);
  } else {
    DSLog(routine) << "********************" << endlog;
    DSLog(routine) << "*** AUTOSEED OFF ***" << endlog;
    DSLog(routine) << "********************" << endlog; 
  
  }
  

  if(DSIO::Get()->IsG4DS()) {
    DSIO::Get()->OpenG4DSFile();
    if( DSIO::Get()->GetG4DSFile().fail()) {
      DSLog(error) << "G4DS file does not exist!" << endlog ;
      DSLog(fatal) << "Fatal " << endlog ;
    }
    DSLog(routine) << "G4DS file " << DSIO::Get()->GetG4DSFileName() << " opened" << endlog ;
    DSG4DSReader::Get()->ReadHeader();
  }

  
  DSLog(routine) << "Random seed: " << HepRandom::getTheSeed() << endlog ;
  
  //if(DSIO::Get()->IsG4DS()) {
  //  DSIO::Get()->OpenG4DSFile();
  //  if( DSIO::Get()->GetG4DSFile().fail()) {
  //    DSLog(error) << "G4DS file does not exist!" << endlog ;
  //    DSLog(fatal) << "Fatal " << endlog ;
  //  }
  //  DSLog(routine) << "G4DS file " << DSIO::Get()->GetG4DSFileName() << " opened" << endlog ;
  //}
  DSIO::Get()->OpenBinaryFile();
  
  DSLog(routine) << "Initialized Binary File: " << DSIO::Get()->GetBinaryFileName() << endlog;      
  //if(!DSIO::Get()->GetIsBinary()) { // Header writing
    int SIZE =  sizeof(HeaderStructure);
    DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&SIZE)   ,sizeof( int  ));
    DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&DSEventHandler::Get()->GetHeader()),SIZE);
    DSIO::Get()->GetBinaryFile().write(reinterpret_cast<char*>(&SIZE)   ,sizeof( int  ));
  //}
  


}
    
void DSRunAction::EndOfRunAction(const G4Run* aRun) {


  if (G4VVisManager::GetConcreteInstance()) 
        G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");

  DSIO::Get()->CloseBinaryFile();
  DSIO::Get()->CloseLogFiles();  
  DSLog(routine) << "Binary File: " << DSIO::Get()->GetBinaryFileName()<< " closed" << endlog;      
 
  timer->Stop();
  DSLog(routine) << "Number of event = " << aRun->GetNumberOfEvent() << " " << *timer << endlog;
}

/*
 * $Log: DSRunAction.cc,v $
 * Revision 1.1  2014/05/07 12:21:04  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.3  2014/03/11 09:54:38  meregaglia
 * Added generator starting from energy deposits
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
