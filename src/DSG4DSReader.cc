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
//
//	DSG4DSReader.cc
//
#include "DSG4DSReader.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSIO.hh"
#include "DSEventHandler.hh"
#include "DSLogger.hh"
#include "DSEventStructure.hh"
//#include "DSOutputVertex.hh"
using namespace std;


DSG4DSReader* DSG4DSReader::me = 0;

// singleton
DSG4DSReader::DSG4DSReader() {
}

DSG4DSReader* DSG4DSReader::Get() {
  if (!me)
    me = new DSG4DSReader();
  return me;
}


G4bool DSG4DSReader::ReadEvent() { 
 
  int BuffDimension;
  DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&BuffDimension),  sizeof ( int ));

  
  DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fEvent),  sizeof ( EventStructureDiskFormat ));  

  DSEventHandler::Get()->SetEventID(fEvent.EventID);
  DSEventHandler::Get()->SetEnergy(fEvent.Energy);
  DSEventHandler::Get()->SetPDG(fEvent.PDG);

  /*
  DSOutputVertex::Get()->SetVertexStructureDiskFormat(fVertex);

  DSOutputVertex::Get()->SetEventID(fVertex.EventID) ;		 
  DSOutputVertex::Get()->SetPDG(fVertex.PDG)	 ; 
  DSOutputVertex::Get()->SetNSequence(fVertex.NSequence) ;  	 
  DSOutputVertex::Get()->SetIsotopeCoinc(fVertex.IsotopeCoinc) ;	
  DSOutputVertex::Get()->SetTime(fVertex.Time)	 ;	 
  DSOutputVertex::Get()->SetEnergy(fVertex.Energy)	 ;	 
  DSOutputVertex::Get()->SetVisEnergy(0) ;	
  DSOutputVertex::Get()->SetPosition(DSOutputVertex::Get()->CopyArrayToVector(fVertex.Position)); 
  DSOutputVertex::Get()->SetBarycenter(DSOutputVertex::Get()->CopyArrayToVector(fVertex.Barycenter) ); 
  DSOutputVertex::Get()->SetDirection(DSOutputVertex::Get()->CopyArrayToVector(fVertex.Direction ) ) ;
  DSOutputVertex::Get()->SetNDaughters(fVertex.NDaughters) ; 	
  DSOutputVertex::Get()->SetNDeposits(0)   ;	 
  DSOutputVertex::Get()->SetNPE(0) ;	 
  DSOutputVertex::Get()->SetMuNPE(0) ;		 
  DSOutputVertex::Get()->SetNPhotons(0) ;		 
  DSOutputVertex::Get()->SetNMuPhotons(0) ; 	
  DSOutputVertex::Get()->SetNUsers(fVertex.NUsers) ;		 
  */
 
  for(int i=0; i<fEvent.NDaughters;i++) {
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fDaughter),  sizeof ( DaughterStructure ));
    //DSOutputVertex::Get()->SetDaughter(fDaughter);  
    //DSOutputVertex::Get()->SetDaughters();  
  }
 
  for(int i = 0; i <fEvent.NDeposits; i++ ) { 
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fDeposit),  sizeof ( DepositStructure ));
    fDeposit.Energy *= MeV ;
    SetDeposit(fDeposit);   
    SetDeposits();   
    //cout<<fDeposit.Energy*MeV<<" "<<fDeposit.Time<<endl; 
    //cout<<" size "<<theDeposits.size()<<endl;
  }
  
  for(int i = 0; i <fEvent.NUsers; i++ ) { 
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fUser),  sizeof ( UserStructure ));
    //DSOutputVertex::Get()->SetUser(fUser);   
    //DSOutputVertex::Get()->SetUsers();   
  }

  for(int i = 0; i <fEvent.NPH; i++ ) { 
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fPhoton),  sizeof ( PhotonStructure ));
    //DSOutputVertex::Get()->SetPhoton(fPhoton);   
    //DSOutputVertex::Get()->SetPhotons();   
  }

  for(int i = 0; i <fEvent.NPE; i++ ) { 
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fPhotonPM),  sizeof ( PhotoElectronStructure ));
    //DSOutputVertex::Get()->SetPhoton(fPhoton);   
    //DSOutputVertex::Get()->SetPhotons();   
  }

 for(int i = 0; i <fEvent.VetoNPE; i++ ) { 
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fVetoPhoton),  sizeof ( PhotoElectronStructure ));
    //DSOutputVertex::Get()->SetMuPhoton(fMuPhoton);   
    //DSOutputVertex::Get()->SetMuPhotons();   
  }

  for(int i = 0; i <fEvent.MuNPE; i++ ) { 
    DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&fMuPhoton),  sizeof ( PhotoElectronStructure ));
    //DSOutputVertex::Get()->SetMuPhoton(fMuPhoton);   
    //DSOutputVertex::Get()->SetMuPhotons();   
  }
  
  int BuffDimension2;
  DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&BuffDimension2),  sizeof ( int ));

 

  if(BuffDimension == BuffDimension2) return true;
  return false;
}
void DSG4DSReader::ReadHeader() {

  int BuffDimension;
  DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&BuffDimension),  sizeof ( int ));
  DSIO::Get()->GetG4DSFile().ignore(BuffDimension);
  int BuffDimension2;
  DSIO::Get()->GetG4DSFile().read(reinterpret_cast<char*>(&BuffDimension2),  sizeof ( int ));
  if(BuffDimension != BuffDimension2)
    DSLog(error) << "Error in reading the G4DS file header!" << endl;

}


void DSG4DSReader::DumpHeader() {


}


void DSG4DSReader::DumpEvent() {


}



void DSG4DSReader::ClearAll(){

  theDeposits.clear();
}
/*
 * $Log: DSG4DSReader.cc,v $
 * Revision 1.2  2015/04/29 14:49:13  dfranco
 * Fixed a bug in the DSGeneratorEnergyDeposit
 *
 * Revision 1.1  2014/05/07 12:21:02  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.1  2014/03/11 09:54:38  meregaglia
 * Added generator starting from energy deposits
 *
 * Revision 1.1  2009-10-20 12:37:37  dfranco
 * Added a new generator in order to read a g4ds output file and begin a simulation
 * from the energy deposits. Useful for external background.
 *
 *
 */
