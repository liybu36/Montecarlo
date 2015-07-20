#include "DSEventHandler.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

DSEventHandler* DSEventHandler::me = 0;

// singleton
DSEventHandler::DSEventHandler(){
  ClearAll();
}

DSEventHandler* DSEventHandler::Get() {
  if (!me) 
    me = new DSEventHandler();
       
  return me;
}



G4ThreeVector DSEventHandler::CopyArrayToVector( float *val) {
   G4ThreeVector vector;
   vector[0] = val[0] ;
   vector[1] = val[1] ;
   vector[2] = val[2] ;
  return vector;

}



void  DSEventHandler::SetPosition(G4ThreeVector  val)	   { 
  for(int i = 0; i < 3; i++) theEventStructure.Position[i]    = val[i]; 
}

void  DSEventHandler::SetCenterOfMass(G4ThreeVector val)     { 
  for(int i = 0; i < 3; i++) theEventStructure.CenterOfMass[i]    = val[i]; 

}
void  DSEventHandler::SetDirection(G4ThreeVector val)	   { 
  for(int i = 0; i < 3; i++) theEventStructure.Direction[i]    = val[i]; 
}

void  DSEventHandler::SetDPosition(G4ThreeVector val)       { 
  for(int i = 0; i < 3; i++) theDaughterStructure.Position[i] = val[i]  ; 
}

void  DSEventHandler::SetDDirection(G4ThreeVector val)      { 
  for(int i = 0; i < 3; i++) theDaughterStructure.Direction[i] = val[i]  ; 
}

void  DSEventHandler::SetDepPosition(G4ThreeVector val) { 
  for(int i = 0; i < 3; i++) theDepositStructure.Position[i]     = val[i]; 
}

void  DSEventHandler::SetPhotonPosition(G4ThreeVector val) { 
  for(int i = 0; i < 3; i++) thePhotonStructure.Position[i]     = val[i]; 
}


void DSEventHandler::DumpHeader() {    
  DSLog(trace) << "Header: N    " << theHeaderStructure.Events	        << endlog;
  DSLog(trace) << "Header: Run  " << theHeaderStructure.Run 	        << endlog;
  DSLog(trace) << "Header: Rate " << theHeaderStructure.Rate	        << endlog; 
  DSLog(trace) << "Header: Det  " << theHeaderStructure.DetectorFlag    << endlog;   
  DSLog(trace) << "Header: Gen  " << theHeaderStructure.LArIndex	<< endlog;   
  DSLog(trace) << "Header: Gen  " << theHeaderStructure.ScintillatorIndex<< endlog;   
  DSLog(trace) << "Header: PY   " << theHeaderStructure.PDG             << endlog;   
  DSLog(trace) << endlog;      
}
void DSEventHandler::DumpEvent(){
  DSLog(trace) << "Event: ID        " << theEventStructure.EventID	    << endlog;
  DSLog(trace) << "Event: PDG       " << theEventStructure.PDG		    << endlog;     
  DSLog(trace) << "Event: T         " << theEventStructure.Time		    << endlog;     
  DSLog(trace) << "Event: E         " << theEventStructure.Energy 	    << endlog; 
  DSLog(trace) << "Event: S1 E      " << theEventStructure.S1Energy	    << endlog; 
  DSLog(trace) << "Event: S2 E      " << theEventStructure.S2Energy	    << endlog; 
  DSLog(trace) << "Event: Veto QE   " << theEventStructure.VetoVisEnergy    << endlog; 
  DSLog(trace) << "Event: TPC QE    " << theEventStructure.MuVisEnergy	    << endlog; 
  DSLog(trace) << "Event: TPC DepE  " << theEventStructure.TPCDepEnergy	    << endlog; 
  DSLog(trace) << "Event: Veto DepE " << theEventStructure.VetoDepEnergy    << endlog; 
  DSLog(trace) << "Event: Mu DepE   " << theEventStructure.MuDepEnergy	    << endlog; 
  DSLog(trace) << "Event: X         " << theEventStructure.Position[0]	    << endlog;         
  DSLog(trace) << "Event: Y         " << theEventStructure.Position[1]	    << endlog;         
  DSLog(trace) << "Event: Z         " << theEventStructure.Position[2]	    << endlog;         
  DSLog(trace) << "Event: BX        " << theEventStructure.CenterOfMass[0]    << endlog;       
  DSLog(trace) << "Event: BY        " << theEventStructure.CenterOfMass[1]    << endlog;       
  DSLog(trace) << "Event: BZ        " << theEventStructure.CenterOfMass[2]    << endlog;       
  DSLog(trace) << "Event: pX        " << theEventStructure.Direction[0]	    << endlog;         
  DSLog(trace) << "Event: pY        " << theEventStructure.Direction[1]	    << endlog;         
  DSLog(trace) << "Event: pZ        " << theEventStructure.Direction[2]	    << endlog;         
  DSLog(trace) << "Event: ND        " << theEventStructure.NDaughters	    << endlog;         
  DSLog(trace) << "Event: NDep      " << theEventStructure.NDeposits	    << endlog;         
  DSLog(trace) << "Event: NUs       " << theEventStructure.NUsers 	    << endlog;         
  DSLog(trace) << "Event: NPE       " << theEventStructure.NPE		    << endlog;  	       
  DSLog(trace) << "Event: MuNPE     " << theEventStructure.MuNPE  	    << endlog;		   
  DSLog(trace) << "Event: VetoNPE   " << theEventStructure.VetoNPE  	    << endlog;		   
  DSLog(trace) << "Event: NPH       " << theEventStructure.NPH	            << endlog;		   
  DSLog(trace) << " " << endlog;      
}

void DSEventHandler::DumpPhotoElectron(){
  for(int i = 0; i < int(theEventStructure.thePhotoElectrons.size()); i++) {
    DSLog(trace) << "PhotoElectron: PMT      " << theEventStructure.thePhotoElectrons[i].PMT      << endlog;
    DSLog(trace) << "PhotoElectron: Time     " << theEventStructure.thePhotoElectrons[i].Time     << endlog;    
  }
}

void DSEventHandler::DumpMuPhotoElectron(){
  for(int i = 0; i < int(theEventStructure.theMuPhotoElectrons.size()); i++) {
    DSLog(trace) << "MuPhotoElectron: PMT      " << theEventStructure.theMuPhotoElectrons[i].PMT      << endlog;
    DSLog(trace) << "MuPhotoElectron: Time     " << theEventStructure.theMuPhotoElectrons[i].Time     << endlog;    
  }
}

void DSEventHandler::DumpVetoPhotoElectron(){
  for(int i = 0; i < int(theEventStructure.theVetoPhotoElectrons.size()); i++) {
    DSLog(trace) << "VetoPhotoElectron: PMT      " << theEventStructure.theVetoPhotoElectrons[i].PMT      << endlog;
    DSLog(trace) << "VetoPhotoElectron: Time     " << theEventStructure.theVetoPhotoElectrons[i].Time     << endlog;    
  }
}


void DSEventHandler::DumpPhoton(){
  for(int i = 0; i < int(theEventStructure.thePhotons.size()); i++) {
    DSLog(trace) << "Photon: VolumeID   " << theEventStructure.thePhotons[i].VolumeID       << endlog;
    DSLog(trace) << "Photon: PID        " << theEventStructure.thePhotons[i].PID            << endlog;
    DSLog(trace) << "Photon: Time       " << theEventStructure.thePhotons[i].Time           << endlog;    
    DSLog(trace) << "Photon: Wavelength " << theEventStructure.thePhotons[i].Wavelength     << endlog;    
    for(int j=0;j<3;j++)
      DSLog(trace) << "Photon: Position      " << j << " " << theEventStructure.thePhotons[i].Position[j]      << endlog;
    DSLog(trace) << endlog;     
  }
}



void DSEventHandler::DumpDaughter(){
  for(int i = 0; i < int(theEventStructure.theDaughters.size()); i++) {
    DSLog(trace) << "Daughter: Id   " << theEventStructure.theDaughters[i].Id           << endlog;
    DSLog(trace) << "Daughter: PDG  " << theEventStructure.theDaughters[i].PDG          << endlog;     
    // DSLog(trace) << "Daughter: PID  " << theEventStructure.theDaughters[i].PID          << endlog;  
    DSLog(trace) << "Daughter: TrackID  " << theEventStructure.theDaughters[i].TrackID          << endlog;  
    DSLog(trace) << "Daughter: ParentTrackID  " << theEventStructure.theDaughters[i].ParentTrackID          << endlog;     
    DSLog(trace) << "Daughter: PROC " << theEventStructure.theDaughters[i].Process      << endlog;  
    DSLog(trace) << "Daughter: T    " << theEventStructure.theDaughters[i].Time	        << endlog;  
    DSLog(trace) << "Daughter: E    " << theEventStructure.theDaughters[i].Energy       << endlog; 
    DSLog(trace) << "Daughter: X    " << theEventStructure.theDaughters[i].Position[0]  << endlog;     
    DSLog(trace) << "Daughter: Y    " << theEventStructure.theDaughters[i].Position[1]  << endlog;     
    DSLog(trace) << "Daughter: Z    " << theEventStructure.theDaughters[i].Position[2]  << endlog;     
    DSLog(trace) << "Daughter: pX   " << theEventStructure.theDaughters[i].Direction[0] << endlog;     
    DSLog(trace) << "Daughter: pY   " << theEventStructure.theDaughters[i].Direction[1] << endlog;     
    DSLog(trace) << "Daughter: pZ   " << theEventStructure.theDaughters[i].Direction[2] << endlog;     
    DSLog(trace) << endlog;     
  }
}

void DSEventHandler::DumpDeposit(){
  for(int i = 0; i < int(theEventStructure.theDeposits.size()); i++) {
    DSLog(trace) << "Deposit: PDG "    << theEventStructure.theDeposits[i].PID           << endlog;
    DSLog(trace) << "Deposit: E   "    << theEventStructure.theDeposits[i].Energy	      << endlog;    
    DSLog(trace) << "Deposit: KineticE   "    << theEventStructure.theDeposits[i].KineticEnergy	      << endlog;    
    DSLog(trace) << "Deposit: TotalE   "    << theEventStructure.theDeposits[i].TotalEnergy	      << endlog;    
    DSLog(trace) << "Deposit: V   "    << theEventStructure.theDeposits[i].Volume	      << endlog;    
    DSLog(trace) << "Deposit: Track   "    << theEventStructure.theDeposits[i].Track	      << endlog;    
    DSLog(trace) << "Deposit: ParentTrack   "    << theEventStructure.theDeposits[i].ParentTrack	      << endlog;    
    DSLog(trace) << "Deposit: Step   " << theEventStructure.theDeposits[i].Step	      << endlog;    
    DSLog(trace) << "Deposit: T   "    << theEventStructure.theDeposits[i].Time	      << endlog;    
    DSLog(trace) << "Deposit: X   "    << theEventStructure.theDeposits[i].Position[0]   << endlog;	
    DSLog(trace) << "Deposit: Y   "    << theEventStructure.theDeposits[i].Position[1]   << endlog; 
    DSLog(trace) << "Deposit: Z   "    << theEventStructure.theDeposits[i].Position[2]   << endlog;	
    DSLog(trace) << endlog;     
  }
}

void DSEventHandler::DumpUser(){
  for(int i = 0; i < int(theEventStructure.theUsers.size()); i++) {
    DSLog(trace) << "User: I1  " << theEventStructure.theUsers[i].UserInt1    << endlog;
    DSLog(trace) << "User: I2  " << theEventStructure.theUsers[i].UserInt2    << endlog;  
    DSLog(trace) << "User: F1  " << theEventStructure.theUsers[i].UserFloat1  << endlog;	 
    DSLog(trace) << "User: F2  " << theEventStructure.theUsers[i].UserFloat2  << endlog; 
    DSLog(trace) << "User: D   " << theEventStructure.theUsers[i].UserDouble  << endlog;	 
    DSLog(trace) << endlog;     
  }
}

void DSEventHandler::DumpAll(){
  DumpHeader();
  DumpEvent();
  DumpPhotoElectron();
  DumpMuPhotoElectron();
  DumpVetoPhotoElectron();
  DumpPhoton();
  DumpDaughter();
  DumpDeposit();
  DumpUser();
}

void DSEventHandler::ClearHeader() {
  theHeaderStructure.Events	     = 0;
  theHeaderStructure.Run	     = 0;
  theHeaderStructure.Rate	     = 0;
  theHeaderStructure.DetectorFlag    = 0;
  theHeaderStructure.LArIndex	     = -1;
  theHeaderStructure.ScintillatorIndex= -1;
  theHeaderStructure.PDG             = 0;
}

void DSEventHandler::ClearEvent()  {
  theEventStructure.EventID	      = 0; 
  theEventStructure.PDG	              = 0;
  theEventStructure.Time	      = 0;
  theEventStructure.Energy	      = 0;
  theEventStructure.S1Energy          = 0;
  theEventStructure.S2Energy          = 0;
  theEventStructure.VetoVisEnergy     = 0;
  theEventStructure.MuVisEnergy	      = 0;
  theEventStructure.TPCDepEnergy      = 0;
  theEventStructure.VetoDepEnergy     = 0;
  theEventStructure.MuDepEnergy       = 0;
  theEventStructure.Position[0]       = 0;
  theEventStructure.Position[1]       = 0;
  theEventStructure.Position[2]       = 0;
  theEventStructure.CenterOfMass[0]   = 0;
  theEventStructure.CenterOfMass[1]   = 0;
  theEventStructure.CenterOfMass[2]   = 0;
  theEventStructure.Direction[0]      = 0;
  theEventStructure.Direction[1]      = 0;
  theEventStructure.Direction[2]      = 0;
  theEventStructure.NDaughters        = 0;
  theEventStructure.NDeposits         = 0;
  theEventStructure.NPE	              = 0;
  theEventStructure.MuNPE	      = 0;
  theEventStructure.VetoNPE	      = 0;
  theEventStructure.NPH	              = 0;
  (theEventStructure.thePhotoElectrons).clear();    
  (theEventStructure.theMuPhotoElectrons).clear();    
  (theEventStructure.theVetoPhotoElectrons).clear();    
  (theEventStructure.theUsers).clear();    
  (theEventStructure.thePhotons).clear();    
  (theEventStructure.theDaughters).clear();
  (theEventStructure.theDeposits).clear();  
}   

void DSEventHandler::ClearPhotoElectron()    {
  thePhotoElectronStructure.PMT         = 0;
  thePhotoElectronStructure.Time 	= 0;
}   

void DSEventHandler::ClearMuPhotoElectron()    {
  theMuPhotoElectronStructure.PMT         = 0;
  theMuPhotoElectronStructure.Time 	  = 0;
}   

void DSEventHandler::ClearVetoPhotoElectron()    {
  theVetoPhotoElectronStructure.PMT         = 0;
  theVetoPhotoElectronStructure.Time 	    = 0;
}   

void DSEventHandler::ClearPhoton()    {
  thePhotonStructure.VolumeID 	 = 0;
  thePhotonStructure.PID  	 = 0;
  thePhotonStructure.Time 	 = 0;
  thePhotonStructure.Wavelength  = 0;
  thePhotonStructure.Position[0] = 0;
  thePhotonStructure.Position[1] = 0;
  thePhotonStructure.Position[2] = 0;
}   


void DSEventHandler::ClearDaughter() {
  theDaughterStructure.Id	       = 0;
  theDaughterStructure.PDG	       = 0;
  //  theDaughterStructure.PID	       = 0;
  theDaughterStructure.TrackID	       = 0;
  theDaughterStructure.ParentTrackID	       = 0;
  theDaughterStructure.Process	       = 0;
  theDaughterStructure.Time	       = 0;
  theDaughterStructure.Energy	       = 0;
  theDaughterStructure.Position[0]     = 0; 
  theDaughterStructure.Position[1]     = 0; 
  theDaughterStructure.Position[2]     = 0; 
  theDaughterStructure.Direction[0]    = 0;
  theDaughterStructure.Direction[1]    = 0;
  theDaughterStructure.Direction[2]    = 0;
}   

void DSEventHandler::ClearDeposit()  {
  theDepositStructure.PID 	   = 0; 

  theDepositStructure.Track        = 0;
  theDepositStructure.ParentTrack        = 0;
  theDepositStructure.Volume 	   = 0; 
  theDepositStructure.Energy	   = 0; 
  theDepositStructure.KineticEnergy	   = 0; 
  theDepositStructure.TotalEnergy	   = 0; 
  theDepositStructure.Step	   = 0; 
  theDepositStructure.Time	   = 0; 
  theDepositStructure.Position[0]  = 0;   
  theDepositStructure.Position[1]  = 0;   
  theDepositStructure.Position[2]  = 0;   
}    
void DSEventHandler::ClearUser()  {
  theUserStructure.UserInt1        = 0;
  theUserStructure.UserInt2        = 0;
  theUserStructure.UserFloat1      = 0;
  theUserStructure.UserFloat2      = 0;
  theUserStructure.UserDouble      = 0;
}    

void DSEventHandler::ClearAll()  {
  ClearEvent();
  ClearPhotoElectron();
  ClearMuPhotoElectron();
  ClearVetoPhotoElectron();
  ClearPhoton();
  ClearDaughter();
  ClearUser();
  ClearDeposit();
}

/*
 * $Log: DSEventHandler.cc,v $
 * Revision 1.4  2014/11/13 16:47:04  dfranco
 * removed variables which were creating conflicts with the previous version of g4ds10
 *
 * Revision 1.3  2014/10/13 18:43:56  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.2  2014/05/08 10:59:19  pagnes
 * Scintillator Index added in binary header
 *
 * Revision 1.1  2014/05/07 12:21:02  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.7  2013/08/06 13:58:20  dfranco
 * Added 3 variables to store the deposited energy in LAr, scintillator, and water. The last two are not yet implemented. g4rooter has been updated with 3 new variables: tpcene, vetoene, and muene
 *
 * Revision 1.6  2013/07/24 09:49:01  dfranco
 * Added S1 and S2 equivalent energies, the material index for each deposit, the command killS1S2 to kill photons and electrons generated by DSLight (after storing the equivalent energies)
 *
 * Revision 1.5  2013/06/25 15:57:38  dfranco
 * Fixed a bug in clearing Users variable
 *
 * Revision 1.4  2013/04/04 09:04:16  dfranco
 * added step length info to the deposit structure
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
