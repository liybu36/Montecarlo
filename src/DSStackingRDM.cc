#include "DSStorage.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSEventHandler.hh"
#include "DSStackingRDMMessenger.hh"
#include "DSStackingRDM.hh"
#include "DSLogger.hh"
#include "G4StackManager.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4ParticleTable.hh"
#include "G4NeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4OpticalPhoton.hh" 
#include "G4GenericIon.hh" 
#include "globals.hh"
#include "CLHEP/Random/RandExponential.h"
#include "DSStorage.hh"

using namespace std;

DSStackingRDM::DSStackingRDM() {

   IsReclassify = false ;
   fMessenger = new DSStackingRDMMessenger(this);
   fLevelZero = false ;
   fPreTrack  = 0 ;
   DSLog(routine) << "RDM Stacking Methode Active" << endlog ;  
 
}


DSStackingRDM::~DSStackingRDM(){;}

G4ClassificationOfNewTrack DSStackingRDM::DSClassifyNewTrack (const G4Track* aTrack) {
  
    
  
  // skip optical photons
  if (aTrack->GetDefinition()->GetPDGEncoding() == 50) return fUrgent; 
  /*cout << "stack " 
       << aTrack->GetDefinition()->GetParticleName() << " "
       << aTrack->GetGlobalTime()/ns << " "
       << aTrack->GetTrackID() << " "
       << aTrack->GetParentID() << " "
       << endl ;       
  */
  // kill neutrinos
  if (aTrack->GetDefinition() == G4NeutrinoE::NeutrinoE())	       return fKill; 
  if (aTrack->GetDefinition() == G4NeutrinoMu::NeutrinoMu())	       return fKill; 
  if (aTrack->GetDefinition() == G4NeutrinoTau::NeutrinoTau())         return fKill; 
  if (aTrack->GetDefinition() == G4AntiNeutrinoE::AntiNeutrinoE())     return fKill; 
  if (aTrack->GetDefinition() == G4AntiNeutrinoMu::AntiNeutrinoMu())   return fKill; 
  if (aTrack->GetDefinition() == G4AntiNeutrinoTau::AntiNeutrinoTau()) return fKill; 
  
  
  
  
/*  cout << "stack " 
       << aTrack->GetDefinition()->GetParticleName() << " "
       << aTrack->GetKineticEnergy ()/MeV << " "
       << aTrack->GetTrackID() << " "
       << aTrack->GetParentID() << " "
       << endl ;       
*/    
       
  if( aTrack->GetDefinition()->GetAtomicNumber () <= 2) {
    for (std::vector<int>::iterator it = fPIDs.begin() ; it != fPIDs.end(); ++it) {
      if(aTrack->GetParentID() == *it) {
	if(aTrack->GetTrackID() == fPreTrack) break ;
        //cout << "Loop " << *it << " " <<aTrack->GetDefinition()->GetParticleName() << endl ;
	int index = -1 ;
	DSEventHandler::Get()->SetDId(int(DSEventHandler::Get()->GetVDaughters().size()));
	//	DSEventHandler::Get()->SetDPID(aTrack->GetParentID());
	DSEventHandler::Get()->SetDParentTrackID(aTrack->GetParentID());
	DSEventHandler::Get()->SetDPDG(aTrack->GetDefinition()->GetPDGEncoding());
	DSEventHandler::Get()->SetDProcess(index);
	DSEventHandler::Get()->SetDTime(aTrack->GetGlobalTime()/ns);
	DSEventHandler::Get()->SetDEnergy(aTrack->GetKineticEnergy ()/MeV);
	DSEventHandler::Get()->SetDPosition(aTrack->GetVertexPosition()/m);
	DSEventHandler::Get()->SetDDirection(aTrack->GetVertexMomentumDirection() );
	DSEventHandler::Get()->SetDaughters();  
	break ;      
      }
    }
  }
  
  fPreTrack = aTrack->GetTrackID();
  
  // save parent identification numbers (fZ fA)
  if(aTrack->GetDefinition()->GetAtomicNumber () > 2 && aTrack->GetParentID() == 0) {
    G4Ions *ion = (G4Ions*)aTrack->GetDefinition();
    fZ = ion->GetAtomicNumber ();
    fA = ion->GetAtomicMass ();
    fPIDs.push_back(aTrack->GetTrackID());
    DSStorage::Get()->SetRealPDGMeanLife(aTrack->GetDefinition()->GetPDGLifeTime());
    G4double coinctime =  CLHEP::RandExponential::shoot(DSStorage::Get()->GetRealPDGMeanLife())
                         +DSStorage::Get()->GetPreAbsTime();
    DSStorage::Get()->SetPreAbsTime(coinctime);
    aTrack->GetDefinition()->SetPDGLifeTime(0.*ns);
    //DSEventHandler::Get()->SetPDG(aTrack->GetDefinition()->GetPDGEncoding());
    DSEventHandler::Get()->SetTime(coinctime);
    //DSEventHandler::Get()->SetPosition(aTrack->GetPosition ()/cm );
    // DSEventHandler::Get()->SetDirection(aTrack->GetMomentumDirection () );
    return fUrgent ; 
  }
  
  if(  aTrack->GetDefinition()->GetAtomicNumber () > 2) {
    G4Ions *ion = (G4Ions*)aTrack->GetDefinition();
    if( ion->GetAtomicNumber () == fZ && ion->GetAtomicMass () == fA)    fPIDs.push_back(aTrack->GetTrackID());
  }


  // kill the stable state of the daughter or save the parent ID
  if(    aTrack->GetDefinition()->GetAtomicNumber () > 2 
      && aTrack->GetParentID() != 0
      && (fZ != aTrack->GetDefinition()->GetAtomicNumber() || fA != aTrack->GetDefinition()->GetAtomicMass())      
   ) {
    G4Ions *ion = (G4Ions*)aTrack->GetDefinition();
    if(ion->GetExcitationEnergy () == 0) return fKill;
    else { 
       fPIDs.push_back(aTrack->GetTrackID());
      return fUrgent ;     
    }
  }
  
  
  //if(aTrack->GetDefinition()->GetAtomicNumber () > 2 && aTrack->GetParentID() != 0) 
  //  return fKill ; 

  DSLog(development) 
   <<   aTrack->GetDefinition()->GetParticleName() << " "
   <<  " gtime: " << aTrack->GetGlobalTime()/ns << " "
   <<  " ltime: " << aTrack->GetLocalTime()/ns  << " "
   <<  " steps: " << aTrack->GetCurrentStepNumber() << " "
   <<  " E: "     <<aTrack->GetKineticEnergy()/eV<< " "
   <<  " ID: "    <<aTrack->GetTrackID() << " "
   <<  " PID: "   <<aTrack->GetParentID() << " " 
  << endlog ;


  if(aTrack->GetParentID() == 0) {
    DSStorage::Get()->SetRealPDGMeanLife(aTrack->GetDefinition()->GetPDGLifeTime());
    DSEventHandler::Get()->SetTime(CLHEP::RandExponential::shoot(DSStorage::Get()->GetRealPDGMeanLife()));
    aTrack->GetDefinition()->SetPDGLifeTime(0*s);
    //cout << aTrack->GetDefinition()->GetPDGLifeTime() << endl; 
    return fUrgent;
  }

  
  
  if(    (aTrack->GetCurrentStepNumber() == 0) 
	 &&  (aTrack->GetDefinition()->GetPDGEncoding() != 50) 
     && (abs(aTrack->GetParentID()) == 1) ) {
    
    for(G4int i = 0; i < G4int( fPDGToBeKilled.size() ); i++ ) 
      if (aTrack->GetDefinition()->GetPDGEncoding() == fPDGToBeKilled[i])   return fKill;
    
    for(G4int i = 0; i < G4int( fLEPDGToBeKilled.size() ); i++ ) 
      if (  aTrack->GetDefinition()->GetPDGEncoding() == fLEPDGToBeKilled[i] 
         && aTrack->GetKineticEnergy() < fLEnergyToBeKilled[i])   return fKill;
    
    DSLog(development)  <<   aTrack->GetDefinition()->GetParticleName() << " "
	  <<  " gtime: " << aTrack->GetGlobalTime() << " "
	  <<  " ltime: " << aTrack->GetLocalTime()/ns  << " "
	  <<  " steps: " << aTrack->GetCurrentStepNumber() << " "
	  <<  " E: " <<aTrack->GetKineticEnergy()/eV<< " "
	  <<  " ID: " <<aTrack->GetTrackID() << " "
	  <<  " PID: " <<aTrack->GetParentID() << " "
	  <<  " Proper time: " <<aTrack->GetProperTime ()
	  << endlog ;
     
  }



  /*
  G4ClassificationOfNewTrack status = fUrgent;

  if(IsReclassify) {
    switch(stage) {
      case 0: // Stage 0 
	// Here, I identified and suspended all photons at their generation point. 
	// If at least a compton does not reach the bulk, all photons are thrown away

	if(aTrack->GetParentID()==0) return fUrgent ;
	 if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhoton()) {
	  if(aTrack->GetGlobalTime()/ns > 200. ) return fKill ;
	  fCounter++;
	  return  fWaiting ;
	}



	IsValid = true ;


	break;

      case 1: // Stage 1

	// Here, you can set some further cut on the photons belonging to the waiting stack
	break; 
    }
  }  */
  
  return fUrgent;
}

void DSStackingRDM::DSNewStage() { 


 
  if(IsReclassify) {
    stage++;
    if(!IsValid)  DSStackClearAll();
    else          DSStackReClassify(); 
  }
  /* candidates
  DSStackAbort() 	       
  DSStackClearAll()	       
  DSStackClearUrgentAndWaiting()
  DSStackClearWaiting()         
  DSStackClearUrgent()	       
  DSStackClearPostponed()       
  DSStackCheckStatus()	       
  DSStackReClassify()	    
  */   
  
  
}

void DSStackingRDM::DSPrepareNewEvent() {  
  fCounter   = 0;
  IsValid    = false ;
  stage      = -1;
  fDaughters = 0;
  fAlpha     = 0;
  fBeta      = 0;
  fGamma     = 0;
  fPIDs.clear();
  //cout << "------ NEW ------- " << endl ;
  //IsTheFirstDaughter = false ;
}

/*
 * Revision 1.1  2013/03/22 13:50:09  dfranco
 * Added radioactive decay and radioactive decay chain generators. They require the correspondent stacking actions. Two mac files are included as examples
 *
 *
 */

