#include "DSStorage.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSEventHandler.hh"
#include "DSStackingRDMChainMessenger.hh"
#include "DSStackingRDMChain.hh"
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
#include "G4Ions.hh" 
#include "globals.hh"
#include "CLHEP/Random/RandExponential.h"
using namespace std;

DSStackingRDMChain::DSStackingRDMChain() {
   changedParticle = new G4ParticleChange;
   IsShort = false ;
   fCounter = 0 ;
   DSStorage::Get()->SetRDMDecay(true);
   IsReclassify = false ;
   fMaxLifeTime = 1.e30*s;//10*24*3600*s ;
   isDaughter   = false ;
   fMessenger = new DSStackingRDMChainMessenger(this);
   isPostponed = false ;
   isFirst     = true ;
   //  fGateTime   = 200.*ns;  //1.*s;
   fGateTime   = 100.*ns;

   DSLog(routine) << "RDMChain Stacking Methode Active" << endlog ;   
}


DSStackingRDMChain::~DSStackingRDMChain(){;}

G4ClassificationOfNewTrack DSStackingRDMChain::DSClassifyNewTrack (const G4Track* aTrack) {
  
  if(fMaxLifeTime == 0) fMaxLifeTime = 1.e30*s ;
  
  
  // skip optical photons
  if (aTrack->GetDefinition()->GetPDGEncoding() == 50) return fUrgent; 
  
  // kill neutrinos 
  if (aTrack->GetDefinition() == G4NeutrinoE::NeutrinoE())	       return fKill; 
  if (aTrack->GetDefinition() == G4NeutrinoMu::NeutrinoMu())	       return fKill; 
  if (aTrack->GetDefinition() == G4NeutrinoTau::NeutrinoTau())         return fKill; 
  if (aTrack->GetDefinition() == G4AntiNeutrinoE::AntiNeutrinoE())     return fKill; 
  if (aTrack->GetDefinition() == G4AntiNeutrinoMu::AntiNeutrinoMu())   return fKill; 
  if (aTrack->GetDefinition() == G4AntiNeutrinoTau::AntiNeutrinoTau()) return fKill; 

/*
  if (aTrack->GetDefinition()->GetPDGEncoding() > 1e6 && aTrack->GetDefinition()->GetAtomicNumber () > 2) {
    G4Ions *ion = (G4Ions*)aTrack->GetDefinition();
    cout << aTrack->GetDefinition()->GetParticleName() << " " 
    << aTrack->GetDefinition()->GetAtomicNumber() << " "
    <<  aTrack->GetDefinition()->GetAtomicMass() << " "
    <<  aTrack->GetParentID() << " "
    <<  "LT " << aTrack->GetDefinition()->GetPDGLifeTime()/s << " "
    << ion->GetExcitationEnergy ()
    <<  " ID: " <<aTrack->GetTrackID() << " "
    <<  " PID: " <<aTrack->GetParentID() << " "
    <<  " gtime: " << aTrack->GetGlobalTime()/s << " "
    <<  " life: " << aTrack->GetDefinition()->GetPDGLifeTime()/s<< " "
    <<  ion->GetExcitationEnergy () << " " 
    << aTrack->GetDefinition()->GetPDGEncoding() << " "
    << mymap[aTrack->GetDefinition()->GetPDGEncoding()] 
    <<  endl ;
  }*/
  
  
  
  if (aTrack->GetDefinition()->GetPDGEncoding() > 1e6 && aTrack->GetDefinition()->GetAtomicNumber () > 2) {
    if(aTrack->GetDefinition()->GetPDGLifeTime()  > 0) 
       mymap[aTrack->GetDefinition()->GetPDGEncoding()] = aTrack->GetDefinition()->GetPDGLifeTime() ;

      
    int pid = aTrack->GetParentID();
    //int tid = aTrack->GetTrackID();
    G4Ions *ion = (G4Ions*)aTrack->GetDefinition();    
    
    //cout << aTrack->GetDefinition()->GetParticleName() << " "  
    //     << pid << " " << tid << endl ;
	 
    if(aTrack->GetDefinition()->GetPDGStable() && ion->GetExcitationEnergy () == 0) {
      //cout << "Uccido lo stabile: " << aTrack->GetDefinition()->GetParticleName() <<  endl ;
      isFirst = true ;
      return fKill ;
    }

    if(mymap[aTrack->GetDefinition()->GetPDGEncoding()] > fMaxLifeTime) {
      //cout << "Uccido il nucleo che vive decisamente troppo: " << aTrack->GetDefinition()->GetParticleName() <<  endl ;
      isFirst = true ;
      return fKill ;
    }



    if(pid == 0) {
      if(isFirst) {
        //cout << "lancio il padre " << aTrack->GetDefinition()->GetParticleName() <<endl ;
        aTrack->GetDefinition()->SetPDGLifeTime(0*s);
        isFirst = false ;
	DSEventHandler::Get()->SetPDG(aTrack->GetDefinition()->GetPDGEncoding());
	DSEventHandler::Get()->SetTime(CLHEP::RandExponential::shoot(mymap[aTrack->GetDefinition()->GetPDGEncoding()]));
	DSEventHandler::Get()->SetPosition(aTrack->GetPosition ()/cm );
	DSEventHandler::Get()->SetDirection(aTrack->GetMomentumDirection () );
        return fUrgent ;
      } else    {    
        //cout << "Uccido il capostipite " <<endl ;
        if(fCounter > 2) { // protezione contro strani isotopi eccitati che non decadono
          //cout << "Salvo il capostipite" << endl ;
          DSEventHandler::Get()->SetPDG(aTrack->GetDefinition()->GetPDGEncoding());
	  DSEventHandler::Get()->SetTime(CLHEP::RandExponential::shoot(mymap[aTrack->GetDefinition()->GetPDGEncoding()]));
	  DSEventHandler::Get()->SetPosition(aTrack->GetPosition ()/cm );
	  DSEventHandler::Get()->SetDirection(aTrack->GetMomentumDirection () );

          return fUrgent ;
        }
        fCounter++;
        return fKill ;      
      }
    }
    
    
    if(pid == -1 && !fDaughters) {
      fCounter = 0 ;
      //cout << "Lancio il figlio " <<aTrack->GetDefinition()->GetParticleName() << endl ;
      fDaughters = 1 ;
      DSEventHandler::Get()->SetPDG(aTrack->GetDefinition()->GetPDGEncoding());
      DSEventHandler::Get()->SetTime(CLHEP::RandExponential::shoot(mymap[aTrack->GetDefinition()->GetPDGEncoding()]));
      DSEventHandler::Get()->SetPosition(aTrack->GetPosition ()/cm );
      DSEventHandler::Get()->SetDirection(aTrack->GetMomentumDirection () );
      return fUrgent ;
    }
    
    if(mymap[aTrack->GetDefinition()->GetPDGEncoding()] > fGateTime) {
      //cout << "Salvo il figlio " <<  aTrack->GetDefinition()->GetParticleName() << endl ;
      aTrack->GetDefinition()->SetPDGLifeTime(0*s);
      isPostponed = true ;
      return fPostpone ;
    }    
    
    if(isPostponed) {
      //cout << "Uccido il figlio del posticipato  " <<  aTrack->GetDefinition()->GetParticleName() << endl ;
      return fKill ;
    }
  
  
  }  
  


/*  
  DSLog(development) << "Stacking " 
    <<   aTrack->GetDefinition()->GetParticleName() << " "
    <<  " gtime: " << aTrack->GetGlobalTime() << " "
    <<  " ptime: " << aTrack->GetProperTime() << " "
    <<  " ltime: " << aTrack->GetLocalTime()/ns  << " "
    <<  " steps: " << aTrack->GetCurrentStepNumber() << " "
    <<  " E: " <<aTrack->GetKineticEnergy()/eV<< " "
    <<  " ID: " <<aTrack->GetTrackID() << " "
    <<  " PID: " <<aTrack->GetParentID() << " "
    << endlog ;
  
  
  
  
  
  
  if( aTrack->GetDefinition()->GetAtomicNumber () <= 2) {
    for (std::vector<int>::iterator it = fPIDs.begin() ; it != fPIDs.end(); ++it) {
      if(aTrack->GetParentID() == *it) {
	if(aTrack->GetTrackID() == fPreTrack) break ;
	int index = -1 ;
	DSEventHandler::Get()->SetDId(int(DSEventHandler::Get()->GetVDaughters().size()));
	DSEventHandler::Get()->SetDPDG(aTrack->GetDefinition()->GetPDGEncoding());
	DSEventHandler::Get()->SetDProcess(index);
	DSEventHandler::Get()->SetDTime(aTrack->GetGlobalTime()/ns);
	DSEventHandler::Get()->SetDEnergy(aTrack->GetKineticEnergy ()/MeV);
	DSEventHandler::Get()->SetDPosition(aTrack->GetVertexPosition()/m);
	DSEventHandler::Get()->SetDDirection(aTrack->GetVertexMomentumDirection() );
	DSEventHandler::Get()->SetDaughters();  
	break;      
      }
    }
  }
 
  fPreTrack = aTrack->GetTrackID();
  
  // save daughters, if not ions
  if( fCounter &&  aTrack->GetDefinition()->GetAtomicNumber () < 2 ) {
    DSStorage::Get()->SetRealPDGMeanLife(aTrack->GetDefinition()->GetPDGLifeTime());
    G4double coinctime =  CLHEP::RandExponential::shoot(DSStorage::Get()->GetRealPDGMeanLife())
                         +DSStorage::Get()->GetPreAbsTime();
    DSStorage::Get()->SetPreAbsTime(coinctime);
    DSEventHandler::Get()->SetDPDG(aTrack->GetDefinition()->GetPDGEncoding());
    DSEventHandler::Get()->SetDTime(coinctime);
    DSEventHandler::Get()->SetDPosition(aTrack->GetPosition ()/cm );
    DSEventHandler::Get()->SetDDirection(aTrack->GetMomentumDirection () );

    aTrack->GetDefinition()->SetPDGLifeTime(0*s);
    fPIDs.push_back(aTrack->GetTrackID());
    return fUrgent ; 
  }
  
  
  // kill stable particle 
  if(aTrack->GetDefinition()->GetAtomicMass () > 4 && aTrack->GetDefinition()->GetPDGStable() ) {
    isDaughter = false;
    return fUrgent ; 
    
  }
  // kill primary if daughter exists
  if(isDaughter) {
    if(aTrack->GetDefinition()->GetAtomicMass () > 4 && aTrack->GetParentID() == 0) {
      return fKill ; 
    }
  }
  
  
  // process particles with excited states
  if(aTrack->GetDefinition()->GetAtomicNumber () > 2) {
    G4Ions *ion = (G4Ions*)aTrack->GetDefinition();    
    if(ion->GetExcitationEnergy () > 0) {
      if(aTrack->GetDefinition()->GetPDGLifeTime() > 1*s)  {
        isDaughter = true ;
        return fPostpone;
      }
      fPIDs.push_back(aTrack->GetTrackID());
      return fUrgent ;   
    }
  }
  

  // primary particle 
  if(!fCounter && aTrack->GetDefinition()->GetAtomicNumber () > 2 && aTrack->GetParentID() == 0) {
    isDaughter = false ;
    G4Ions *ion = (G4Ions*)aTrack->GetDefinition();
    fZ = ion->GetAtomicNumber ();
    fA = ion->GetAtomicMass ();
    DSStorage::Get()->SetRealPDGMeanLife(aTrack->GetDefinition()->GetPDGLifeTime());
    G4double coinctime =  CLHEP::RandExponential::shoot(DSStorage::Get()->GetRealPDGMeanLife())
                         +DSStorage::Get()->GetPreAbsTime();
    DSStorage::Get()->SetPreAbsTime(coinctime);
    DSEventHandler::Get()->SetPDG(aTrack->GetDefinition()->GetPDGEncoding());
    DSEventHandler::Get()->SetTime(coinctime);
    DSEventHandler::Get()->SetPosition(aTrack->GetPosition ()/cm );
    DSEventHandler::Get()->SetDirection(aTrack->GetMomentumDirection () );

    aTrack->GetDefinition()->SetPDGLifeTime(0*s);
    fPIDs.push_back(aTrack->GetTrackID());
    fCounter = 1 ;
    return fUrgent;
  }
   
   
   
   
  // postpone daughter 
  if(    aTrack->GetDefinition()->GetAtomicMass () > 4   
        && aTrack->GetParentID() != 0
      && (fZ != aTrack->GetDefinition()->GetAtomicNumber() || fA != aTrack->GetDefinition()->GetAtomicMass())      
   ) {
    G4Ions *ion = (G4Ions*)aTrack->GetDefinition();
    cout << "aaaa " << ion->GetParticleName() << endl ; 
    //if(ion->GetExcitationEnergy () > 0) return fUrgent ;
    
    fZ = ion->GetAtomicNumber ();
    fA = ion->GetAtomicMass ();
    if(aTrack->GetDefinition()->GetPDGLifeTime() < fMaxLifeTime ) {
      aTrack->GetDefinition()->SetPDGLifeTime(0*s);
      isDaughter = true ;  
      return fPostpone; 
    } else {
      isDaughter = false ;  
      return fKill;      
    }
    
  }
   
*/  



  return fUrgent;
}

void DSStackingRDMChain::DSNewStage() { 

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
  //DSStackCheckStatus();
  //DSStackReClassify();	    
  //DSStackCheckStatus();
  
  return ;
}

void DSStackingRDMChain::DSPrepareNewEvent() {  
  //fCounter         = 0;
  IsValid          = false ;
  stage            = 0;
  isSecondDaughter = false;
  fDaughters       = 0;
  isPostponed      = false ;
  //isFirst          = false ;
  fAlpha           = 0;
  fBeta            = 0;
  fGamma           = 0;
  fPIDs.clear();
}

/*
 * Revision 1.1  2013/03/22 13:50:09  dfranco
 * Added radioactive decay and radioactive decay chain generators. They require the correspondent stacking actions. Two mac files are included as examples
 * *
 */

