#ifndef DSStackingRDM_h
#define DSStackingRDM_h 1
#include "G4ClassificationOfNewTrack.hh"
#include "G4UserStackingAction.hh"
#include "DSVStackingAction.hh"
#include "DSStackingAction.hh"
#include <vector>

using namespace std ;
class G4StackManager;
class G4Track;
class DSStackingRDMMessenger;

class DSStackingRDM : public DSVStackingAction
{
  public:
      DSStackingRDM();
      virtual ~DSStackingRDM();

  public: // with description
//---------------------------------------------------------------
// vitual methods to be implemented by user
//---------------------------------------------------------------
//
      virtual G4ClassificationOfNewTrack 
        DSClassifyNewTrack(const G4Track* aTrack);
//
//    Reply G4ClassificationOfNewTrack determined by the
//  newly coming G4Track.
//
//    enum G4ClassificationOfNewTrack
//    {
//      fUrgent,    // put into the urgent stack
//      fWaiting,   // put into the waiting stack
//      fPostpone,  // postpone to the next event
//      fKill       // kill without stacking
//    };
//
//    The parent_ID of the track indicates the origin of it.
//                
//    G4int parent_ID = aTrack->get_parentID();
//   
//      parent_ID = 0 : primary particle
//                > 0 : secondary particle
//                < 0 : postponed from the previous event
//
//---------------------------------------------------------------
//
      virtual void DSNewStage();
//
//    This method is called by G4StackManager when the urgentStack
//  becomes empty and contents in the waitingStack are transtered
//  to the urgentStack.
//    Note that this method is not called at the begining of each
//  event, but "PrepareNewEvent" is called.
//
//    In case re-classification of the stacked tracks is needed,
//  use the following method to request to G4StackManager.
//
//    stackManager->ReClassify();
//
//  All of the stacked tracks in the waitingStack will be re-classified 
//  by "ClassifyNewTrack" method.
//    To abort current event, use the following method.
//
//    stackManager->clear();
//
//  Note that this way is valid and safe only for the case it is called
//  from this user class. The more global way of event abortion is
//
//    G4UImanager * UImanager = G4UImanager::GetUIpointer();
//    UImanager->ApplyCommand("/event/abort");
//
//---------------------------------------------------------------
//
      virtual void DSPrepareNewEvent();
//
//    This method is called by G4StackManager at the begining of
//  each event.
//    Be careful that the urgentStack and the waitingStack of 
//  G4StackManager are empty at this moment, because this method
//  is called before accepting primary particles. Also, note that
//  the postponeStack of G4StackManager may have some postponed
//  tracks.
//
//---------------------------------------------------------------
      void   SetIsReclassify(G4bool val) { IsReclassify = val;  }      
      G4bool GetIsReclassify()           { return IsReclassify; }      
      
    
      void   KillParticles( G4int val ) { fPDGToBeKilled.push_back(val) ;}    
      void   KillLEParticles( G4int val, G4double ene ) { fLEPDGToBeKilled.push_back(val); fLEnergyToBeKilled.push_back(ene);}    
      
    private:
      
      DSStackingRDMMessenger  *fMessenger ;
      
      G4bool            IsReclassify ;
      G4int             fCounter, stage;
      G4bool            IsValid ;
      G4int             fDaughters;  
      G4int             fAlpha, fBeta, fGamma;
      G4bool            IsTheFirstDaughter;
      vector<G4int>     fPDGToBeKilled ;
      vector<G4int>     fLEPDGToBeKilled;
      vector<G4double>  fLEnergyToBeKilled;
      G4double          fAngle;
      G4bool            fLevelZero;
      G4int             fA;
      G4int             fZ;
      std::vector<int>  fPIDs;
      G4int             fPreTrack ;
      
     
     
      
};

#endif

/*
 * $Log: DSStackingRDM.hh,v $
 * Revision 1.1  2014/05/07 12:20:55  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.5  2013/07/07 09:42:20  dfranco
 * Added a default stacking for killing particles from long living nucleus decays
 *
 *
 */
