#ifndef DSStackingRDMChain_h
#define DSStackingRDMChain_h 1
#include "G4ClassificationOfNewTrack.hh"
#include "G4UserStackingAction.hh"
#include "G4ParticleChange.hh"
#include "DSVStackingAction.hh"
#include "DSStackingAction.hh"
#include <vector>
#include <map>

class G4StackManager;
class G4Track;
class DSStackingRDMChainMessenger;
// class description:
//
//  This is the base class of one of the user's optional action classes.
// This class gives the hooks for G4StackManager which controls the stacks
// of G4Track objects.
//

class DSStackingRDMChain : public DSVStackingAction
{
  public:
      DSStackingRDMChain();
      virtual ~DSStackingRDMChain();

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
      
      
      void   SetMaxLifeTime(G4double ltime) {  fMaxLifeTime = ltime;  }      
      
    private:
      
      DSStackingRDMChainMessenger  *fMessenger ;
      
      G4bool            IsReclassify ;
      G4bool            IsShort ;
      G4int             fCounter; 
      G4int             fNSequence; 
      G4int             stage;
      G4bool            IsValid ;
      G4bool            isDaughter; 
      G4bool            isSecondDaughter; 
      G4int             fAlpha;
      G4int             fBeta; 
      G4int             fGamma; 
      G4int             fDaughters;
      G4double          fRadius;
      G4double          fMaxLifeTime;
      G4ParticleChange* changedParticle;
      G4int             fA;
      G4int             fZ;
      std::vector<int>  fPIDs;
      G4int             fPreTrack ;
      G4bool            isPostponed;
      G4bool            isFirst ;
      G4double          fGateTime ;
      std::map<long int,double> mymap;
};

#endif

/*
 * Revision 1.1  2013/03/22 13:50:09  dfranco
 * Added radioactive decay and radioactive decay chain generators. They require the correspondent stacking actions. Two mac files are included as examples
 *
 *
 */
