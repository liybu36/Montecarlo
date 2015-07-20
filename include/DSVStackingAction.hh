//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DSVStackingAction.hh,v 1.1 2014/05/07 12:20:56 dfranco Exp $
// GEANT4 tag $Name:  $
//

#ifndef DSVStackingAction_h
#define DSVStackingAction_h 1
#include "globals.hh"
//#include "G4ClassificationOfNewTrack.hh"
#include "G4UserStackingAction.hh"
#include "DSStackingAction.hh"
#include "G4UImanager.hh"
#include "G4StackingMessenger.hh"

class G4StackManager ;

class DSVStackingAction  {
  public:

      DSVStackingAction();
  
      virtual ~DSVStackingAction();

  public: // virtual methods

      virtual G4ClassificationOfNewTrack  DSClassifyNewTrack(const G4Track* aTrack) = 0;

      virtual void DSNewStage() = 0;

      virtual void DSPrepareNewEvent() = 0;
      
 
  public: // concrete methods

      inline void DSStackAbort()                 { UImanager->ApplyCommand("/event/abort");         }
      inline void DSStackClearAll()              { UImanager->ApplyCommand("/event/stack/clear 2"); }
      inline void DSStackClearUrgentAndWaiting() { UImanager->ApplyCommand("/event/stack/clear 1"); }
      inline void DSStackClearWaiting()          { UImanager->ApplyCommand("/event/stack/clear 0"); }
      inline void DSStackClearUrgent()           { UImanager->ApplyCommand("/event/stack/clear -1");}
      inline void DSStackClearPostponed()        { UImanager->ApplyCommand("/event/stack/clear -2");}
      inline void DSStackCheckStatus()           { UImanager->ApplyCommand("/event/stack/status");}
      inline void DSStackReClassify()            { UImanager->ApplyCommand("/event/stack/clear -3");}
    

   private:
   
      G4StackManager *fManager;  
      G4UImanager * UImanager ;
     

};

#endif



/*
 * Revision 1.3  2013/03/22 13:50:09  dfranco
 * Added radioactive decay and radioactive decay chain generators. They require the correspondent stacking actions. Two mac files are included as examples
 *
 */
