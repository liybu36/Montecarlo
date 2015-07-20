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
// $Id: DSStackingAction.hh,v 1.1 2014/05/07 12:20:54 dfranco Exp $
// GEANT4 tag $Name:  $
//

#ifndef DSStackingAction_h
#define DSStackingAction_h 1

#include "G4ClassificationOfNewTrack.hh"
#include "G4UserStackingAction.hh"
#include "G4StackManager.hh"

class DSStackingActionMessenger;
class G4StackManager;
class G4Track;
class DSVStackingAction;
// class description:
//
//  This is the base class of one of the user's optional action classes.
// This class gives the hooks for G4StackManager which controls the stacks
// of G4Track objects.
//

class DSStackingAction : public G4UserStackingAction  {

  public:
      DSStackingAction();
      virtual ~DSStackingAction();

      void SetDSStackingAction(DSVStackingAction* stack) { fStacking = stack;}
//
      virtual G4ClassificationOfNewTrack  ClassifyNewTrack(const G4Track* aTrack) ;
//
      virtual void NewStage() ;
//
      virtual void PrepareNewEvent() ;

               
  private: 
  
      DSStackingActionMessenger *fMessenger;
      DSVStackingAction         *fStacking;

	 

};

#endif


/*
 * $Log: DSStackingAction.hh,v $
 * Revision 1.1  2014/05/07 12:20:54  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.5  2013/07/07 09:42:20  dfranco
 * Added a default stacking for killing particles from long living nucleus decays
 *
 *
 */

