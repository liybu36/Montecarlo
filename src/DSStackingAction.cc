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
// $Id: DSStackingAction.cc,v 1.1 2014/05/07 12:21:05 dfranco Exp $
// GEANT4 tag $Name:  $
//
#include "DSStackingAction.hh"
#include "DSVStackingAction.hh"
#include "DSStackingActionMessenger.hh"
#include "G4StackManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


using namespace std;

DSStackingAction::DSStackingAction() {
   fMessenger = new DSStackingActionMessenger(this);
}


DSStackingAction::~DSStackingAction(){
    delete fMessenger;
}


G4ClassificationOfNewTrack  DSStackingAction::ClassifyNewTrack(const G4Track* aTrack) {
    
  return fStacking->DSClassifyNewTrack(aTrack);
}

void DSStackingAction::NewStage() {
  fStacking->DSNewStage();
  //stackManager->ReClassify();
  //return ;
}

void DSStackingAction::PrepareNewEvent() {
  fStacking->DSPrepareNewEvent();

}


/*
 * $Log: DSStackingAction.cc,v $
 * Revision 1.1  2014/05/07 12:21:05  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.5  2013/07/07 09:42:16  dfranco
 * Added a default stacking for killing particles from long living nucleus decays
 *
 *
 */
