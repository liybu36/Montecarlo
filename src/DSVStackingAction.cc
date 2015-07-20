#include "DSVStackingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4StackManager.hh"
#include "DSEventHandler.hh"
#include "G4ProcessManager.hh"
//---------------------------------------------------------------------------//

DSVStackingAction::DSVStackingAction() {
   UImanager = G4UImanager::GetUIpointer();
  //fManager = new G4StackManager;
  //fManager = DSStackingAction::GetStackManager();
}

//---------------------------------------------------------------------------//

//DSVStackingAction::DSVStackingAction(const DSVStackingAction & other)
//{;}

//---------------------------------------------------------------------------//

DSVStackingAction::~DSVStackingAction(){;}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

/*
 * $Log: DSVStackingAction.cc,v $
 * Revision 1.1  2014/05/07 12:21:06  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
