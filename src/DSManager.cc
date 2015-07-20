#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "DSLogger.hh"
#include "DSManagerMessenger.hh"
//---------------------------------------------------------------------------//

#include "Randomize.hh"
#include <time.h>
#include <iostream>

//---------------------------------------------------------------------------//

#include "DSManager.hh"
#include "DSIO.hh"

//---------------------------------------------------------------------------//
DSManager* DSManager::Manager = 0;


DSManager::DSManager() {
  Manager         = this;   
  fDSMessenger    = new DSManagerMessenger(this);
}

DSManager* DSManager::Get() {
  if (!Manager)
    Manager = new DSManager();
  return Manager;
}

//---------------------------------------------------------------------------//

//DSManager::DSManager(const DSManager & other)
//{;}

//---------------------------------------------------------------------------//

DSManager::~DSManager()
{
  delete fDSMessenger;
}



/*
 * $Log: DSManager.cc,v $
 * Revision 1.1  2014/05/07 12:21:04  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
