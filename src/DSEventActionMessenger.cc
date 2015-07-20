//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//

#include "DSEventActionMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "CLHEP/Random/Random.h"
#include <iostream>
#include <time.h>
#include <sys/times.h>
#include "DSEventAction.hh"
#include "DSRunAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcommand.hh"
#include "globals.hh"
#include "DSLogger.hh"


DSEventActionMessenger::DSEventActionMessenger(DSEventAction* EvAct)
:eventAction(EvAct) { 
}


DSEventActionMessenger::~DSEventActionMessenger() {
}


void DSEventActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue){ 
}
/*
 * $Log: DSEventActionMessenger.cc,v $
 * Revision 1.1  2014/05/07 12:21:02  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
