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
// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the DS mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// ********************************************************************

#ifndef DSRunActionMessenger_h
#define DSRunActionMessenger_h 1
#include "DSIO.hh"

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include <sstream>
#include "globals.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

class DSRunAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class DSRunActionMessenger: public G4UImessenger {

  public:
    DSRunActionMessenger(DSRunAction*);
   ~DSRunActionMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    DSRunAction*	       runAction;	    
    G4UIcmdWithABool*	       fSetAutoSeedCmd;     
    G4UIcmdWithAnInteger*      fHEPRandomSeedCmd;   
    G4UIcmdWithAString*	       fSetFileNameCmd;     
    G4UIcmdWithAString*	       fSetG4DSNameCmd;    
    G4UIcmdWithAString*	       fSetGenebNameCmd;    
    G4UIcmdWithAString*	       fSetGenebDecayCmd;   
    G4UIcmdWithAString*	       fSetCommentCmd;   
    G4UIdirectory*	       fDirectory;	    
    G4UIcmdWithABool*	       fSetIsBinaryCmd;     
    G4UIcmdWithAnInteger*      fRunCmd; 	    
    G4UIcmdWithADoubleAndUnit* fRateCmd;
    //G4UIcmdWithADoubleAndUnit* fLiveTimeCmd;
    G4UIcmdWithABool*	       fWriteDepositsCmd;
    G4UIcmdWithABool*	       fWriteIsotopesCmd;
    G4UIcmdWithABool*	       fWriteEBCmd;
};

#endif
/*
 * $Log: DSRunActionMessenger.hh,v $
 * Revision 1.1  2014/05/07 12:20:54  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
