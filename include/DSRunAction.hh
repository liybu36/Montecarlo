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

#ifndef DSRunAction_h
#define DSRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"
#include "G4UIcmdWithAString.hh"

class G4Timer;
class G4Run;
class DSRunActionMessenger;

class DSRunAction : public G4UserRunAction
{
  public:

    DSRunAction();
    virtual ~DSRunAction();

  public:

    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);
    

    
  public:
    inline void SetAutoSeed (const G4bool val)    {autoSeed    = val;}

  private:

    DSRunActionMessenger* runMessenger;
    G4bool                autoSeed;
    G4Timer*              timer;
};

#endif 
/*
 * $Log: DSRunAction.hh,v $
 * Revision 1.1  2014/05/07 12:20:54  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
