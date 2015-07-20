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
 
#ifndef DSEventAction_h
#define DSEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;
class G4Timer;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DSEventActionMessenger;



class DSEventAction : public G4UserEventAction
{
  public:
    DSEventAction();
    virtual ~DSEventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    

  private:
    DSEventActionMessenger*     fMessenger;
    G4Timer*                    fTimer;
    G4int                       fTotNPE;        
    G4int                       fTotNPH;        
};


#endif

    
/*
 * $Log: DSEventAction.hh,v $
 * Revision 1.2  2015/04/28 10:17:38  pagnes
 * bug fixed
 *
 * Revision 1.1  2014/05/07 12:20:51  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
