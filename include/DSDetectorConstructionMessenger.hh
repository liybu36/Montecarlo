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

#ifndef DSDetectorConstructionMessenger_h
#define DSDetectorConstructionMessenger_h 1


#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
using namespace std;

class G4UIcmdWith3VectorAndUnit;
class DSDetectorConstruction;
class G4UIcmdWithADoubleAndUnit;
class DSDetectorConstructionMessenger: public G4UImessenger {

  public:
    DSDetectorConstructionMessenger(DSDetectorConstruction*);
   ~DSDetectorConstructionMessenger();

    void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    
    DSDetectorConstruction*        fDetectorConstruction;
    G4UIdirectory*                 fDirectory;
    G4UIcmdWithAnInteger*          fConfigurationCmd ;
    G4UIcmdWithAString*            fScintillatorCmd ;
    G4UIcmdWith3VectorAndUnit*     fSourceCmd ;
    G4UIcmdWithABool*              fExtLArScintillatingCmd;
    G4UIcmdWithADouble*            fVetoYieldFactorCmd;
    G4UIcmdWithADoubleAndUnit*     fHolderRadiusCmd;
    G4UIcmdWithADoubleAndUnit*     fHolderZCmd;
    G4UIcmdWithADoubleAndUnit*     fHolderPhiCmd;
};
#endif
/*
 * $Log: DSDetectorConstructionMessenger.hh,v $
 * Revision 1.6  2015/01/14 16:58:42  dfranco
 * added source vial and correspondent commands by Laura and Bianca; manual updated
 *
 * Revision 1.5  2014/11/21 10:19:06  dfranco
 * added a command to scale the veto scintillation yield factor and fixed the visible energy variable in the veto
 *
 * Revision 1.4  2014/11/20 15:32:12  dfranco
 * added a command to remove scintillation process from liquid argon between TPC and cryostat
 *
 * Revision 1.3  2014/11/06 17:39:51  dfranco
 * added source and fixed a bug: removed gaseous nitrogen from the neutron veto
 *
 * Revision 1.2  2014/05/07 14:27:26  dfranco
 * fixed some bugs and added GdScintillator
 *
 * Revision 1.1  2014/05/07 12:20:51  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:39  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
