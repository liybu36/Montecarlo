//---------------------------------------------------------------------------//
//    Adapted by davide.franco@mi.infn.it from the MaGe code:
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//                                                                           //
//                                                                           //
//                         MaGe Simulation                                   //
//                                                                           //
//      This code implementation is the intellectual property of the         //
//      MAJORANA and Gerda Collaborations. It is based on Geant4, an         //
//      intellectual property of the RD44 GEANT4 collaboration.              //
//                                                                           //
//                        *********************                              //
//                                                                           //
//    Neither the authors of this software system, nor their employing       //
//    institutes, nor the agencies providing financial support for this      //
//    work  make  any representation or  warranty, express or implied,       //
//    regarding this software system or assume any liability for its use.    //
//    By copying, distributing or modifying the Program (or any work based   //
//    on on the Program) you indicate your acceptance of this statement,     //
//    and all its terms.                                                     //
//                                                                           //
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//---------------------------------------------------------------------------//

#ifndef DSGeneratorNeutronsAtGSMessenger_h
#define DSGeneratorNeutronsAtGSMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DSGeneratorNeutronsAtGS;

class G4ParticleTable;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

class DSGeneratorNeutronsAtGSMessenger: public G4UImessenger {
  
public:
  DSGeneratorNeutronsAtGSMessenger(DSGeneratorNeutronsAtGS *fPtclGun);
  ~DSGeneratorNeutronsAtGSMessenger();
  
  void SetNewValue(G4UIcommand *command, G4String newValues);
 
  
private:
  DSGeneratorNeutronsAtGS *generator;
  G4UIdirectory* fDirectory;
  G4UIcmdWithADoubleAndUnit* radiusCmd;
  G4UIcmdWithADoubleAndUnit* heightCmd;
  G4UIcmdWithAString* directionCmd;
  G4UIcmdWithABool* fissionCmd;
};

#endif
/*
 * Revision 1.1  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *

 */
