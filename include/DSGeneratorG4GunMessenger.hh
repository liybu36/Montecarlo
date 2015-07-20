
#ifndef DSGeneratorG4GunMessenger_h
#define DSGeneratorG4GunMessenger_h 1
#include "DSIO.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4ParticleGun.hh"

class DSPrimaryGeneratorAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class DSGeneratorG4Gun;

using namespace std;

class DSGeneratorG4GunMessenger: public G4UImessenger {

  public:
    DSGeneratorG4GunMessenger(DSGeneratorG4Gun*  );
   ~DSGeneratorG4GunMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
  
    void IonLevelCommand(G4String );
    void IonCommand(G4String );


    DSGeneratorG4Gun*                    generator;   
    G4UIdirectory*                       fDirectory;
    G4UIcmdWith3VectorAndUnit*           fPositionCmd; 
    G4UIcmdWith3Vector*                  fDirectionCmd;
    G4UIcmdWithADoubleAndUnit*           fEnergyCmd;
    G4UIcmdWithAString*                  fParticleCmd;
    G4UIcommand*                         fIonCmd;
    G4UIcommand*                         fIonLvlCmd;
    
        
    G4bool   fShootIon; 
    G4int    fAtomicNumber;
    G4int    fAtomicMass;
    G4int    fIonCharge;
    G4double fIonExciteEnergy;
    G4int    fIonEnergyLevel;
    
};

#endif
/*
 * $Log: DSGeneratorG4GunMessenger.hh,v $
 * Revision 1.2  2014/05/21 10:28:19  dfranco
 * added the possibility to shoot ions
 *
 * Revision 1.1  2014/05/07 12:20:52  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
