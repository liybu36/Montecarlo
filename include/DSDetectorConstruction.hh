#ifndef DSDETECTORCONSTRUCTION_H
#define DSDETECTORCONSTRUCTION_H

#include <vector>

#include "G4VUserDetectorConstruction.hh" 
#include "G4ThreeVector.hh"

class   G4VSolid;
class	G4VPhysicalVolume;
class	G4PVPlacement;

class	G4LogicalVolume;
class	G4LogicalBorderSurface;
class	G4OpticalSurface;
class   G4LogicalSkinSurface;

class	G4Box;
class	G4Tubs;
class	G4Sphere;
class	G4Cons;
class	G4Polycone;
class	G4BREPSolidPCone;
class	G4CSGSolid;
class	G4Material;
class   DSDetectorConstructionMessenger;
class	G4SubtractionSolid;

class	DSMaterial;

class   DSDetectorConstruction : public G4VUserDetectorConstruction {

  public:

    DSDetectorConstruction();
   ~DSDetectorConstruction();

    
    void UpdateGeometry();
    void PrintDetectorParameters();
    
    void  SetDetectorConfiguration(G4int val)  { fDetectorConfiguration = val;  }
    G4int SetDetectorConfiguration()           { return fDetectorConfiguration; }
    
    void SetIsSource(bool val)                 { fIsSource = val ;}

  private:
    G4VPhysicalVolume* Construct();
    G4VPhysicalVolume* ConstructDetector();
    void DefineSizes();

    DSDetectorConstructionMessenger *fMessenger;
    DSMaterial *dsmaterial;
    
    G4int                       fDetectorConfiguration;

//_____World_____________________________________________

    G4double           		fWorldSizeX;
    G4double           		fWorldSizeY;
    G4double			fWorldSizeZ;
    G4Material*			fWorldMaterial;
    G4Box*			fSolidWorld;
    G4LogicalVolume*		fLogicWorld;
    G4VPhysicalVolume*		fPhysicWorld;
    
    bool                        fIsSource ;



};

#endif
/*
 * $Log: DSDetectorConstruction.hh,v $
 * Revision 1.2  2014/11/06 17:39:51  dfranco
 * added source and fixed a bug: removed gaseous nitrogen from the neutron veto
 *
 * Revision 1.1  2014/05/07 12:20:50  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:39  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
