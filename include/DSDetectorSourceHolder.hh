#ifndef DSDetectorSourceHolder_H
#define DSDetectorSourceHolder_H

#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
//#include "G4LogicalBorderSurface.hh"
//#include "G4OpticalSurface.hh"
#include "G4UnionSolid.hh"

class DSMaterial;

class DSDetectorSourceHolder 
{	
public:
  DSDetectorSourceHolder(G4VPhysicalVolume*);
  ~DSDetectorSourceHolder();
  G4VPhysicalVolume* GetDetectorComponent() {
    return fPhysicTubeBody;
  }

private:
  //  void DefineSurfaces();
  DSMaterial* dsmaterial;
  G4VPhysicalVolume* fMotherVolume;
  G4Tubs* fSolidTubeBody;
  G4LogicalVolume* fLogicTubeBody;
  G4VPhysicalVolume* fPhysicTubeBody;
  G4Tubs* fSolidTubeEndCap;
  G4LogicalVolume* fLogicTubeEndCap;
  G4VPhysicalVolume* fPhysicTubeEndCap;
  G4VPhysicalVolume* fPhysicTubeFrontCap;
  G4Tubs* fSolidTubeScint;
  G4LogicalVolume* fLogicTubeScint;
  G4VPhysicalVolume* fPhysicTubeScint;
  G4Tubs* fSolidTubeLead;
  G4LogicalVolume* fLogicTubeLead;
  G4VPhysicalVolume* fPhysicTubeLead;

  G4Sphere* fSolidSphere;
  G4LogicalVolume* fLogicSphere;
  G4VPhysicalVolume* fPhysicSphere;

  /*G4OpticalSurface*           fOpElectropolishedStainlessSteelSurface;
    G4OpticalSurface*           fOpUntreatedStainlessSteelSurface;
    G4OpticalSurface*           fOpAluminumFoilSurface;
    G4OpticalSurface*           fOpLumirrorSurface;
    G4OpticalSurface*           fOpPMTBackSurface;
    G4LogicalBorderSurface*     fSteelInnerSurface;
    G4LogicalBorderSurface*     fSteelInnerSurfaceFlip;
    G4LogicalBorderSurface*     fSteelOuterSurface;
    G4LogicalBorderSurface*     fCryoSheathSurface;
    G4LogicalBorderSurface*      fTopFlangeSurface;
    */
};

#endif
