//
//  DSDetectorDeployment.h
//  
//
//  Created by HaoQian on 9/4/14.
//
//

#ifndef ____DSDetectorDeployment__
#define ____DSDetectorDeployment__

#include <iostream>
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4UnionSolid.hh"

class DSMaterial;

class DSDetectorDeployment
{
public:
	DSDetectorDeployment(G4VPhysicalVolume*);
	~DSDetectorDeployment();
	G4VPhysicalVolume* GetDetectorComponent() {
		return fPhysicpipe;
	}
    
private:
	void DefineSurfaces();
	DSMaterial* dsmaterial;
	G4VPhysicalVolume* fMotherVolume;
	G4Tubs* fSolidpipeBody;
	G4LogicalVolume* fLogicpipe;
	G4VPhysicalVolume* fPhysicpipe;
	G4Tubs* fSolidpipeScintBody;
	G4LogicalVolume* fLogicpipeScint;
        G4VPhysicalVolume* fPhysicpipeScint;


  /*    G4OpticalSurface*           fOpElectropolishedStainlessSteelSurface;
     G4OpticalSurface*           fOpUntreatedStainlessSteelSurface;
     G4OpticalSurface*           fOpAluminumFoilSurface;
     G4OpticalSurface*           fOpLumirrorSurface;
     G4LogicalBorderSurface*     fpipeInnerSurface;
     G4LogicalBorderSurface*     fpipeOuterSurface;
  */
};
#endif /* defined(____DSDetectorDeployment__) */
