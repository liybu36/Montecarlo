#ifndef DSDetectorNeutronVetoG3_H
#define DSDetectorNeutronVetoG3_H

#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4SubtractionSolid.hh"
#include "G4OpticalSurface.hh"
#include "G4UnionSolid.hh"

class	DSMaterial;

class   DSDetectorNeutronVetoG3  {

  public:

    DSDetectorNeutronVetoG3(G4VPhysicalVolume*);
   ~DSDetectorNeutronVetoG3();

    G4VPhysicalVolume* GetDetectorComponent()  { return  fPhysicBScintillator;}
    

  private:
    void                        DefineSurfaces();
    
    DSMaterial*                 dsmaterial;

    G4VPhysicalVolume*          fMotherVolume;

    G4Tubs*                     fSubtrCylinderTop; 
    G4Tubs*                     fSubtrCylinderBottom; 

    G4Orb*                      fSolidSteelVessel_0;
    G4SubtractionSolid*         fSolidSteelVessel_1;
    G4SubtractionSolid*         fSolidSteelVessel;
    G4LogicalVolume*            fLogicSteelVessel;
    G4VPhysicalVolume*          fPhysicSteelVessel;
    G4OpticalSurface*           fOpElectropolishedStainlessSteelSurface;
    G4OpticalSurface*           fOpUntreatedStainlessSteelSurface;
    G4OpticalSurface*           fOpAluminumFoilSurface;
    G4OpticalSurface*           fOpLumirrorSurface;
    G4OpticalSurface*           fOpPMTBackSurface;
    G4LogicalBorderSurface*     fSteelInnerSurface;
    G4LogicalBorderSurface*     fSteelInnerSurfaceFlip;
    G4LogicalBorderSurface*     fSteelOuterSurface;
    G4LogicalBorderSurface*     fCryoSheathSurface;
    G4LogicalBorderSurface*     fTopFlangeSurface;

    G4Tubs*                     fSolidTopFlange;
    G4LogicalVolume*            fLogicTopFlange;
    G4VPhysicalVolume*          fPhysicTopFlange;

    G4Tubs*                     fSolidBottomFlange;
    G4LogicalVolume*            fLogicBottomFlange;
    G4VPhysicalVolume*          fPhysicBottomFlange;

    G4Orb*                      fSolidBScintillator_0;
    G4SubtractionSolid*         fSolidBScintillator_1;
    G4SubtractionSolid*         fSolidBScintillator;
    G4LogicalVolume*            fLogicBScintillator;
    G4VPhysicalVolume*          fPhysicBScintillator;

    G4Tubs*                     fSolidCryoSheath;
    G4LogicalVolume*            fLogicCryoSheath;
    G4VPhysicalVolume*          fPhysicCryoSheath;

    //The Scintillator Cell
    G4Tubs*                     fSolidCellBody;
    G4Tubs*                     fSolidCellStem;
    G4UnionSolid*               fSolidCell;
    G4LogicalVolume*            fLogicCell;
    G4VPhysicalVolume*          fPhysicCell;

    G4Tubs*                     fSolidScintBody;
    G4Tubs*                     fSolidScintStem;
    G4UnionSolid*               fSolidScint;
    G4LogicalVolume*            fLogicScint;
    G4VPhysicalVolume*          fPhysicScint;

    G4Tubs*                     fSolidBubble;
    G4LogicalVolume*            fLogicBubble;
    G4VPhysicalVolume*          fPhysicBubble;

};

#endif
