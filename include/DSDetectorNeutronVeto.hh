#ifndef DSDetectorNeutronVeto_H
#define DSDetectorNeutronVeto_H

#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"

class	DSMaterial;

class   DSDetectorNeutronVeto  {

  public:

    DSDetectorNeutronVeto(G4VPhysicalVolume*);
   ~DSDetectorNeutronVeto();

    G4VPhysicalVolume* GetDetectorComponent()  { return  fPhysicBScintillator;}
    
  private:
    void                        DefineSurfaces();
    
    DSMaterial*                 dsmaterial;

    G4VPhysicalVolume*          fMotherVolume;
 
    G4Box*			fSolidNeutronVeto;
    G4LogicalVolume*		fLogicNeutronVeto;
    G4VPhysicalVolume*		fPhysicNeutronVeto;
   
    G4Orb*                      fSolidSteelVessel;
    G4LogicalVolume*		fLogicSteelVessel;
    G4VPhysicalVolume*		fPhysicSteelVessel;
    G4OpticalSurface*           fOpElectropolishedStainlessSteelSurface;
    G4OpticalSurface*           fOpUntreatedStainlessSteelSurface;
    G4OpticalSurface*           fOpAluminumFoilSurface;
    G4OpticalSurface*           fOpLumirrorSurface;
    G4OpticalSurface*           fOpPMTBackSurface;
    G4LogicalBorderSurface*     fSteelInnerSurface;
    G4LogicalBorderSurface*     fSteelInnerSurfaceFlip;
    G4LogicalBorderSurface*     fSteelOuterSurface;
    G4LogicalBorderSurface*     fCryoSheathSurface;
    G4LogicalBorderSurface*      fTopFlangeSurface;

    G4Orb*                      fSolidBScintillator;
    G4LogicalVolume*		fLogicBScintillator;
    G4VPhysicalVolume*		fPhysicBScintillator;

    G4Tubs*                     fSolidCryoSheath;
    G4LogicalVolume*            fLogicCryoSheath;
    G4VPhysicalVolume*          fPhysicCryoSheath;

    G4Sphere*                   fSolidTopFlange;
    G4LogicalVolume*            fLogicTopFlange;
    G4VPhysicalVolume*          fPhysicTopFlange;

    //The Scintillator Cell
    G4Tubs*                     fSolidCellBody;
    G4Tubs*                     fSolidCellStem;
    G4UnionSolid*               fSolidCell;
    G4LogicalVolume*	      fLogicCell;
    G4VPhysicalVolume*	      fPhysicCell;

    G4Tubs*                     fSolidScintBody;
    G4Tubs*                     fSolidScintStem;
    G4UnionSolid*               fSolidScint;
    G4LogicalVolume*            fLogicScint;
    G4VPhysicalVolume*          fPhysicScint;

    G4Tubs*                     fSolidBubble;
    G4LogicalVolume*            fLogicBubble;
  G4VPhysicalVolume*          fPhysicBubble;
  
  //Organ pipe
  G4Tubs*                      fSolidpipeBody;
  G4SubtractionSolid*          fpipeBody;
  G4LogicalVolume*             fLogicpipe;
  G4VPhysicalVolume*           fPhysicpipe;

  G4Tubs*                      fSolidpipeScintBody;
  G4SubtractionSolid*          fpipeScintBody;
  G4LogicalVolume*             fLogicpipeScint;
  G4VPhysicalVolume*           fPhysicpipeScint;

  G4Tubs*                      fSolidSurfaceHole;
  G4IntersectionSolid*          fSurfaceHole;
  //  G4SubtractionSolid*          fSurfaceHole;
  G4LogicalVolume*             fLogicSurfaceHole;
  G4VPhysicalVolume*           fPhysicSurfaceHole;

  G4LogicalBorderSurface*     fpipeInnerSurface;
  G4LogicalBorderSurface*     fpipeInnerSurfaceFlip;
  G4LogicalBorderSurface*     fpipeOuterSurface;


};

#endif
/*
 * $Log: DSDetectorNeutronVeto.hh,v $
 * Revision 1.1  2014/05/07 12:20:51  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.7  2014/04/18 16:21:40  swesterd
 * fixed an overlap in the G2 detector
 *
 * Revision 1.6  2013/08/27 04:07:00  swesterd
 * some fine tuning of the boron scintillator kB and scint yield, and some modifications to the DSG2 geometry
 *
 * Revision 1.5  2013/08/05 03:13:51  swesterd
 * some fine tuning of bscint and veto parameters
 *
 * Revision 1.4  2013/05/27 23:59:00  swesterd
 * added a (currently commented out) Lumirror sheath to the cryostat and introduced DSOpBoundaryProcess to try to figure out why the boundaries are being screwy, with some edits so that it can handle constant and vector properties with freaking out
 *
 * Revision 1.3  2013/05/07 23:06:26  swesterd
 * added optical boundaries and Lumirror in the veto
 *
 * Revision 1.2  2013/03/22 14:09:39  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
