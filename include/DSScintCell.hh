#ifndef DSScintCell_H
#define DSScintCell_H

#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4Para.hh"
#include "G4Trap.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"


class	DSMaterial;

class   DSScintCell  {
  
public:
  
  DSScintCell(G4VPhysicalVolume*);
  ~DSScintCell();

  G4VPhysicalVolume* GetDetectorComponent()  { return  fPhysicCell; }
  
  
private:
  
  void  DefineSurfaces();
    
  
  G4VPhysicalVolume*          fMotherVolume;
    
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

  G4Tubs*                     fSolidPMT;
  G4LogicalVolume*            fLogicPMT;
  G4VPhysicalVolume*          fPhysicPMT;

  G4OpticalSurface*           fOpLumirrorSurface;
  G4LogicalBorderSurface*     fLumirrorSurface;

  G4OpticalSurface*           fOpPMTSurface;
  G4LogicalBorderSurface*     fPMTSurface;

};

#endif
