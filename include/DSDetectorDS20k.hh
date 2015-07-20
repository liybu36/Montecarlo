#ifndef DSDetectorDS20k_H
#define DSDetectorDS20k_H

#include "G4LogicalVolume.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4VPhysicalVolume.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"


class	DSMaterial;

class   DSDetectorDS20k  {

public:
  
  DSDetectorDS20k(G4VPhysicalVolume*);
  ~DSDetectorDS20k();
  
  
  //  G4VPhysicalVolume*          GetDetectorComponent()  { return  fPhysicOuterLiqArgon; }
  //G4VPhysicalVolume*          GetDetectorComponent()  { return  fPhysicCryoOuter; }

private:
  void                        DefineSurfaces();
  
  
  double                      GetOctagonInnerRadius(double );
  double                      GetOctagonOuterRadius(double );
 
  G4Polyhedra                 *fSolidGrid ; 
  G4LogicalVolume             *fLogicGrid; 
  
  G4VPhysicalVolume*          fMotherVolume;
  G4VPhysicalVolume *fPhysicActiveLAr ; 
  G4VPhysicalVolume *fPhysicSiPmTop ; 
  G4VPhysicalVolume *fPhysicSiPmBottom ; 
  G4VPhysicalVolume * fPhysicGrid ; 
  G4VPhysicalVolume *fPhysicGasPocket     ; 
  G4VPhysicalVolume *fPhysicTPBSide     ; 
  G4VPhysicalVolume *fPhysicTPBTop     ; 
  G4VPhysicalVolume *fPhysicTPBBottom     ; 
  G4VPhysicalVolume *fPhysicLArLayer;
  G4VPhysicalVolume *fPhysicTopWindow  ;  
  G4VPhysicalVolume *fPhysicBotWinwow ;   
  G4VPhysicalVolume *fPhysicTeflonBottom;
/*
// Cryostats
  G4Ellipsoid*                fSolidCryoOuterTop;
  G4Tubs*                     fSolidCryoOuterBody;  
  G4UnionSolid*               fSolidCryoOuter;
  G4Ellipsoid*                fSolidCryoOuterBot;
  G4UnionSolid*               fSolidCryoOuterFull;
  G4LogicalVolume*            fLogicCryoOuter;
  G4VPhysicalVolume*          fPhysicCryoOuter;

  G4Ellipsoid*                fSolidVacTop;
  G4Tubs*                     fSolidVacBody;
  G4UnionSolid*               fSolidVac;
  G4Ellipsoid*                fSolidVacBot;
  G4UnionSolid*               fSolidVacFull;
  G4LogicalVolume*            fLogicVac;
  G4VPhysicalVolume*          fPhysicVac;

  G4Ellipsoid*                fSolidCryoInnerTop;
  G4Tubs*                     fSolidCryoInnerBody;
  G4UnionSolid*               fSolidCryoInner;
  G4Ellipsoid*                fSolidCryoInnerBot;
  G4UnionSolid*               fSolidCryoInnerFull;
  G4LogicalVolume*            fLogicCryoInner;
  G4VPhysicalVolume*          fPhysicCryoInner;

  //Liquid Argon
  G4Ellipsoid*                fSolidLArTop;
  G4Tubs*                     fSolidLArBody;
  G4UnionSolid*               fSolidLAr;
  G4Ellipsoid*                fSolidLArBot;
  G4UnionSolid*               fSolidLArFull;
  G4LogicalVolume*            fLogicLAr;
  G4VPhysicalVolume*          fPhysicLAr;

  //Gas Argon at top
  G4Ellipsoid*                fSolidGAr;
  G4LogicalVolume*            fLogicGAr;
  G4VPhysicalVolume*          fPhysicGAr;

  //TPC
  G4Tubs*                     fSolidTeflon;
  G4LogicalVolume*            fLogicTeflon;
  G4VPhysicalVolume*          fPhysicTeflon;

  G4Tubs*                     fSolidRings;
  G4LogicalVolume*            fLogicRings;
  G4VPhysicalVolume*          fPhysicRings;

  G4Tubs*                     fSolidLowerWindow;
  G4LogicalVolume*            fLogicLowerWindow;
  G4VPhysicalVolume*          fPhysicLowerWindow;

  G4Tubs*                     fSolidUpperWindow;
  G4LogicalVolume*            fLogicUpperWindow;
  G4VPhysicalVolume*          fPhysicUpperWindow;

  G4Tubs*                     fSolidGasPocket;
  G4LogicalVolume*            fLogicGasPocket;
  G4VPhysicalVolume*          fPhysicGasPocket;
 
  // PMT Assembly
  G4Tubs*                     fSolidPMTAssemblyTub;
  G4IntersectionSolid*        fSolidPMTAssemblyTop;
  G4IntersectionSolid*        fSolidPMTAssemblyBottom;
  G4LogicalVolume*            fLogicPMTAssemblyTop;
  G4VPhysicalVolume*          fPhysicPMTAssemblyTop;
  G4LogicalVolume*            fLogicPMTAssemblyBottom;
  G4VPhysicalVolume*          fPhysicPMTAssemblyBottom;
  
  G4Tubs*                     fSolidTeflonCapDisk;
  G4LogicalVolume*            fLogicTeflonCapDiskTop;
  G4LogicalVolume*            fLogicTeflonCapDiskBottom;
  G4VPhysicalVolume*          fPhysicTeflonCapDiskTop;
  G4VPhysicalVolume*          fPhysicTeflonCapDiskBottom;

  //Optical Boundaries
  //G4OpticalSurface*           fOpElectropolishedStainlessSteelSurface;
  G4LogicalBorderSurface*     fDS5kOuterSurface;
  G4LogicalBorderSurface*     fDS5kOuterSurface2;
 */
  
};

#endif
