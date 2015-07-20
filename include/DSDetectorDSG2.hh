#ifndef DSDetectorDSG2_H
#define DSDetectorDSG2_H

#include "G4LogicalVolume.hh"
#include "G4Polycone.hh"
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

class   DSDetectorDSG2  {

public:
  
  DSDetectorDSG2(G4VPhysicalVolume*);
  ~DSDetectorDSG2();
  
  
  //  G4VPhysicalVolume*          GetDetectorComponent()  { return  fPhysicOuterLiqArgon; }
  G4VPhysicalVolume*          GetDetectorComponent()  { return  fPhysicCryoOuter; }

private:
  void                        DefineSurfaces();

  G4VPhysicalVolume*          fMotherVolume;
  
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
  G4LogicalBorderSurface*     fDSG2OuterSurface;
  G4LogicalBorderSurface*     fDSG2OuterSurface2;

  
  /*
  // Cathode Window
  G4Tubs*                     fSolidCathodeWindow;
  G4LogicalVolume*            fLogicCathodeWindow;
  G4VPhysicalVolume*          fPhysicCathodeWindow;
  
  G4Tubs*                     fSolidITO;
  G4LogicalVolume*            fLogicITO;
  G4VPhysicalVolume*          fPhysicCathodeITOTop;
  G4VPhysicalVolume*          fPhysicCathodeITOBottom;
  
  // Teflon Reflector
  G4Polycone*                 fSolidReflector;
  G4LogicalVolume*            fLogicReflector;
  G4VPhysicalVolume*          fPhysicReflector;
  
  // TPB Layer 
  G4Polycone*                 fSolidTPB;
  G4LogicalVolume*            fLogicTPB;
  G4VPhysicalVolume*          fPhysicTPB;
  
  // Active LAr 
  G4Polycone*                 fSolidActiveLAr;
  G4LogicalVolume*            fLogicActiveLAr;
  G4VPhysicalVolume*          fPhysicActiveLAr;
  
  // GasPocket
  G4Tubs*                     fSolidGasPocket;
  G4LogicalVolume*            fLogicGasPocket;
  G4VPhysicalVolume*          fPhysicGasPocket;
  
  // Grid
  G4Tubs*                     fSolidGrid;
  G4LogicalVolume*            fLogicGrid;
  G4VPhysicalVolume*          fPhysicGrid;
  
  // Diving Bell
  G4Tubs*                     fSolidBellTop;
  G4Tubs*                     fSolidBellRim;
  G4LogicalVolume*            fLogicBellTop;
  G4LogicalVolume*            fLogicBellRim;
  G4VPhysicalVolume*          fPhysicBellTop;
  G4VPhysicalVolume*          fPhysicBellRim;
  
  G4VPhysicalVolume*          fPhysicBellITOTop;
  G4VPhysicalVolume*          fPhysicBellITOBottom;
    
  // Field Rings
  G4Tubs*                     fSolidFieldRings;
  G4LogicalVolume*            fLogicFieldRings;
  G4VPhysicalVolume*          fPhysicFieldRings;
  */  
  
  
  
/*    // Optical surfaces
      G4OpticalSurface*           fOpWindowLArSurface;
    G4LogicalBorderSurface*     fTopWindowLArSurfxace;
    G4LogicalBorderSurface*     fBottomWindowLArSurface;
    G4MaterialPropertiesTable*  fWLArSurfProp;

    G4OpticalSurface*           fOpWindowTPBSurface;
    G4LogicalBorderSurface*     fTopWindowTPBSurface;
    G4LogicalBorderSurface*     fBottomWindowTPBSurface;
    G4MaterialPropertiesTable*  fWindowTPBSurfProp;

    G4OpticalSurface*           fOpTPBGArSurface;
    G4LogicalBorderSurface*     fTPBGArSurface;
    G4MaterialPropertiesTable*  fTPBGArSurfProp;

    G4OpticalSurface*           fOpTPBLArSurface;
    G4LogicalBorderSurface*     fTPBLArSurface;
    G4MaterialPropertiesTable*  fTPBLArSurfProp;

    G4OpticalSurface*           fOpReflectorTPBSurface;
    G4LogicalBorderSurface*     fReflectorTPBSurface;
    G4MaterialPropertiesTable*  fReflTPBSurfProp;

    G4OpticalSurface*           fOpGArLArSurface;
    G4LogicalBorderSurface*     fGArLArSurface;
    G4MaterialPropertiesTable*  fGArLArSurfProp;

    G4OpticalSurface*           fOpLArTeflonSurface;
    G4LogicalBorderSurface*     fLArReflectorSurface; 
    G4LogicalBorderSurface*     fLArPMTTopAssemblySurface;
    G4LogicalBorderSurface*     fLArPMTBottomAssemblySurface;
    G4LogicalBorderSurface*     fLArTeflonSupportSurface;
    G4MaterialPropertiesTable*  fLArTeflonSurfProp;

    G4OpticalSurface*           fOpElectropolishedStainlessSteelSurface;
    G4LogicalBorderSurface*     fDSG2OuterSurface;

    G4OpticalSurface*           fOpUntreatedStainlessSteelSurface;
    G4LogicalBorderSurface*     fTrunk0Surface;
    G4LogicalBorderSurface*     fTrunk1Surface;
    G4LogicalBorderSurface*     fTrunk2Surface;
    G4LogicalBorderSurface*     fTrunk3Surface;
    G4LogicalBorderSurface*     fTrunk4Surface;
    G4LogicalBorderSurface*     fTrunk5Surface;
    G4LogicalBorderSurface*     fTrunk6Surface;
    G4LogicalBorderSurface*     fTrunkArSurface;
   */
};

#endif
