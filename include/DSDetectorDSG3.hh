#ifndef DSDetectorDSG3_H
#define DSDetectorDSG3_H

#include "G4LogicalVolume.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"


class	DSMaterial;

class   DSDetectorDSG3  {

  public:

    DSDetectorDSG3(G4VPhysicalVolume*);
   ~DSDetectorDSG3();


    G4VPhysicalVolume*          GetDetectorComponent()  { return  fPhysicOuterLiqArgon; }

  private:
    void                        DefineSurfaces();

    G4VPhysicalVolume*          fMotherVolume;
    
    // Cryostats
    G4Polycone*                 fSolidDSOuterCryostat;
    G4LogicalVolume*            fLogicDSOuterCryostat;
    G4VPhysicalVolume*          fPhysicDSOuterCryostat;

    G4Polycone*                 fSolidVacuumCryostat;
    G4LogicalVolume*            fLogicVacuumCryostat;
    G4VPhysicalVolume*          fPhysicVacuumCryostat;

    G4Polycone*                 fSolidInnerCryostat;
    G4LogicalVolume*            fLogicInnerCryostat;
    G4VPhysicalVolume*          fPhysicInnerCryostat;

    G4Polycone*                 fSolidGasArgon;
    G4LogicalVolume*            fLogicGasArgon;
    G4VPhysicalVolume*          fPhysicGasArgon;

    G4Polycone*                 fSolidOuterLiqArgon;
    G4LogicalVolume*            fLogicOuterLiqArgon;
    G4VPhysicalVolume*          fPhysicOuterLiqArgon;

    // PMT Assembly
    G4Tubs*                     fSolidPMTAssembly;
    G4LogicalVolume*            fLogicPMTAssemblyTop;
    G4VPhysicalVolume*          fPhysicPMTAssemblyTop;
    G4LogicalVolume*            fLogicPMTAssemblyBot;
    G4VPhysicalVolume*          fPhysicPMTAssemblyBot;

    G4Tubs*                     fSolidCapDisk;
    G4LogicalVolume*            fLogicCapDiskTop;
    G4LogicalVolume*            fLogicCapDiskBot;
    G4VPhysicalVolume*          fPhysicCapDiskTop;
    G4VPhysicalVolume*          fPhysicCapDiskBot;

    G4Tubs*                     fSolidAnodeWindow;
    G4LogicalVolume*            fLogicAnodeWindow;
    G4VPhysicalVolume*          fPhysicAnodeWindow;

    G4Tubs*                     fSolidFieldRings;
    G4LogicalVolume*            fLogicFieldRings;
    G4VPhysicalVolume*          fPhysicFieldRings;

    G4Tubs*                     fSolidFieldCage;
    G4LogicalVolume*            fLogicFieldCage;
    G4VPhysicalVolume*          fPhysicFieldCage;

    G4Tubs*                     fSolidTPB;
    G4LogicalVolume*            fLogicTPB;
    G4VPhysicalVolume*          fPhysicTPB;

    G4Tubs*                     fSolidGasPocket;
    G4LogicalVolume*            fLogicGasPocket;
    G4VPhysicalVolume*          fPhysicGasPocket;

    G4Tubs*                     fSolidActiveLAr;
    G4LogicalVolume*            fLogicActiveLAr;
    G4VPhysicalVolume*          fPhysicActiveLAr;

    G4Tubs*                     fSolidGrid;
    G4LogicalVolume*            fLogicGrid;
    G4VPhysicalVolume*          fPhysicGrid;

    G4Tubs*                     fSolidCathodeWindow;
    G4LogicalVolume*            fLogicCathodeWindow;
    G4VPhysicalVolume*          fPhysicCathodeWindow;



/*    // Optical surfaces
    G4OpticalSurface*           fOpWindowLArSurface;
    G4LogicalBorderSurface*     fTopWindowLArSurfxace;
    G4LogicalBorderSurface*     fBottomWindowLArSurface;
    G4LogicalBorderSurface*     fLArTopWindowSurfxace;
    G4LogicalBorderSurface*     fLArBottomWindowSurface;
    G4MaterialPropertiesTable*  fWLArSurfProp;

    G4OpticalSurface*           fOpWindowTPBSurface;
    G4LogicalBorderSurface*     fTopWindowTPBSurface;
    G4LogicalBorderSurface*     fBottomWindowTPBSurface;
    G4LogicalBorderSurface*     fTPBTopWindowSurface;
    G4LogicalBorderSurface*     fTPBBottomWindowSurface;
    G4MaterialPropertiesTable*  fWindowTPBSurfProp;

    G4OpticalSurface*           fOpTPBGArSurface;
    G4LogicalBorderSurface*     fTPBGArSurface;
    G4LogicalBorderSurface*     fGArTPBSurface;
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

    G4OpticalSurface*           fOpGridLArSurface;
    G4LogicalBorderSurface*     fGridLArSurface;
    G4MaterialPropertiesTable*  fGridLArSurfProp;

    G4OpticalSurface*           fOpLArGridSurface;
    G4LogicalBorderSurface*     fLArGridSurface;
    G4MaterialPropertiesTable*  fLArGridSurfProp;

    G4OpticalSurface*           fOpLArTeflonSurface;
    G4LogicalBorderSurface*     fLArReflectorSurface; 
    G4LogicalBorderSurface*     fLArPMTTopAssemblySurface;
    G4LogicalBorderSurface*     fLArPMTBottomAssemblySurface;
    G4LogicalBorderSurface*     fLArTeflonSupportSurface;
    G4MaterialPropertiesTable*  fLArTeflonSurfProp;

    G4OpticalSurface*           fOpElectropolishedStainlessSteelSurface;
    G4LogicalBorderSurface*     fDSG3OuterSurface;

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
