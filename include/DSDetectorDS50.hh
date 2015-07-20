#ifndef DSDetectorDS50_H
#define DSDetectorDS50_H

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

class   DSDetectorDS50  {

  public:

    DSDetectorDS50(G4VPhysicalVolume*);
   ~DSDetectorDS50();


    G4VPhysicalVolume*          GetDetectorComponent()  { return  fPhysicOuterLiqArgon; }

  private:
    void                        DefineSurfaces();

    G4VPhysicalVolume*          fMotherVolume;
    
    // Trunks
    G4Tubs*                     fSolidTrunkBottom;
    G4Tubs*                     fSolidTrunkBottomAr;
    G4Tubs*                     fSolidTrunkMiddle;
    G4Tubs*                     fSolidTrunkMiddleAr;
    G4Tubs*                     fSolidTrunkTop;
    G4Tubs*                     fSolidTrunkTopAr;
    G4UnionSolid*               fSolidTrunkBotMid;
    G4UnionSolid*               fSolidTrunkBotMidAr;
    G4UnionSolid*               fSolidTrunk;
    G4UnionSolid*               fSolidTrunkAr;
    G4LogicalVolume*            fLogicTrunk;
    G4LogicalVolume*            fLogicTrunkAr;
    G4VPhysicalVolume*          fPhysicTrunk0;
    G4VPhysicalVolume*          fPhysicTrunk1;
    G4VPhysicalVolume*          fPhysicTrunk2;
    G4VPhysicalVolume*          fPhysicTrunk3;
    G4VPhysicalVolume*          fPhysicTrunk4;
    G4VPhysicalVolume*          fPhysicTrunk5;
    G4VPhysicalVolume*          fPhysicTrunk6;
    G4VPhysicalVolume*          fPhysicTrunkAr;

    // Cryostats
    G4Polycone*                 fSolidDS50_0;
    G4SubtractionSolid*         fSolidDS50_1;
    G4SubtractionSolid*         fSolidDS50_2;
    G4SubtractionSolid*         fSolidDS50_3;
    G4SubtractionSolid*         fSolidDS50_4;
    G4SubtractionSolid*         fSolidDS50_5;
    G4SubtractionSolid*         fSolidDS50_6;
    G4SubtractionSolid*         fSolidDS50;
    G4LogicalVolume*            fLogicDS50;
    G4VPhysicalVolume*          fPhysicDS50;

    G4Polycone*                 fSolidVacuumCryostat_0;
    G4SubtractionSolid*         fSolidVacuumCryostat_1;
    G4SubtractionSolid*         fSolidVacuumCryostat_2;
    G4SubtractionSolid*         fSolidVacuumCryostat_3;
    G4SubtractionSolid*         fSolidVacuumCryostat_4;
    G4SubtractionSolid*         fSolidVacuumCryostat_5;
    G4SubtractionSolid*         fSolidVacuumCryostat_6;
    G4SubtractionSolid*         fSolidVacuumCryostat;
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

    // TPC Teflon Support
    G4Tubs*                     fSolidTeflonSupport;
    G4LogicalVolume*            fLogicTeflonSupport;
    G4VPhysicalVolume*          fPhysicTeflonSupport;

    G4Tubs*                     fSolidInnerLiqArgon;
    G4LogicalVolume*            fLogicInnerLiqArgon;
    G4VPhysicalVolume*          fPhysicInnerLiqArgon;

    G4Tubs*                     fSolidInnerLiqArgonBetween;
    G4LogicalVolume*            fLogicInnerLiqArgonBetween;
    G4VPhysicalVolume*          fPhysicInnerLiqArgonBetween1;
    G4VPhysicalVolume*          fPhysicInnerLiqArgonBetween2;
    G4VPhysicalVolume*          fPhysicInnerLiqArgonBetween3;

    // PMT Assembly
    G4Tubs*                     fSolidPMTAssembly;
    G4LogicalVolume*            fLogicPMTAssemblyTop;
    G4VPhysicalVolume*          fPhysicPMTAssemblyTop;
    G4LogicalVolume*            fLogicPMTAssemblyBottom;
    G4VPhysicalVolume*          fPhysicPMTAssemblyBottom;

    G4Tubs*                     fSolidTeflonCapDisk;
    G4LogicalVolume*            fLogicTeflonCapDiskTop;
    G4LogicalVolume*            fLogicTeflonCapDiskBottom;
    G4VPhysicalVolume*          fPhysicTeflonCapDiskTop;
    G4VPhysicalVolume*          fPhysicTeflonCapDiskBottom;


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

    // LArLayer
    G4Tubs*                     fSolidLArLayer;
    G4LogicalVolume*            fLogicLArLayer;
    G4VPhysicalVolume*          fPhysicLArLayer;

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
    G4LogicalBorderSurface*     fDS50OuterSurface;

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
/*
 * $Log: DSDetectorDS50.hh,v $
 * Revision 1.6  2015/01/17 11:31:46  pagnes
 * PAr model added form optical tuning
 *
 * Revision 1.5  2014/06/10 13:33:39  meyers
 * Make fused silica-LAr optical surface bi-directional
 *
 * Revision 1.4  2014/06/09 21:13:15  meyers
 * Minor adjustments to Fused silica --> TPB optical model and make it bi-directional
 *
 * Revision 1.3  2014/06/09 14:08:15  meyers
 * Make GAr-TPB optical surface bi-directional
 *
 * Revision 1.2  2014/06/09 13:42:27  meyers
 * Add declarations for previous mod -- not really needed, I think
 *
 * Revision 1.1  2014/05/07 12:20:51  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.15  2013/08/05 10:56:16  perassos
 * Grid added to DS50; GridSteel Rindex defined; GridSteel set as Grid Material
 *
 * Revision 1.14  2013/06/10 13:55:40  dfranco
 * preliminary (but not tuned) working optics for DS50
 *
 * Revision 1.13  2013/06/05 23:03:28  swesterd
 * moved optical boundary MPTs to DSMaterial and gave the trunks optical boundary properties consistent with untreated stainless steel
 *
 * Revision 1.12  2013/05/30 12:34:56  dfranco
 * Fixing the optical properties of the TPC. Not yet concluded
 *
 * Revision 1.11  2013/05/27 23:59:00  swesterd
 * added a (currently commented out) Lumirror sheath to the cryostat and introduced DSOpBoundaryProcess to try to figure out why the boundaries are being screwy, with some edits so that it can handle constant and vector properties with freaking out
 *
 * Revision 1.10  2013/05/26 03:23:01  swesterd
 * added the data file for veto PMTs and optical boundary properties for the outside of the cryostat
 *
 * Revision 1.9  2013/05/06 14:59:52  perassos
 * Updates on the TPC surface properties and geometry
 *
 * Revision 1.8  2013/05/02 12:44:01  perassos
 * A Few Changes to the TPC Optical Processes
 *
 * Revision 1.7  2013/04/30 14:44:17  perassos
 * Added Boundaries Optical Properties - still to be defined correctly
 *
 * Revision 1.6  2013/04/26 09:32:19  dfranco
 * Added DefineSurfaces method
 *
 * Revision 1.5  2013/04/19 16:10:09  perassos
 * Added ITO and TPB layers
 *
 * Revision 1.4  2013/03/26 17:13:58  perassos
 * PMT DS50 Updated
 *
 * Revision 1.3  2013/03/22 14:09:39  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
