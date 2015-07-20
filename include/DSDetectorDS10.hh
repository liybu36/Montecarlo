#ifndef DSDetectorDS10_H
#define DSDetectorDS10_H

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"


class DSDetectorDS10 {

  public:

    DSDetectorDS10(G4VPhysicalVolume*);
   ~DSDetectorDS10(); 


  private:
    void                     DefineSurfaces();


    G4VPhysicalVolume*       fMotherVolume;

    // Dewar
    G4Tubs*                  fSolidDewar;
    G4LogicalVolume*         fLogicDewar;
    G4VPhysicalVolume*       fPhysicDewar;

    G4Tubs*                  fSolidLArBath;
    G4LogicalVolume*         fLogicLArBath;
    G4VPhysicalVolume*       fPhysicLArBath;

    // Rods
    G4Tubs*                  fSolidSupportRod;
    G4LogicalVolume*         fLogicSupportRod;
    G4VPhysicalVolume*       fPhysicSupportRod0;
    G4VPhysicalVolume*       fPhysicSupportRod1;
    G4VPhysicalVolume*       fPhysicSupportRod2;
    G4VPhysicalVolume*       fPhysicSupportRod3;

    // Vessel
    G4Tubs*                  fSolidInnerVesselWall;
    G4LogicalVolume*         fLogicInnerVesselWall;
    G4VPhysicalVolume*       fPhysicInnerVesselWall;

    G4Tubs*                  fSolidInnerVesselWindow;
    G4LogicalVolume*         fLogicInnerVesselWindow;
    G4VPhysicalVolume*       fPhysicInnerVesselWindowTop;
    G4VPhysicalVolume*       fPhysicInnerVesselWindowBottom;

    // Gas Pocket and Active LAr
    G4Tubs*                  fSolidGasPocket;
    G4LogicalVolume*         fLogicGasPocket;
    G4VPhysicalVolume*       fPhysicGasPocket;

    G4Tubs*                  fSolidInnerLAr;
    G4LogicalVolume*         fLogicInnerLAr;
    G4VPhysicalVolume*       fPhysicInnerLAr; 

    // Grid
    G4Tubs*                  fSolidGrid;
    G4LogicalVolume*         fLogicGrid;
    G4VPhysicalVolume*       fPhysicGrid;

    // Active Volume Boundaries
    G4Tubs*                  fSolidThreeMLowerFoil;
    G4Tubs*                  fSolidThreeMMiddleFoil;
    G4Tubs*                  fSolidThreeMTopFoil;
    G4LogicalVolume*         fLogicThreeMLowerFoil; 
    G4LogicalVolume*         fLogicThreeMMiddleFoil;
    G4LogicalVolume*         fLogicThreeMTopFoil;
    G4VPhysicalVolume*       fPhysicThreeMLowerFoil; 
    G4VPhysicalVolume*       fPhysicThreeMMiddleFoil;
    G4VPhysicalVolume*       fPhysicThreeMTopFoil;

    // 3M foil support rings
    G4Tubs*                  fSolidSupportRing;
    G4Tubs*                  fSolidTopLiqSupportRing;
    G4Tubs*                  fSolidTopGasSupportRing;
    G4LogicalVolume*         fLogicSupportRing;
    G4LogicalVolume*         fLogicTopLiqSupportRing;
    G4LogicalVolume*         fLogicTopGasSupportRing;
    G4VPhysicalVolume*       fPhysicSupportRing0;
    G4VPhysicalVolume*       fPhysicSupportRing1;
    G4VPhysicalVolume*       fPhysicSupportRing2;
    G4VPhysicalVolume*       fPhysicTopLiqSupportRing;
    G4VPhysicalVolume*       fPhysicTopGasSupportRing;

    // TPB
    G4Tubs*                  fSolidTPBLowerLateral;
    G4Tubs*                  fSolidTPBMiddleLateral;
    G4Tubs*                  fSolidTPBTopLateral;
    G4Tubs*                  fSolidTPBBases;

    G4LogicalVolume*         fLogicTPBLowerLateral;
    G4LogicalVolume*         fLogicTPBMiddleLateral;
    G4LogicalVolume*         fLogicTPBTopLateral;
    G4LogicalVolume*         fLogicTPBBottomBase;
    G4LogicalVolume*         fLogicTPBTopBase;

    G4VPhysicalVolume*       fPhysicTPBLowerLateral;
    G4VPhysicalVolume*       fPhysicTPBMiddleLateral;
    G4VPhysicalVolume*       fPhysicTPBTopLateral;
    G4VPhysicalVolume*       fPhysicTPBBottomBase;
    G4VPhysicalVolume*       fPhysicTPBTopBase;

    // Field Rings and Kapton cover
    G4Tubs*                  fSolidFieldRing;
    G4LogicalVolume*         fLogicFieldRing;
    G4VPhysicalVolume*       fPhysicFieldRings;

    G4Tubs*                  fSolidKaptonBand;
    G4LogicalVolume*         fLogicKaptonBand;
    G4VPhysicalVolume*       fPhysicKaptonBand;

    // Compression plate
    G4Tubs*                  fSolidCompressionPlateTmp;
    G4Tubs*                  fSolidPMTMold;
    G4SubtractionSolid*      fSolidCompressionPlate[7];
    G4LogicalVolume*         fLogicCompressionPlate;
    G4VPhysicalVolume*       fPhysicCompressionPlateTop;
    G4VPhysicalVolume*       fPhysicCompressionPlateBottom;

    // Rods between the Compression plates
    G4Tubs*                  fSolidSteelRod;
    G4LogicalVolume*         fLogicSteelRod;
    G4VPhysicalVolume*       fPhysicThisSteelRod[12];

    // PMTs
    G4Tubs*                  fSolidPMTWindow;
    G4Tubs*                  fSolidPMTBody;
    G4Tubs*                  fSolidPMTVacuum;
    G4LogicalVolume*         fLogicPMTBody;
    G4LogicalVolume*         fLogicPMTVacuum;
    G4LogicalVolume*         fLogicPMTWindow[14];
    G4VPhysicalVolume*       fPhysicPMTWindow[14];
    G4VPhysicalVolume*       fPhysicPMTBody[14];
    G4VPhysicalVolume*       fPhysicPMTVacuum;
   
    // Bubbler Tubes
    G4Tubs*                  fSolidBubblerTube;
    G4LogicalVolume*         fLogicBubblerTube;
    G4VPhysicalVolume*       fPhysicBubblerTube1;
    G4VPhysicalVolume*       fPhysicBubblerTube2;



};


#endif
