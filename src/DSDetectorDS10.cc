#include <iostream>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4UIcommand.hh"

#include "DSDetectorDS10.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSStorage.hh"



using namespace std;



DSDetectorDS10::DSDetectorDS10(G4VPhysicalVolume* myMotherVolume){

  fMotherVolume = myMotherVolume;

  const double myTwoPi = 2*M_PI*rad;
  G4bool myCheckOverlap = DSStorage::Get()->GetCheckOverlap();

  DSLog(routine) << " Constructing DS10 Geometry " << endlog;



  // Trunk
  // Note: treat the trunk as a unique volume instead of dividing it into sub-volumes




  // Copper Dewar + LAr bath
  G4double myDewarWall_t  = 1 * cm;
  G4double myLArBath_h = 127 * cm;     // h e d sono stati invertiti per coerenza
  G4double myLArBath_d = 45.5 * cm;
  G4double myPorts_h = 4.59 * cm;         // From ports_vol = 8135 cm3 (Ref. M)

  G4ThreeVector myZeros( 0, 0, 0 );                               // TEMPORARY: At the center of the Mother Volume? 
  G4ThreeVector myLArBathPos( 0, 0, -myPorts_h/2. );

  fSolidDewar  = new G4Tubs( "Dewar_Solid", 0, myLArBath_d/2. + myDewarWall_t, myLArBath_h/2. + myPorts_h/2. + myDewarWall_t, 0, myTwoPi );
  fLogicDewar  = new G4LogicalVolume( fSolidDewar, DSMaterial::Get()->GetMetalCopper(), "Dewar_Logic" );
  fPhysicDewar = new G4PVPlacement( 0, myZeros, "Dewar", fLogicDewar, fMotherVolume, false, 0, myCheckOverlap ); 

  fSolidLArBath  = new G4Tubs( "LArBath_Solid", 0, myLArBath_d/2., myLArBath_h/2., 0, myTwoPi );
  fLogicLArBath  = new G4LogicalVolume( fSolidLArBath, DSMaterial::Get()->GetLiquidArgon(), "LArBath_Logic" );
  fPhysicLArBath = new G4PVPlacement( 0, myLArBathPos, "LArBath", fLogicLArBath, fPhysicDewar, false, 0, myCheckOverlap );




  // 4 Steel Rods keeping the TPC in place in the dewar
  G4double mySupportRod_d = 3 * cm;
  G4double mySupportRod_t = 0.15 * cm;
  G4double mySupportRod_h = 22.94 * 2.54 * cm;   // 58.2676 cm              //Ref. K
  G4double myCompressionPlate_d = 36.8 * cm;

  G4double mySupportRodX = ( myCompressionPlate_d/2. - mySupportRod_d/2. ) / sqrt( 2. );
  G4double mySupportRodZ =   myLArBath_h/2. - mySupportRod_h/2.; 
  G4double myOffCenter_x = 2.125 * cm;
  G4ThreeVector myXShift( myOffCenter_x, 0, 0);
  G4ThreeVector mySupportRodPos( mySupportRodX, mySupportRodX, mySupportRodZ );


  fSolidSupportRod   = new G4Tubs( "SupportRod_Solid", mySupportRod_d/2. - mySupportRod_t, mySupportRod_d/2., mySupportRod_h/2., 0, myTwoPi );
  fLogicSupportRod   = new G4LogicalVolume( fSolidSupportRod, DSMaterial::Get()->GetStainlessSteel(), "SupportRod_Logic");  
  fPhysicSupportRod0 = new G4PVPlacement( 0, mySupportRodPos + myXShift, "SupportRod", fLogicSupportRod, fPhysicLArBath, true, 0, myCheckOverlap );
  mySupportRodPos.rotateZ( 90*deg );
  fPhysicSupportRod1 = new G4PVPlacement( 0, mySupportRodPos + myXShift, "SupportRod", fLogicSupportRod, fPhysicLArBath, true, 1, myCheckOverlap );
  mySupportRodPos.rotateZ( 90*deg );
  fPhysicSupportRod2 = new G4PVPlacement( 0, mySupportRodPos + myXShift, "SupportRod", fLogicSupportRod, fPhysicLArBath, true, 2, myCheckOverlap );
  mySupportRodPos.rotateZ( 90*deg );
  fPhysicSupportRod3 = new G4PVPlacement( 0, mySupportRodPos + myXShift, "SupportRod", fLogicSupportRod, fPhysicLArBath, true, 3, myCheckOverlap );



  // Inner Vessel: Cylindrical Acrylic Walls + Two Fused Silica Windows Disks (Top and Bottom)
  G4double myLAr_d = 24.1 * cm;
  G4double myLAr_h = 21.5 * cm;
  G4double myGAr_h = 2.0 * cm;
  G4double myCompressionPlate_h = 1.27 * cm;
  G4double myInnerVesselWindow_t = 1.27 * cm;
  G4double myInnerVesselWindow_d = 33.0 * cm;
  G4double myInnerVesselWall_t = 1.9 * cm;
  G4double myPMTInnerVessel_sep = 0.1 * cm;     // Ref. Light Yield in Dark Side 10

  G4double myInnerVesselWall_r = myLAr_d/2. + myInnerVesselWall_t; 
  G4double myInnerVesselWall_h = myLAr_h + myGAr_h;
  G4double myInnerVesselWall_z = myLArBath_h/2. - mySupportRod_h - myCompressionPlate_h - myPMTInnerVessel_sep - myInnerVesselWindow_t - myInnerVesselWall_h/2.;

  G4ThreeVector myInnerVesselWallPos( myOffCenter_x, 0, myInnerVesselWall_z );
  G4ThreeVector myIVDiskTopPos    = myInnerVesselWallPos + G4ThreeVector( 0, 0, myInnerVesselWall_h/2. + myInnerVesselWindow_t/2. );
  G4ThreeVector myIVDiskBottomPos = myInnerVesselWallPos - G4ThreeVector( 0, 0, myInnerVesselWall_h/2. + myInnerVesselWindow_t/2. );


  fSolidInnerVesselWall  = new G4Tubs( "InnerVesselWall_Solid", 0, myInnerVesselWall_r, myInnerVesselWall_h/2., 0, myTwoPi );
  fLogicInnerVesselWall  = new G4LogicalVolume( fSolidInnerVesselWall, DSMaterial::Get()->GetAcrylic(), "InnerVesselWall_Logic" );
  fPhysicInnerVesselWall = new G4PVPlacement( 0, myInnerVesselWallPos, "InnerVesselWall", fLogicInnerVesselWall, fPhysicLArBath, false, 0, myCheckOverlap );


  fSolidInnerVesselWindow        = new G4Tubs( "InnerVesselWindow_Solid", 0, myInnerVesselWindow_d/2., myInnerVesselWindow_t/2., 0, myTwoPi );
  fLogicInnerVesselWindow        = new G4LogicalVolume( fSolidInnerVesselWindow, DSMaterial::Get()->GetFusedSilica(), "InnerVesselWindow_Logic" );
  fPhysicInnerVesselWindowTop    = new G4PVPlacement( 0, myIVDiskTopPos, "InnerVesselWindowTop",    fLogicInnerVesselWindow, fPhysicLArBath, false, 0, myCheckOverlap );
  fPhysicInnerVesselWindowBottom = new G4PVPlacement( 0, myIVDiskBottomPos, "InnerVesselWindowBottom", fLogicInnerVesselWindow, fPhysicLArBath, false, 0, myCheckOverlap );




  // Gas Pocket and Inner LAr
  G4ThreeVector myGasPocketPos( 0, 0,  myInnerVesselWall_h/2. - myGAr_h/2. );
  G4ThreeVector myInnerLArPos( 0, 0, -myInnerVesselWall_h/2. + myLAr_h/2. );

  fSolidGasPocket  = new G4Tubs( "GasPocket_Solid", 0, myLAr_d/2., myGAr_h/2., 0, myTwoPi );
  fLogicGasPocket  = new G4LogicalVolume( fSolidGasPocket, DSMaterial::Get()->GetGaseousArgon(), "GasPocket_Logic" );
  fPhysicGasPocket = new G4PVPlacement( 0, myGasPocketPos, "GasPocket", fLogicGasPocket, fPhysicInnerVesselWall, false, 0, myCheckOverlap );

  fSolidInnerLAr  = new G4Tubs( "InnerLAr_Solid", 0, myLAr_d/2., myLAr_h/2., 0, myTwoPi );
  fLogicInnerLAr  = new G4LogicalVolume( fSolidInnerLAr, DSMaterial::Get()->GetLiquidArgon(), "InnerLAr_Logic" );
  fPhysicInnerLAr = new G4PVPlacement( 0, myInnerLArPos, "InnerLAr", fLogicInnerLAr, fPhysicInnerVesselWall, false, 111, myCheckOverlap );



  // Set the z coordinate of the LAr - GAr interface, necessary for S2 generation in DSLightX
  G4double myLArGArBoundaryPosZ = ( myGasPocketPos + myInnerVesselWallPos + myLArBathPos ).z() - myGAr_h/2.;
  DSStorage::Get()->SetLArGArBoundaryPosZ( myLArGArBoundaryPosZ );



  // The definition of the LAr region here is needed to set the range cuts for this volume 
  // to a smaller value with respect to the rest of the detector ( see DSPhysicsList::SetCuts() )
  G4Region* fLArRegion = new G4Region("LAr_Logic");
  fLogicInnerLAr->SetRegion( fLArRegion );
  fLArRegion->AddRootLogicalVolume(fLogicInnerLAr);         



  // Grid
  // Optical Trasparency for normally incident light: 89%
  G4double myGrid_h = 0.01 * cm;
  G4double myGridZBelowLArSurface = 0.5 * cm;

  G4ThreeVector myGridPos( 0, 0, myLAr_h/2. - myGridZBelowLArSurface);

  //  Grid is simulated as a dielectric. Hence, instead using steel, we use FusedSilica.
   
  fSolidGrid  = new G4Tubs( "Grid_Solid", 0, myLAr_d/2., myGrid_h/2., 0, myTwoPi );
  //fLogicGrid  = new G4LogicalVolume( fSolidGrid, DSMaterial::Get()->GetFusedSilica(), "Grid_Logic" );
  fLogicGrid  = new G4LogicalVolume( fSolidGrid, DSMaterial::Get()->GetGridSteel(), "Grid_Logic" );
  fPhysicGrid = new G4PVPlacement( 0, myGridPos, "Grid", fLogicGrid, fPhysicInnerLAr, false, 0, myCheckOverlap );



  
  // Active Volume: 3M foils + TPB (simulated as surface properties) + Support Rings
  G4double myThreeMFoil_t = 0.0001 * cm;
  G4double myActiveVol_d  = 21.0 * cm;
  G4double myThreeMLowerFoil_h  = myLAr_h - myGridZBelowLArSurface - myGrid_h/2.;
  G4double myThreeMMiddleFoil_h = myGridZBelowLArSurface - myGrid_h/2.;

  G4ThreeVector myThreeMLowerFoilPos ( 0, 0, -myLAr_h/2. + myThreeMLowerFoil_h/2. );
  G4ThreeVector myThreeMMiddleFoilPos( 0, 0,  myLAr_h/2. - myThreeMMiddleFoil_h/2. );


  fSolidThreeMLowerFoil  = new G4Tubs( "ThreeMLowerFoil_Solid",  myActiveVol_d/2., myActiveVol_d/2. + myThreeMFoil_t, myThreeMLowerFoil_h/2.,  0, myTwoPi );  
  fSolidThreeMMiddleFoil = new G4Tubs( "ThreeMMiddleFoil_Solid", myActiveVol_d/2., myActiveVol_d/2. + myThreeMFoil_t, myThreeMMiddleFoil_h/2., 0, myTwoPi );  
  fSolidThreeMTopFoil    = new G4Tubs( "ThreeMTopFoil_Solid",    myActiveVol_d/2., myActiveVol_d/2. + myThreeMFoil_t, myGAr_h/2.,              0, myTwoPi );  

  fLogicThreeMLowerFoil  = new G4LogicalVolume( fSolidThreeMLowerFoil,  DSMaterial::Get()->GetThreeMFoil(), "ThreeMLowerFoil_Logic");
  fLogicThreeMMiddleFoil = new G4LogicalVolume( fSolidThreeMMiddleFoil, DSMaterial::Get()->GetThreeMFoil(), "ThreeMMiddleFoil_Logic");
  fLogicThreeMTopFoil    = new G4LogicalVolume( fSolidThreeMTopFoil,    DSMaterial::Get()->GetThreeMFoil(), "ThreeMTopFoil_Logic");

  fPhysicThreeMLowerFoil  = new G4PVPlacement( 0, myThreeMLowerFoilPos,  "ThreeMLowerFoil",  fLogicThreeMLowerFoil,  fPhysicInnerLAr,  false, 0, myCheckOverlap );
  fPhysicThreeMMiddleFoil = new G4PVPlacement( 0, myThreeMMiddleFoilPos, "ThreeMMiddleFoil", fLogicThreeMMiddleFoil, fPhysicInnerLAr,  false, 0, myCheckOverlap );
  fPhysicThreeMTopFoil    = new G4PVPlacement( 0, myZeros,               "ThreeMTopFoil",    fLogicThreeMTopFoil,    fPhysicGasPocket, false, 0, myCheckOverlap );


  
  G4double mySupportRing_id    = myActiveVol_d + 2*myThreeMFoil_t;
  G4double mySupportRing_od    = mySupportRing_id + 2.936875 * cm;
  G4double mySupportRing_h     = 0.782 * cm; 
  //G4double myTopSupportRing_od = mySupportRing_id + 3.413125 * cm;
  G4double myTopSupportRing_od = mySupportRing_od;
  G4double myTopSupportRing_h  = 1.3589 * cm;
  //G4double myTopLiqSupportRing_h = myTopSupportRing_h/2.;
  G4double myTopLiqSupportRing_h = myGridZBelowLArSurface - myGrid_h;

  G4ThreeVector mySupportRingPos0 = G4ThreeVector( 0, 0, -myLAr_h/2. + mySupportRing_h/2. );
  G4ThreeVector mySupportRingPos1 = mySupportRingPos0 + G4ThreeVector( 0, 0, mySupportRing_h + 3.0*2.54*cm );
  G4ThreeVector mySupportRingPos2 = mySupportRingPos1 + G4ThreeVector( 0, 0, mySupportRing_h + 4.1*2.54*cm );
  G4ThreeVector myTopLiqSupportRing( 0, 0, myLAr_h/2. - myTopLiqSupportRing_h/2. );
  G4ThreeVector myTopGasSupportRing( 0, 0, myGAr_h/2. - myTopSupportRing_h/2. ); 


  fSolidSupportRing       = new G4Tubs( "SupportRing_Solid",       mySupportRing_id/2., mySupportRing_od/2.,    mySupportRing_h/2.,       0, myTwoPi );
  fSolidTopLiqSupportRing = new G4Tubs( "TopLiqSupportRing_Solid", mySupportRing_id/2., myTopSupportRing_od/2., myTopLiqSupportRing_h/2., 0, myTwoPi );
  fSolidTopGasSupportRing = new G4Tubs( "TopGasSupportRing_Solid", mySupportRing_id/2., myTopSupportRing_od/2., myTopSupportRing_h/2.,    0, myTwoPi );

  fLogicSupportRing       = new G4LogicalVolume( fSolidSupportRing, DSMaterial::Get()->GetTeflon(), "SupportRing_Logic" ); 
  fLogicTopLiqSupportRing = new G4LogicalVolume( fSolidTopLiqSupportRing, DSMaterial::Get()->GetTeflon(), "TopLiqSupportRing_Logic");
  fLogicTopGasSupportRing = new G4LogicalVolume( fSolidTopGasSupportRing, DSMaterial::Get()->GetTeflon(), "TopGasSupportRing_Logic");

  fPhysicSupportRing0      = new G4PVPlacement( 0, mySupportRingPos0, "SupportRingBottom", fLogicSupportRing, fPhysicInnerLAr, false, 0, myCheckOverlap ); 
  fPhysicSupportRing1      = new G4PVPlacement( 0, mySupportRingPos1, "SupportRingMiddle", fLogicSupportRing, fPhysicInnerLAr, false, 0, myCheckOverlap ); 
  fPhysicSupportRing2      = new G4PVPlacement( 0, mySupportRingPos2, "SupportRingTop",    fLogicSupportRing, fPhysicInnerLAr, false, 0, myCheckOverlap ); 
  fPhysicTopLiqSupportRing = new G4PVPlacement( 0, myTopLiqSupportRing, "TopLiqSupportRing", fLogicTopLiqSupportRing, fPhysicInnerLAr, false, 0, myCheckOverlap ); 
  fPhysicTopGasSupportRing = new G4PVPlacement( 0, myTopGasSupportRing, "TopGasSupportRing", fLogicTopGasSupportRing, fPhysicGasPocket, false, 0, myCheckOverlap ); 



  // TPB
  G4double myTPBThickness = 0.1*mm;
   
  G4ThreeVector myTPBBottomBasePos( 0, 0, -myLAr_h/2. + myTPBThickness/2. );
  G4ThreeVector myTPBTopBasePos   ( 0, 0,  myGAr_h/2. - myTPBThickness/2. );

  fSolidTPBLowerLateral  = new G4Tubs( "TPBLowerLateral_Solid", myActiveVol_d/2. - myTPBThickness, myActiveVol_d/2., myThreeMLowerFoil_h/2.,  0, myTwoPi );
  fSolidTPBMiddleLateral = new G4Tubs( "TPBLowerLateral_Solid", myActiveVol_d/2. - myTPBThickness, myActiveVol_d/2., myThreeMMiddleFoil_h/2., 0, myTwoPi );
  fSolidTPBTopLateral    = new G4Tubs( "TPBLowerLateral_Solid", myActiveVol_d/2. - myTPBThickness, myActiveVol_d/2., myGAr_h/2.,              0, myTwoPi );
  fSolidTPBBases         = new G4Tubs( "TPBBases_Solid",        0, myActiveVol_d/2. - myTPBThickness, myTPBThickness/2., 0, myTwoPi ); 

  fLogicTPBLowerLateral  = new G4LogicalVolume( fSolidTPBLowerLateral,  DSMaterial::Get()->GetTPB(), "TPBLowerLateral_Logic" );
  fLogicTPBMiddleLateral = new G4LogicalVolume( fSolidTPBMiddleLateral, DSMaterial::Get()->GetTPB(), "TPBMiddleLateral_Logic" );
  fLogicTPBTopLateral    = new G4LogicalVolume( fSolidTPBTopLateral,    DSMaterial::Get()->GetTPB(), "TPBTopLateral_Logic" );
  fLogicTPBTopBase       = new G4LogicalVolume( fSolidTPBBases,         DSMaterial::Get()->GetTPB(), "TPBTopBase_Logic" );
  fLogicTPBBottomBase    = new G4LogicalVolume( fSolidTPBBases,         DSMaterial::Get()->GetTPB(), "TPBBottomBase_Logic" );

  fPhysicTPBLowerLateral  = new G4PVPlacement( 0, myThreeMLowerFoilPos,  "TPBLowerLateral",  fLogicTPBLowerLateral,  fPhysicInnerLAr,  false, 0, myCheckOverlap ); 
  fPhysicTPBMiddleLateral = new G4PVPlacement( 0, myThreeMMiddleFoilPos, "TPBMiddleLateral", fLogicTPBMiddleLateral, fPhysicInnerLAr,  false, 0, myCheckOverlap ); 
  fPhysicTPBTopLateral    = new G4PVPlacement( 0, myZeros,               "TPBTopLateral",    fLogicTPBTopLateral,    fPhysicGasPocket, false, 0, myCheckOverlap ); 
  fPhysicTPBBottomBase    = new G4PVPlacement( 0, myTPBBottomBasePos,    "TPBBottomBase",    fLogicTPBBottomBase,    fPhysicInnerLAr,  false, 0, myCheckOverlap );
  fPhysicTPBTopBase       = new G4PVPlacement( 0, myTPBTopBasePos,       "TPBTopBase",       fLogicTPBTopBase,       fPhysicGasPocket, false, 0, myCheckOverlap );



  // Field Rings and Kapton Band covering them
  // both are simulated as a unique hollow cylinder
  G4double myCopperRing_t = 0.0068 * mm;               //  2. * 305*g/m2/DSMaterial::Get()->GetMetalCopper()->GetDensity();        //Ref. F  
  G4double myKapton_t     = 0.0076 * mm;               //  3. * 0.001 * 2.54 * cm;                        //Ref. F

  G4double myKapton_rin   = myInnerVesselWall_r + myCopperRing_t;

  fSolidFieldRing   = new G4Tubs( "FieldRings_Solid", myInnerVesselWall_r, myInnerVesselWall_r + myCopperRing_t, myInnerVesselWall_h/2., 0, myTwoPi ); 
  fLogicFieldRing   = new G4LogicalVolume( fSolidFieldRing, DSMaterial::Get()->GetMetalCopper(), "FieldRings_Logic" );
  fPhysicFieldRings = new G4PVPlacement( 0, myInnerVesselWallPos, "FieldRings", fLogicFieldRing, fPhysicLArBath, false, 0, myCheckOverlap ); 

  fSolidKaptonBand  = new G4Tubs( "KaptonBand_Solid", myKapton_rin, myKapton_rin + myKapton_t, myInnerVesselWall_h/2., 0, myTwoPi );
  fLogicKaptonBand  = new G4LogicalVolume( fSolidKaptonBand, DSMaterial::Get()->GetKapton(), "KaptonBand_Logic" );
  fPhysicKaptonBand = new G4PVPlacement( 0, myInnerVesselWallPos, "KaptonBand", fLogicKaptonBand, fPhysicLArBath, false, 0, myCheckOverlap );



  // Compression Plate (StainlessSteel Plates at the top and bottom of the TPC holding it together)
  // the space among the PMTs is filled with 1.3 cm thick PTFE reflectors 
  // and the exposed stainlesssteel areas are covered with 3M foils
  G4double myPMT_d = 7.6 * cm;                                         //  3.0 * 2.54 * cm;   Ref. Light Yield in Dark Side 10 

  G4ThreeVector myComprPlTopPos    = myIVDiskTopPos    + G4ThreeVector( 0, 0, myInnerVesselWindow_t/2. + myPMTInnerVessel_sep + myCompressionPlate_h/2. ); 
  G4ThreeVector myComprPlBottomPos = myIVDiskBottomPos - G4ThreeVector( 0, 0, myInnerVesselWindow_t/2. + myPMTInnerVessel_sep + myCompressionPlate_h/2. ); 

 
  fSolidCompressionPlateTmp = new G4Tubs( "CompressionPlateTmp_Solid", 0, myCompressionPlate_d/2., myCompressionPlate_h/2., 0, myTwoPi );
  fSolidPMTMold             = new G4Tubs( "PMTMold_Solid", 0, myPMT_d/2., myCompressionPlate_h/2., 0, myTwoPi );

  fSolidCompressionPlate[0] = new G4SubtractionSolid( "CompressionPlate0_Solid", fSolidCompressionPlateTmp, fSolidPMTMold, 0, myZeros ); 
  G4ThreeVector mySubtrPMT( myPMT_d, 0, 0 );
  for( G4int i = 1; i < 7; i++ ){
    G4String mySubtrSolidName = "CompressionPlate" + G4UIcommand::ConvertToString( i ) + "_Solid";
    fSolidCompressionPlate[i] = new G4SubtractionSolid( mySubtrSolidName, fSolidCompressionPlate[i-1], fSolidPMTMold, 0, mySubtrPMT );
    mySubtrPMT.rotateZ( 60*deg );
  } 

  fLogicCompressionPlate = new G4LogicalVolume( fSolidCompressionPlate[6], DSMaterial::Get()->GetStainlessSteel(), "CompressionPlate_Logic" );

  fPhysicCompressionPlateTop    = new G4PVPlacement( 0, myComprPlTopPos, "CompressionPlateTop", fLogicCompressionPlate, fPhysicLArBath, false, 0, myCheckOverlap );
  fPhysicCompressionPlateBottom = new G4PVPlacement( 0, myComprPlBottomPos, "CompressionPlateBottom", fLogicCompressionPlate, fPhysicLArBath, false, 0, myCheckOverlap );

  

  // PMTs
  // they are approximated with a cylinder of diameter equal to the photocathode window diameter
  // and formed by three volumes: body (Kovar), photocathode window (FusedSilica), vacuum interior
  G4double myQuartz_t = 7.*70.*g/DSMaterial::Get()->GetQuartz()->GetDensity()/(myTwoPi/2.*(3*myPMT_d/2.)*(3*myPMT_d/2.));
  G4double myPMT_h = 12.3 * cm;                                        // Hamamatsu Specs
  G4double myPMT_t = 0.150 * cm;                                       //Ref. C,D,E,H,I
  G4double myPMT_rin   = myPMT_d/2. - myPMT_t/2.;
  G4double myPMT_rout  = myPMT_d/2.;
  G4double myPMTBody_h = myPMT_h - myQuartz_t;

  G4ThreeVector myPMTVacuumPos( 0, 0, - myPMT_t/2. );
  G4ThreeVector myPMTBodyPos[14];   
  G4ThreeVector myPMTWindowPos[14];

  myPMTWindowPos[0] = myComprPlTopPos    + G4ThreeVector( 0, 0, -myCompressionPlate_h/2. + myQuartz_t/2. );
  myPMTWindowPos[7] = myComprPlBottomPos - G4ThreeVector( 0, 0, -myCompressionPlate_h/2. + myQuartz_t/2. );
  myPMTBodyPos[0] = myPMTWindowPos[0] + G4ThreeVector( 0, 0, myQuartz_t/2. + myPMTBody_h/2. );
  myPMTBodyPos[7] = myPMTWindowPos[7] - G4ThreeVector( 0, 0, myQuartz_t/2. + myPMTBody_h/2. );

  for( G4int i = 1; i < 7; i++ ) { 
    myPMTBodyPos[ i ]     = myPMTBodyPos[ 0 ] + mySubtrPMT;
    myPMTBodyPos[ 7+i ]   = myPMTBodyPos[ 7 ] + mySubtrPMT;
    myPMTWindowPos[ i ]   = myPMTWindowPos[ 0 ] + mySubtrPMT;
    myPMTWindowPos[ 7+i ] = myPMTWindowPos[ 7 ] + mySubtrPMT;
    mySubtrPMT.rotateZ( 60*deg );
  }

  fSolidPMTWindow = new G4Tubs( "PMTWindow_Solid", 0, myPMT_rout, myQuartz_t/2.,               0, myTwoPi );
  fSolidPMTBody   = new G4Tubs( "PMTBody_Solid",   0, myPMT_rout, myPMTBody_h/2.,              0, myTwoPi );
  fSolidPMTVacuum = new G4Tubs( "PMTVacuum_Solid", 0, myPMT_rin,  myPMTBody_h/2. - myPMT_t/2., 0, myTwoPi );
  
  fLogicPMTBody   = new G4LogicalVolume( fSolidPMTBody,   DSMaterial::Get()->GetKovar(),  "PMTBody_Logic" );
  fLogicPMTVacuum = new G4LogicalVolume( fSolidPMTVacuum, DSMaterial::Get()->GetVacuum(), "PMTVacuum_Logic" ); 


  G4RotationMatrix* myRot = 0;
  G4RotationMatrix* rotX180 = new G4RotationMatrix;
  rotX180->rotateX( M_PI * rad );

  for( G4int i = 0; i < 14; i++ ){

    if( i == 7 ) myRot = rotX180;

    G4String myWindowLogicName  = "PMTWindow" + G4UIcommand::ConvertToString( i ) + "_Logic";
    G4String myWindowPhysicName = "TPMT_"     + G4UIcommand::ConvertToString( i );
    G4String myBodyPhysicName   = "PMTBody_"  + G4UIcommand::ConvertToString( i );

    fLogicPMTWindow[i]  = new G4LogicalVolume( fSolidPMTWindow, DSMaterial::Get()->GetBialkali(), myWindowLogicName );
   
    // IMPORTANT: set here the index of the cathode material
    DSStorage::Get()->SetPMTMaterialIndex(fLogicPMTWindow[i]->GetMaterial()->GetIndex());

    fPhysicPMTWindow[i] = new G4PVPlacement( 0,     myPMTWindowPos[i], myWindowPhysicName, fLogicPMTWindow[i], fPhysicLArBath, false, 0, myCheckOverlap );
    fPhysicPMTBody[ i ] = new G4PVPlacement( myRot, myPMTBodyPos[i],   myBodyPhysicName, fLogicPMTBody, fPhysicLArBath, false, 0, myCheckOverlap ); 

  }

  fPhysicPMTVacuum = new G4PVPlacement( 0, myPMTVacuumPos, "PMTVacuum", fLogicPMTVacuum, fPhysicPMTBody[0], false, 0, myCheckOverlap );



  // 12 Steel Rods keeping the TPC together
  G4double mySteelRod_d = 0.95 * cm;
  G4double mySteelRod_h = myInnerVesselWall_h + 2*myInnerVesselWindow_t + 2*myPMTInnerVessel_sep;

  G4ThreeVector myRodDispl( myCompressionPlate_d/2. - mySteelRod_d/2., 0, 0 );  

  fSolidSteelRod = new G4Tubs( "SteelRod", 0, mySteelRod_d/2., mySteelRod_h/2., 0, myTwoPi );
  fLogicSteelRod = new G4LogicalVolume( fSolidSteelRod, DSMaterial::Get()->GetStainlessSteel(), "SteelRod_Logic" );
  
  for( G4int i = 0; i < 12; i++ ){
    G4ThreeVector mySteelRodPos = myInnerVesselWallPos + myRodDispl;  
    fPhysicThisSteelRod[i] = new G4PVPlacement( 0, mySteelRodPos, "SteelRod", fLogicSteelRod, fPhysicLArBath, true, i, myCheckOverlap ); 
    myRodDispl.rotateZ( 30*deg );
  }



  // Bubbler Tubes
  G4double myBubbler_t = 0.3 * cm;
  G4double myBubbler_d = 5 * cm;                                      //Ref. G
  G4double myBubbler_h = 16.5 * cm;                                   //Ref. G

  G4ThreeVector myBubblerTube1Pos = myInnerVesselWallPos - G4ThreeVector( myCompressionPlate_d/2. + myBubbler_d/2.,  myBubbler_d/2., 0 );
  G4ThreeVector myBubblerTube2Pos = myInnerVesselWallPos - G4ThreeVector( myCompressionPlate_d/2. + myBubbler_d/2., -myBubbler_d/2., 0 );


  fSolidBubblerTube  = new G4Tubs( "BubblerTube_Solid", myBubbler_d/2. - myBubbler_t, myBubbler_d/2., myBubbler_h/2., 0, myTwoPi );
  fLogicBubblerTube  = new G4LogicalVolume( fSolidBubblerTube, DSMaterial::Get()->GetTeflon(), "BubblerTube_Logic" );  

  fPhysicBubblerTube1 = new G4PVPlacement( 0, myBubblerTube1Pos, "BubblerTube", fLogicBubblerTube, fPhysicLArBath, true, 0, myCheckOverlap );
  fPhysicBubblerTube2 = new G4PVPlacement( 0, myBubblerTube2Pos, "BubblerTube", fLogicBubblerTube, fPhysicLArBath, true, 1, myCheckOverlap );


  DefineSurfaces();

}


DSDetectorDS10::~DSDetectorDS10(){
  ;
}


void DSDetectorDS10::DefineSurfaces(){


  ////////////////////////////////////////
  // GAr - LAr
  ////////////////////////////////////////
  G4OpticalSurface *fOpGArLArSurface = new G4OpticalSurface("OpGArLArSurface");
  new G4LogicalBorderSurface("GArLArSurface", fPhysicGasPocket, fPhysicInnerLAr, fOpGArLArSurface); 
  fOpGArLArSurface->SetType( dielectric_dielectric );
  fOpGArLArSurface->SetModel( unified );
  fOpGArLArSurface->SetFinish( polished );
  G4MaterialPropertiesTable *fGArLArSurfProp = new G4MaterialPropertiesTable();
  fOpGArLArSurface->SetMaterialPropertiesTable( fGArLArSurfProp );


  ////////////////////////////////////////
  // LAr - StainlessSteel (supportRod + compression plates)
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArSteelSurface = new G4OpticalSurface("OpLArSteelSurface");
  new G4LogicalBorderSurface("LArSteelSurfaceSupp0", fPhysicLArBath,  fPhysicSupportRod0, fOpLArSteelSurface); 
  new G4LogicalBorderSurface("LArSteelSurfaceSupp1", fPhysicLArBath,  fPhysicSupportRod1, fOpLArSteelSurface); 
  new G4LogicalBorderSurface("LArSteelSurfaceSupp2", fPhysicLArBath,  fPhysicSupportRod2, fOpLArSteelSurface); 
  new G4LogicalBorderSurface("LArSteelSurfaceSupp3", fPhysicLArBath,  fPhysicSupportRod3, fOpLArSteelSurface); 
  new G4LogicalBorderSurface("LArSteelSurfacePlateTop", fPhysicLArBath,  fPhysicCompressionPlateTop, fOpLArSteelSurface); 
  new G4LogicalBorderSurface("LArSteelSurfacePlateBot", fPhysicLArBath,  fPhysicCompressionPlateBottom, fOpLArSteelSurface); 
  fOpLArSteelSurface->SetType( dielectric_metal );
  fOpLArSteelSurface->SetModel( unified );
  fOpLArSteelSurface->SetFinish( polished );
  fOpLArSteelSurface->SetSigmaAlpha(0.1);                          
  G4MaterialPropertiesTable *fLArSteelSurfProp = new G4MaterialPropertiesTable();
  G4double LArSteelENE[2] = {0.1*eV, 20.0*eV};
  G4double LArSteelREF[2] = {0.99, 0.99};
  G4double LArSteelEFF[2] = {0.00, 0.00};
  fLArSteelSurfProp->AddProperty("REFLECTIVITY", LArSteelENE, LArSteelREF, 2);			 
  fLArSteelSurfProp->AddProperty("EFFICIENCY",   LArSteelENE, LArSteelEFF, 2);		       
  fOpLArSteelSurface->SetMaterialPropertiesTable( fLArSteelSurfProp );


  ////////////////////////////////////////
  // LAr - GridSteel 
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArGridSurface = new G4OpticalSurface("OpLArGridSurface");
  new G4LogicalBorderSurface("LArGridSurfacePlateBot", fPhysicInnerLAr,  fPhysicGrid, fOpLArGridSurface); 
  fOpLArGridSurface->SetType( dielectric_dielectric );
  fOpLArGridSurface->SetModel( unified );
  fOpLArGridSurface->SetFinish( polished );
  //fOpLArGridSurface->SetSigmaAlpha(0.1);                          
  G4MaterialPropertiesTable *fLArGridSurfProp = new G4MaterialPropertiesTable();
  //G4double LArGridENE[2] = {0.1*eV, 20.0*eV};
  //G4double LArGridREF[2] = {0.11, 0.11};   // 89% absorption
  //G4double LArGridEFF[2] = {0.00, 0.00};
  //fLArGridSurfProp->AddProperty("REFLECTIVITY", LArGridENE, LArGridREF, 2);			 
  //fLArGridSurfProp->AddProperty("EFFICIENCY",   LArGridENE, LArGridEFF, 2);		       
  fOpLArGridSurface->SetMaterialPropertiesTable( fLArGridSurfProp );


  ////////////////////////////////////////
  // LArBath - Copper (dewar)
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArCopperSurface = new G4OpticalSurface("OpLArCopperSurface");
  new G4LogicalBorderSurface("LArCopperSurface", fPhysicLArBath, fPhysicDewar, fOpLArCopperSurface); 
  fOpLArCopperSurface->SetType( dielectric_metal );
  fOpLArCopperSurface->SetModel( unified );
  fOpLArCopperSurface->SetFinish( polished );
  fOpLArCopperSurface->SetSigmaAlpha(0.1);                          
  G4MaterialPropertiesTable *fLArCopperSurfProp = new G4MaterialPropertiesTable();
  G4double LArCopperENE[2] = {0.1*eV, 20.0*eV};
  G4double LArCopperREF[2] = {0.99, 0.99};
  G4double LArCopperEFF[2] = {0.00, 0.00};
  fLArCopperSurfProp->AddProperty("REFLECTIVITY", LArCopperENE, LArCopperREF, 2);			 
  fLArCopperSurfProp->AddProperty("EFFICIENCY",   LArCopperENE, LArCopperEFF, 2);		       
  fOpLArCopperSurface->SetMaterialPropertiesTable( fLArCopperSurfProp );


  ////////////////////////////////////////
  // TPB - GAr
  ////////////////////////////////////////
  G4OpticalSurface *fOpTPBGArSurface = new G4OpticalSurface("OpTPBGArSurface");
  new G4LogicalBorderSurface("TPBGArSurface", fPhysicGasPocket, fPhysicTPBTopLateral, fOpTPBGArSurface );
  new G4LogicalBorderSurface("TPBGArSurface", fPhysicGasPocket, fPhysicTPBTopBase,    fOpTPBGArSurface );
  fOpTPBGArSurface->SetType( dielectric_dielectric );
  fOpTPBGArSurface->SetModel( unified );
  fOpTPBGArSurface->SetFinish( ground );
  G4MaterialPropertiesTable *fTPBGArSurfProp = new G4MaterialPropertiesTable();
  fOpTPBGArSurface->SetMaterialPropertiesTable( fTPBGArSurfProp );


  ////////////////////////////////////////
  // LAr(bath) - Kapton
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArKaptonSurface = new G4OpticalSurface("OpLArKaptonSurface");
  new G4LogicalBorderSurface("LArKaptonSurface", fPhysicLArBath, fPhysicKaptonBand, fOpLArKaptonSurface );
  fOpLArKaptonSurface->SetType( dielectric_metal );
  fOpLArKaptonSurface->SetModel( unified );
  fOpLArKaptonSurface->SetFinish( polished );
  fOpLArKaptonSurface->SetSigmaAlpha(0.1);                          
  G4MaterialPropertiesTable *fLArKaptonSurfProp = new G4MaterialPropertiesTable();
  G4double LArKaptonENE[2] = {0.1*eV, 20.0*eV};
  G4double LArKaptonREF[2] = {0.99, 0.99};
  G4double LArKaptonEFF[2] = {0.00, 0.00};
  fLArKaptonSurfProp->AddProperty("REFLECTIVITY", LArKaptonENE, LArKaptonREF, 2);
  fLArKaptonSurfProp->AddProperty("EFFICIENCY",   LArKaptonENE, LArKaptonEFF, 2);
  fOpLArKaptonSurface->SetMaterialPropertiesTable( fLArKaptonSurfProp );


  ////////////////////////////////////////
  // TPB - Acrylic(inner vessel)
  ////////////////////////////////////////
  G4OpticalSurface *fOpTPBVesselSurface = new G4OpticalSurface("OpTPBVesselSurface");
  new G4LogicalBorderSurface("TPBVesselSurfaceTop",    fPhysicTPBTopBase,       fPhysicInnerVesselWindowTop,    fOpTPBVesselSurface );
  new G4LogicalBorderSurface("TPBVesselSurfaceBot",    fPhysicTPBBottomBase,    fPhysicInnerVesselWindowBottom, fOpTPBVesselSurface );
  new G4LogicalBorderSurface("TPBVesselSurfaceTopLat", fPhysicTPBTopLateral,    fPhysicInnerVesselWall,         fOpTPBVesselSurface );
  new G4LogicalBorderSurface("TPBVesselSurfaceMidLat", fPhysicTPBMiddleLateral, fPhysicInnerVesselWall,         fOpTPBVesselSurface );
  new G4LogicalBorderSurface("TPBVesselSurfaceLowLat", fPhysicTPBLowerLateral,  fPhysicInnerVesselWall,         fOpTPBVesselSurface );
  fOpTPBVesselSurface->SetType( dielectric_dielectric );
  fOpTPBVesselSurface->SetModel( unified );
  fOpTPBVesselSurface->SetFinish( ground );
  G4MaterialPropertiesTable *fTPBVesselSurfProp = new G4MaterialPropertiesTable();
  fOpTPBVesselSurface->SetMaterialPropertiesTable( fTPBVesselSurfProp );


  ////////////////////////////////////////
  // LAr - Acrylic(inner vessel)
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArVesselSurface = new G4OpticalSurface("OpLArVesselSurface");
  new G4LogicalBorderSurface("LArVesselWall",        fPhysicInnerLAr,   fPhysicInnerVesselWall,         fOpLArVesselSurface );
  new G4LogicalBorderSurface("LArVesselSurfaceTop",  fPhysicLArBath,    fPhysicInnerVesselWindowTop,    fOpLArVesselSurface );
  new G4LogicalBorderSurface("LArVesselSurfaceWall", fPhysicLArBath,    fPhysicInnerVesselWindowBottom, fOpLArVesselSurface );
  fOpLArVesselSurface->SetType( dielectric_dielectric );
  fOpLArVesselSurface->SetModel( unified );
  fOpLArVesselSurface->SetFinish( ground );
  G4MaterialPropertiesTable *fLArVesselSurfProp = new G4MaterialPropertiesTable();
  fOpLArVesselSurface->SetMaterialPropertiesTable( fLArVesselSurfProp );


  ////////////////////////////////////////
  // GAr - Acrylic(inner vessel)
  ////////////////////////////////////////
  G4OpticalSurface *fOpGArVesselSurface = new G4OpticalSurface("OpGArVesselSurface");
  new G4LogicalBorderSurface("GArVessel",        fPhysicGasPocket,   fPhysicInnerVesselWall,  fOpGArVesselSurface );
  fOpGArVesselSurface->SetType( dielectric_dielectric );
  fOpGArVesselSurface->SetModel( unified );
  fOpGArVesselSurface->SetFinish( ground );
  G4MaterialPropertiesTable *fGArVesselSurfProp = new G4MaterialPropertiesTable();
  fOpGArVesselSurface->SetMaterialPropertiesTable( fGArVesselSurfProp );


  ////////////////////////////////////////
  // TPB - 3Mfoil (reflector)
  ////////////////////////////////////////
  G4OpticalSurface *fOpTPBThreeMSurface = new G4OpticalSurface("OpTPBThreeMSurface");
  new G4LogicalBorderSurface("TPBThreeMSurfaceTop", fPhysicTPBTopLateral,     fPhysicThreeMTopFoil,    fOpTPBThreeMSurface );
  new G4LogicalBorderSurface("TPBThreeMSurfaceMid", fPhysicTPBMiddleLateral,  fPhysicThreeMMiddleFoil, fOpTPBThreeMSurface );
  new G4LogicalBorderSurface("TPBThreeMSurfaceBot", fPhysicTPBLowerLateral,   fPhysicThreeMLowerFoil,  fOpTPBThreeMSurface );
  fOpTPBThreeMSurface->SetType( dielectric_metal );
  fOpTPBThreeMSurface->SetModel( unified );
  fOpTPBThreeMSurface->SetFinish( ground );
  fOpTPBThreeMSurface->SetSigmaAlpha(0.5);
  G4double TPBThreeMENE[2] = {0.1*eV, 20.0*eV};
  G4double TPBThreeMREF[2] = {0.99, 0.99};
  G4double TPBThreeMEFF[2] = {0.00, 0.00};
  G4MaterialPropertiesTable *fTPBThreeMSurfProp = new G4MaterialPropertiesTable();
  fTPBThreeMSurfProp->AddProperty("REFLECTIVITY", TPBThreeMENE, TPBThreeMREF, 2);			 
  fTPBThreeMSurfProp->AddProperty("EFFICIENCY",   TPBThreeMENE, TPBThreeMEFF, 2);		       
  fOpTPBThreeMSurface->SetMaterialPropertiesTable( fTPBThreeMSurfProp );


  ////////////////////////////////////////
  // LAr - 3Mfoil (reflector)
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArThreeMSurface = new G4OpticalSurface("OpLArThreeMSurface");
  new G4LogicalBorderSurface("LArThreeMSurfaceTop", fPhysicInnerLAr,  fPhysicThreeMTopFoil,    fOpLArThreeMSurface );
  new G4LogicalBorderSurface("LArThreeMSurfaceMid", fPhysicInnerLAr,  fPhysicThreeMMiddleFoil, fOpLArThreeMSurface );
  new G4LogicalBorderSurface("LArThreeMSurfaceBot", fPhysicInnerLAr,  fPhysicThreeMLowerFoil,  fOpLArThreeMSurface );
  fOpLArThreeMSurface->SetType( dielectric_metal );
  fOpLArThreeMSurface->SetModel( unified );
  fOpLArThreeMSurface->SetFinish( ground );
  fOpLArThreeMSurface->SetSigmaAlpha(0.5);
  G4double LArThreeMENE[2] = {0.1*eV, 20.0*eV};
  G4double LArThreeMREF[2] = {0.97, 0.97};
  G4double LArThreeMEFF[2] = {0.00, 0.00};
  G4MaterialPropertiesTable *fLArThreeMSurfProp = new G4MaterialPropertiesTable();
  fLArThreeMSurfProp->AddProperty("REFLECTIVITY", LArThreeMENE, LArThreeMREF, 2);			 
  fLArThreeMSurfProp->AddProperty("EFFICIENCY",   LArThreeMENE, LArThreeMEFF, 2);		       
  fOpLArThreeMSurface->SetMaterialPropertiesTable( fLArThreeMSurfProp );


  ////////////////////////////////////////
  // Teflon - 3Mfoil   [NOT NEEDED ?]
  ////////////////////////////////////////
  G4OpticalSurface *fOpTeflonThreeMSurface = new G4OpticalSurface("OpTeflonThreeMSurface");
  new G4LogicalBorderSurface("TeflonThreeMSurfaceGas", fPhysicThreeMTopFoil, fPhysicTopGasSupportRing, fOpTeflonThreeMSurface );
  new G4LogicalBorderSurface("TeflonThreeMSurfaceLiq", fPhysicThreeMTopFoil, fPhysicTopLiqSupportRing, fOpTeflonThreeMSurface );
  new G4LogicalBorderSurface("TeflonThreeMSurfaceR0",  fPhysicThreeMLowerFoil, fPhysicSupportRing0,    fOpTeflonThreeMSurface );  // da verificare !!!
  new G4LogicalBorderSurface("TeflonThreeMSurfaceR1",  fPhysicThreeMMiddleFoil, fPhysicSupportRing1,   fOpTeflonThreeMSurface );// da verificare !!!
  new G4LogicalBorderSurface("TeflonThreeMSurfaceR2",  fPhysicThreeMTopFoil, fPhysicSupportRing2,      fOpTeflonThreeMSurface ); // da verificare !!!
  fOpTeflonThreeMSurface->SetType( dielectric_metal );
  fOpTeflonThreeMSurface->SetModel( unified );
  fOpTeflonThreeMSurface->SetFinish( ground );
  fOpTeflonThreeMSurface->SetSigmaAlpha(0.1);
  G4double TeflonThreeMENE[2] = {0.1*eV, 20.0*eV};
  G4double TeflonThreeMREF[2] = {0.99, 0.99};
  G4double TeflonThreeMEFF[2] = {0.00, 0.00};
  G4MaterialPropertiesTable *fTeflonThreeMSurfProp = new G4MaterialPropertiesTable();
  fTeflonThreeMSurfProp->AddProperty("REFLECTIVITY", TeflonThreeMENE, TeflonThreeMREF, 2);			 
  fTeflonThreeMSurfProp->AddProperty("EFFICIENCY",   TeflonThreeMENE, TeflonThreeMEFF, 2);		       
  fOpTeflonThreeMSurface->SetMaterialPropertiesTable( fTeflonThreeMSurfProp );


  ////////////////////////////////////////
  // Acrylic - Copper
  ////////////////////////////////////////
  G4OpticalSurface *fOpAcrylicCopperSurface = new G4OpticalSurface("OpAcrylicCopperSurface");
  new G4LogicalBorderSurface("AcrylicCopperSurface", fPhysicInnerVesselWall, fPhysicFieldRings, fOpAcrylicCopperSurface );
  fOpAcrylicCopperSurface->SetType( dielectric_metal );
  fOpAcrylicCopperSurface->SetModel( unified );
  fOpAcrylicCopperSurface->SetFinish( ground );
  fOpAcrylicCopperSurface->SetSigmaAlpha(0.1);
  G4double AcrylicCopperENE[2] = {0.1*eV, 20.0*eV};
  G4double AcrylicCopperREF[2] = {0.99, 0.99};
  G4double AcrylicCopperEFF[2] = {0.00, 0.00};
  G4MaterialPropertiesTable *fAcrylicCopperSurfProp = new G4MaterialPropertiesTable();
  fAcrylicCopperSurfProp->AddProperty("REFLECTIVITY", AcrylicCopperENE, AcrylicCopperREF, 2);			 
  fAcrylicCopperSurfProp->AddProperty("EFFICIENCY",   AcrylicCopperENE, AcrylicCopperEFF, 2);		       
  fOpAcrylicCopperSurface->SetMaterialPropertiesTable( fAcrylicCopperSurfProp );


  ////////////////////////////////////////
  // LAr - Teflon 
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArTeflonSurface     = new G4OpticalSurface("OpLArTeflonSurface");
  new G4LogicalBorderSurface("LArTeflonSurfaceTop",   fPhysicInnerLAr, fPhysicTopLiqSupportRing, fOpLArTeflonSurface); 
  new G4LogicalBorderSurface("LArTeflonSurfaceR0",    fPhysicInnerLAr, fPhysicSupportRing0,      fOpLArTeflonSurface); 
  new G4LogicalBorderSurface("LArTeflonSurfaceR1",    fPhysicInnerLAr, fPhysicSupportRing1,      fOpLArTeflonSurface); 
  new G4LogicalBorderSurface("LArTeflonSurfaceR2",    fPhysicInnerLAr, fPhysicSupportRing2,      fOpLArTeflonSurface); 
  new G4LogicalBorderSurface("LArTeflonSurfaceBub1",  fPhysicLArBath,  fPhysicBubblerTube1,      fOpLArTeflonSurface); 
  new G4LogicalBorderSurface("LArTeflonSurfaceBub2",  fPhysicLArBath,  fPhysicBubblerTube2,      fOpLArTeflonSurface); 
  fOpLArTeflonSurface->SetType( dielectric_metal );
  fOpLArTeflonSurface->SetModel( unified );
  fOpLArTeflonSurface->SetFinish( polished );
  fOpLArTeflonSurface->SetSigmaAlpha(0.1);                          
  G4MaterialPropertiesTable *fLArTeflonSurfProp = new G4MaterialPropertiesTable();
  G4double LArTeflonENE[2] = {0.1*eV, 20.0*eV};
  G4double LArTeflonREF[2] = {0.99, 0.99};
  G4double LArTeflonEFF[2] = {0.00, 0.00};
  fLArTeflonSurfProp->AddProperty("REFLECTIVITY", LArTeflonENE, LArTeflonREF, 2);			 
  fLArTeflonSurfProp->AddProperty("EFFICIENCY",   LArTeflonENE, LArTeflonEFF, 2);		       
  fOpLArTeflonSurface->SetMaterialPropertiesTable( fLArTeflonSurfProp );


  ////////////////////////////////////////
  // GAr - Teflon (reflector)
  ////////////////////////////////////////
  G4OpticalSurface *fOpGArTeflonSurface     = new G4OpticalSurface("OpGArTeflonSurface");
  new G4LogicalBorderSurface("GArTeflonSurfaceTop",   fPhysicGasPocket, fPhysicTopGasSupportRing, fOpGArTeflonSurface); 
  fOpGArTeflonSurface->SetType( dielectric_metal );
  fOpGArTeflonSurface->SetModel( unified );
  fOpGArTeflonSurface->SetFinish( polished );
  fOpGArTeflonSurface->SetSigmaAlpha(0.1);                          
  G4MaterialPropertiesTable *fGArTeflonSurfProp = new G4MaterialPropertiesTable();
  G4double GArTeflonENE[2] = {0.1*eV, 20.0*eV};
  G4double GArTeflonREF[2] = {0.99, 0.99};
  G4double GArTeflonEFF[2] = {0.00, 0.00};
  fGArTeflonSurfProp->AddProperty("REFLECTIVITY", GArTeflonENE, GArTeflonREF, 2);			 
  fGArTeflonSurfProp->AddProperty("EFFICIENCY",   GArTeflonENE, GArTeflonEFF, 2);		       
  fOpGArTeflonSurface->SetMaterialPropertiesTable( fGArTeflonSurfProp );


  ////////////////////////////////////////
  // ITO /////
  ////////////////////////////////////////
			 
  G4MaterialPropertiesTable *fITOSurfProp = new G4MaterialPropertiesTable();
  fITOSurfProp->AddConstProperty("DOITO",1);
  
  
  ////////////////////////////////////////
  // LAr (inner) - TPB
  ////////////////////////////////////////
  G4OpticalSurface *fOpTPBLArSurface = new G4OpticalSurface("OpTBPLArSurface");
  new G4LogicalBorderSurface("TPBLArSurfaceInnerBotBas", fPhysicInnerLAr, fPhysicTPBBottomBase,    fOpTPBLArSurface );
  new G4LogicalBorderSurface("TPBLArSurfaceInnerLowLat", fPhysicInnerLAr, fPhysicTPBLowerLateral,  fOpTPBLArSurface );
  new G4LogicalBorderSurface("TPBLArSurfaceInnerMidLat", fPhysicInnerLAr, fPhysicTPBMiddleLateral, fOpTPBLArSurface );
  new G4LogicalBorderSurface("TPBLArSurfaceInnerTopLat", fPhysicInnerLAr, fPhysicTPBTopLateral,    fOpTPBLArSurface );   // not needed ?
  new G4LogicalBorderSurface("TPBLArSurfaceInnerTopBas", fPhysicInnerLAr, fPhysicTPBTopBase,       fOpTPBLArSurface );  // not needed ?
  fOpTPBLArSurface->SetType( dielectric_dielectric );
  fOpTPBLArSurface->SetModel( unified );
  fOpTPBLArSurface->SetFinish( ground );                           
  fOpTPBLArSurface->SetSigmaAlpha(0.1);                          
  G4MaterialPropertiesTable *fTPBLArSurfProp = new G4MaterialPropertiesTable();
  fOpTPBLArSurface->SetMaterialPropertiesTable( fTPBLArSurfProp );

  
  ////////////////////////////////////////
  // PMT surfaces
  ////////////////////////////////////////
 
  // Photocathode - LAr (bath)
  G4OpticalSurface *fOpPMTLArSurface = new G4OpticalSurface("OpPMTLArSurface");
  G4LogicalBorderSurface*     fPMTLArSurface[14];
  G4MaterialPropertiesTable*  fPMTLArSurfProp;
   for(int i = 0; i < 14; i++)  fPMTLArSurface[i] = new G4LogicalBorderSurface("PMTLArSurface", fPhysicLArBath, fPhysicPMTWindow[i], fOpPMTLArSurface );  
  fOpPMTLArSurface->SetType( dielectric_dielectric );
  fOpPMTLArSurface->SetModel( unified );
  fOpPMTLArSurface->SetFinish( polished );  
  fPMTLArSurfProp = new G4MaterialPropertiesTable();
  fPMTLArSurfProp->AddConstProperty("REFLECTIVITY", 0.2);
  fPMTLArSurfProp->AddConstProperty("EFFICIENCY",   0.0);
  fOpPMTLArSurface->SetMaterialPropertiesTable( fPMTLArSurfProp );
  
  
  // LAr (bath) - Kovar (PMT body) 
  G4OpticalSurface *fOpLArKovarSurface = new G4OpticalSurface("OpLArKovarSurface");
  G4LogicalBorderSurface*     fLArKovarSurface[14];
   for(int i = 0; i < 14; i++)  fLArKovarSurface[i] = new G4LogicalBorderSurface("LArKovarSurface", fPhysicLArBath, fPhysicPMTBody[i], fOpLArKovarSurface );  
  fOpLArKovarSurface->SetType( dielectric_metal );
  fOpLArKovarSurface->SetModel( unified );
  fOpLArKovarSurface->SetFinish( polished );
  fOpLArKovarSurface->SetSigmaAlpha(0.1);                          
  G4MaterialPropertiesTable *fLArKovarSurfProp = new G4MaterialPropertiesTable();
  G4double LArKovarENE[2] = {0.1*eV, 20.0*eV};
  G4double LArKovarREF[2] = {1.00, 1.00};
  G4double LArKovarEFF[2] = {0.00, 0.00};
  fLArKovarSurfProp->AddProperty("REFLECTIVITY", LArKovarENE, LArKovarREF, 2);			 
  fLArKovarSurfProp->AddProperty("EFFICIENCY",   LArKovarENE, LArKovarEFF, 2);		       


 ;
}
/*
 * $Log: DSDetectorDS10.cc,v $
 * Revision 1.2  2014/07/25 14:10:09  perassos
 * Improved handling of the LArGArBoundaryZ variable
 *
 * Revision 1.1  2014/05/07 12:21:01  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.11  2014/01/21 10:50:49  perassos
 * Range cuts set to 1um in LAr and to 1 mm elsewhere
 *
 * Revision 1.10  2013/08/06 13:59:47  dfranco
 * Added 3 variables to store the deposited energy in LAr, scintillator, and water. The last two are not yet implemented. g4rooter has been updated with 3 new variables: tpcene, vetoene, and muene
 *
 * Revision 1.9  2013/08/05 10:56:14  perassos
 * Grid added to DS50; GridSteel Rindex defined; GridSteel set as Grid Material
 *
 * Revision 1.8  2013/06/24 13:05:55  dfranco
 * TPC QE values were filled twice: once in the standard way, the second deconvoluting the reflections. The second has been commented
 *
 * Revision 1.7  2013/06/24 12:45:54  dfranco
 * Improved optical surface properties in DS10
 *
 * Revision 1.6  2013/06/22 10:28:24  dfranco
 * added loogger to DSDetectorDS10.cc
 *
 *
 */
