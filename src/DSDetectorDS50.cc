#include <iostream>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "DSDetectorDS50.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSStorage.hh"
#include "DSParameters.hh"
#include "DSDetectorPMTDS50.hh"

#include "G4DisplacedSolid.hh"

using namespace std;

////////////////////////////////////////////////
//////        Detector Description      ////////
////////////////////////////////////////////////
/*

From the center (0,0,0) to the top
 
LiquidArgon     ActiveLAr          27.537 cm 
GaseousArgon    GasPocket          28.474 cm 
TPB             TPB                28.484 cm 
FusedSilica     BellTop            28.802 cm 
LiquidArgon     InnerLiqArgon      29.304 cm 
LiquidArgon     PMTLArDisk_0       29.404 cm 
Bialkali        TPMT_0             29.7627 cm

Back of the PMT:

Vacuum          PMTHeadVac_0       32.853 cm    
Vacuum          PMTJoinVac_0       32.908 cm    
Vacuum          PMTBodyVac_0       36.365 cm 
Vacuum          PMTTopVac_0        41.592 cm 
Kovar           PMTTop_0           41.647 cm 
LiquidArgon     PMTAssemblyTop     43.528 cm 
LiquidArgon     OuterLiquidArgon   44.4 cm 
GaseousArgon    GaseousArgon       62.141 cm 
StainlessSteel  InnerCryostat      62.6 cm  
Vacuum          VacuumCryostat     63.4 cm  
GaseousArgon    TrunkAr            123.4 cm  
Air             World              126.301 cm
StainlessSteel  Trunk6             126.584 cm
GaseousArgon    TrunkAr            131.406 cm
StainlessSteel  Trunk6             131.689 cm
Air             World              2000 cm  

From the center (0,0,0) to the left

LiquidArgon     ActiveLAr          17.77 cm  
TPB             TPB                17.78 cm  
Teflon          Reflector          20.32 cm  
MetalCopper     FieldRings         20.7367 cm  
LiquidArgon     InnerLiqArgon      21.59 cm  
Teflon          TeflonSupport      23.495 cm  
LiquidArgon     OuterLiquidArgon   25.197 cm  
StainlessSteel  InnerCryostat      25.65 cm  
Vacuum          VacuumCryostat     31.6983 cm  
StainlessSteel  OuterCryostat      32.1 cm  


*/
////////////////////////////////////////////////



DSDetectorDS50::DSDetectorDS50(G4VPhysicalVolume *myMotherVolume) {
  
  fMotherVolume = myMotherVolume;

  const double myTwoPi = 2*M_PI*rad;
  G4bool   myCheckOverlap   = DSStorage::Get()->GetCheckOverlap();

  DSLog(routine) << " Constructing DS50 Geometry" << endlog ;


  G4double myCryostatShiftZ = DSParameters::Get()->GetCryostatShiftZ();
  G4ThreeVector myZeros( 0., 0., 0.);
  

  // ----------------------- //
  // ------   Trunks  ------ //
  // ----------------------- //

  G4double mySSSRadius               = DSParameters::Get()->GetSSSRadius()*mm;
  G4double myTrunkBottomFaceZ        = DSParameters::Get()->GetTrunkShiftZ()*mm;
  G4double myTrunk_d                 = DSParameters::Get()->GetTrunkDiameter()*mm;       
  G4double myTrunk_thickness         = DSParameters::Get()->GetTrunkThickness()*mm;
  G4double myTrunkBotTopZoffset      = DSParameters::Get()->GetTrunkTopBottomOffsetZ()*mm;
  G4double myTrunkBotTopXoffset      = DSParameters::Get()->GetTrunkTopBottomOffsetX()*mm;
  G4double myBottomTrunk_h           = DSParameters::Get()->GetTrunkBottomHeight()*mm;
  G4double myMiddleTrunk_h           = DSParameters::Get()->GetTrunkMiddleHeight()*mm;      
  G4double myTrunkDistToCryostatAxis = DSParameters::Get()->GetTrunkDistToCryostatAxis()*mm; 
  G4double myTopTrunk_h              = ( mySSSRadius - ( myCryostatShiftZ + myTrunkBottomFaceZ + myBottomTrunk_h + myTrunkBotTopZoffset ))*0.90; 


  G4ThreeVector myTrunk_3v[7];
  G4ThreeVector myTrunk_pos( myTrunkDistToCryostatAxis, 0, myCryostatShiftZ + myTrunkBottomFaceZ + myBottomTrunk_h/2. );
  myTrunk_3v[0] = myTrunk_pos;
  myTrunk_3v[1] = myTrunk_pos.rotateZ( myTwoPi/6. );
  myTrunk_3v[2] = myTrunk_pos.rotateZ( myTwoPi/6. );
  myTrunk_3v[3] = myTrunk_pos.rotateZ( myTwoPi/6. );
  myTrunk_3v[4] = myTrunk_pos.rotateZ( myTwoPi/6. );
  myTrunk_3v[5] = myTrunk_pos.rotateZ( myTwoPi/6. );
  myTrunk_3v[6] = G4ThreeVector( 0, 0, myCryostatShiftZ + myTrunkBottomFaceZ + myBottomTrunk_h/2. );

  G4ThreeVector myTrunkShiftMid( myTrunkBotTopXoffset/2., 0, myBottomTrunk_h/2. + myTrunkBotTopZoffset/2. );
  G4ThreeVector myTrunkShiftTop( myTrunkBotTopXoffset, 0, myBottomTrunk_h/2. + myTrunkBotTopZoffset + myTopTrunk_h/2. );
  G4RotationMatrix* rotY45  = new G4RotationMatrix;
  rotY45 ->rotateY( myTwoPi/8.);

  fSolidTrunkBottom   = new G4Tubs( "BottomTrunk_Solid", 0, myTrunk_d/2., myBottomTrunk_h/2., 0, myTwoPi );
  fSolidTrunkBottomAr = new G4Tubs( "BottomTrunkAr_Solid", 0, myTrunk_d/2. - myTrunk_thickness, myBottomTrunk_h/2., 0, myTwoPi );

  fSolidTrunkMiddle   = new G4Tubs( "MiddleTrunk_Solid", 0, myTrunk_d/2., myMiddleTrunk_h/2., 0, myTwoPi );
  fSolidTrunkMiddleAr = new G4Tubs( "MiddleTrunkAr_Solid", 0, myTrunk_d/2. - myTrunk_thickness, myMiddleTrunk_h/2., 0, myTwoPi );

  fSolidTrunkTop   = new G4Tubs( "TopTrunk_Solid", 0, myTrunk_d/2., myTopTrunk_h/2., 0, myTwoPi );
  fSolidTrunkTopAr = new G4Tubs( "TopTrunkAr_Solid", 0, myTrunk_d/2. - myTrunk_thickness, myTopTrunk_h/2., 0, myTwoPi );

  fSolidTrunkBotMid   = new G4UnionSolid( "TrunkBotMid_Solid", fSolidTrunkBottom, fSolidTrunkMiddle, rotY45, myTrunkShiftMid );
  fSolidTrunkBotMidAr = new G4UnionSolid( "TrunkBotMidAr_Solid", fSolidTrunkBottomAr, fSolidTrunkMiddleAr, rotY45, myTrunkShiftMid );

  fSolidTrunk   = new G4UnionSolid( "Trunk_Solid", fSolidTrunkBotMid, fSolidTrunkTop, 0, myTrunkShiftTop );
  fSolidTrunkAr = new G4UnionSolid( "TrunkAr_Solid", fSolidTrunkBotMidAr, fSolidTrunkTopAr, 0, myTrunkShiftTop );

  fLogicTrunk   = new G4LogicalVolume( fSolidTrunk, DSMaterial::Get()->GetStainlessSteel(), "Trunk_Logic");
  fLogicTrunkAr = new G4LogicalVolume( fSolidTrunkAr, DSMaterial::Get()->GetGaseousArgon(), "TrunkAr_Logic");

  fPhysicTrunk0 = new G4PVPlacement( 0, myTrunk_3v[0], "Trunk0", fLogicTrunk, fMotherVolume, true, 0, myCheckOverlap );  
  fPhysicTrunk1 = new G4PVPlacement( 0, myTrunk_3v[1], "Trunk1", fLogicTrunk, fMotherVolume, true, 1, myCheckOverlap );  
  fPhysicTrunk2 = new G4PVPlacement( 0, myTrunk_3v[2], "Trunk2", fLogicTrunk, fMotherVolume, true, 2, myCheckOverlap );  
  fPhysicTrunk3 = new G4PVPlacement( 0, myTrunk_3v[3], "Trunk3", fLogicTrunk, fMotherVolume, true, 3, myCheckOverlap );  
  fPhysicTrunk4 = new G4PVPlacement( 0, myTrunk_3v[4], "Trunk4", fLogicTrunk, fMotherVolume, true, 4, myCheckOverlap );  
  fPhysicTrunk5 = new G4PVPlacement( 0, myTrunk_3v[5], "Trunk5", fLogicTrunk, fMotherVolume, true, 5, myCheckOverlap );  
  fPhysicTrunk6 = new G4PVPlacement( 0, myTrunk_3v[6], "Trunk6", fLogicTrunk, fMotherVolume, true, 6, myCheckOverlap );  
  
  fPhysicTrunkAr = new G4PVPlacement( 0, myZeros, "TrunkAr", fLogicTrunkAr, fPhysicTrunk0, true, 0, myCheckOverlap );  
  


  // --------------------------- //
  // ------   Cryostats   ------ //
  // --------------------------- //

  
  //Load cryostat profiles                   // NOTE:  Each cryostat profile is evaluated in the ref frame with z = 0 at the level of the internal cryostat ear
  G4double myOuterCryostatZ[ 30 ];
  G4double myOuterCryostatRout[ 30 ];

  G4double myVacuumCryostatZ[ 30 ];
  G4double myVacuumCryostatRout[ 30 ];

  G4double myInnerCryostatZ[ 30 ];
  G4double myInnerCryostatRout[ 30 ];

  G4double myGasArgonZ[ 30 ];
  G4double myGasArgonRout[ 30 ];

  G4double myLiqArgonZ[ 30 ];
  G4double myLiqArgonRout[ 30 ];

  G4double myRmin[300];
  
  

  G4int    myNumPointsOuterCryo = 0;
  G4int    myNumPointsVacCryo   = 0;
  G4int    myNumPointsInnerCryo = 0;
  G4int    myNumPointsGasArgon  = 0;
  G4int    myNumPointsLiqArgon  = 0;

  int i = 0;

  for(i = 0; i<300;++i) myRmin[i] = 0. ;
  
  // Outer Cryostat - Outer Surface
  DSIO::Get()->GetStreamDSCryostatProfile() >> myNumPointsOuterCryo;
  for( i = 0; i < myNumPointsOuterCryo; i++ ) { 
    DSIO::Get()->GetStreamDSCryostatProfile() >> myOuterCryostatZ[i]  >> myOuterCryostatRout[i];
    myOuterCryostatZ[i] += myCryostatShiftZ;
  }

  // Outer Cryostat - Inner Surface
  DSIO::Get()->GetStreamDSCryostatProfile() >> myNumPointsVacCryo;
  for( i = 0; i < myNumPointsVacCryo;   i++ ) {  
    DSIO::Get()->GetStreamDSCryostatProfile() >> myVacuumCryostatZ[i] >> myVacuumCryostatRout[i];
    myVacuumCryostatZ[i] += myCryostatShiftZ;
  }

  // Inner Cryostat - Outer Surface
  DSIO::Get()->GetStreamDSCryostatProfile() >> myNumPointsInnerCryo;
  for( i = 0; i < myNumPointsInnerCryo; i++ ) {
    DSIO::Get()->GetStreamDSCryostatProfile() >> myInnerCryostatZ[i]  >> myInnerCryostatRout[i];
    myInnerCryostatZ[i] += myCryostatShiftZ;
  }

  // GAr
  DSIO::Get()->GetStreamDSCryostatProfile() >> myNumPointsGasArgon;
  for( i = 0; i < myNumPointsGasArgon;  i++ ) {   
    DSIO::Get()->GetStreamDSCryostatProfile() >> myGasArgonZ[i]  >> myGasArgonRout[i];
    myGasArgonZ[i] += myCryostatShiftZ;
  }

  // LAr
  DSIO::Get()->GetStreamDSCryostatProfile() >> myNumPointsLiqArgon;
  for( i = 0; i < myNumPointsLiqArgon;  i++ ) {
    DSIO::Get()->GetStreamDSCryostatProfile() >> myLiqArgonZ[i]  >> myLiqArgonRout[i];
    myLiqArgonZ[i] += myCryostatShiftZ;
  }
  
  DSIO::Get()->CloseStreamDSCryostatProfile();




  // Outer Cryostat
  fSolidDS50_0 = new G4Polycone( "OuterCryostat0_Solid", 0, myTwoPi, myNumPointsOuterCryo,myOuterCryostatZ, myRmin, myOuterCryostatRout  );
  fSolidDS50_1 = new G4SubtractionSolid( "OuterCryostat1_Solid", fSolidDS50_0, fSolidTrunk, 0, myTrunk_3v[0] );
  fSolidDS50_2 = new G4SubtractionSolid( "OuterCryostat2_Solid", fSolidDS50_1, fSolidTrunk, 0, myTrunk_3v[1] );
  fSolidDS50_3 = new G4SubtractionSolid( "OuterCryostat3_Solid", fSolidDS50_2, fSolidTrunk, 0, myTrunk_3v[2] );
  fSolidDS50_4 = new G4SubtractionSolid( "OuterCryostat4_Solid", fSolidDS50_3, fSolidTrunk, 0, myTrunk_3v[3] );
  fSolidDS50_5 = new G4SubtractionSolid( "OuterCryostat5_Solid", fSolidDS50_4, fSolidTrunk, 0, myTrunk_3v[4] );
  fSolidDS50_6 = new G4SubtractionSolid( "OuterCryostat6_Solid", fSolidDS50_5, fSolidTrunk, 0, myTrunk_3v[5] );
  fSolidDS50   = new G4SubtractionSolid( "OuterCryostat_Solid",  fSolidDS50_6, fSolidTrunk, 0, myTrunk_3v[6] );
  fLogicDS50   = new G4LogicalVolume( fSolidDS50, DSMaterial::Get()->GetStainlessSteel(), "OuterCryostat_Logic" );
  fPhysicDS50  = new G4PVPlacement( 0, myZeros, "OuterCryostat", fLogicDS50, fMotherVolume, false, 0, myCheckOverlap );  //  No translations applied as the shift in z
                                                                                                                         //  is accounted for when profiles are loaded
   
  // Vacuum Region
  //fSolidVacuumCryostat_0 = new G4Polycone( "VacuumCryostat0_Solid", 0, myTwoPi, myNumPointsVacCryo, myVacuumCryostatRout, myVacuumCryostatZ );
  fSolidVacuumCryostat_0 = new G4Polycone( "VacuumCryostat0_Solid", 0, myTwoPi, myNumPointsVacCryo,  myVacuumCryostatZ, myRmin, myVacuumCryostatRout);
  fSolidVacuumCryostat_1 = new G4SubtractionSolid( "VacuumCryostat1_Solid", fSolidVacuumCryostat_0, fSolidTrunk, 0, myTrunk_3v[0] );
  fSolidVacuumCryostat_2 = new G4SubtractionSolid( "VacuumCryostat2_Solid", fSolidVacuumCryostat_1, fSolidTrunk, 0, myTrunk_3v[1] );
  fSolidVacuumCryostat_3 = new G4SubtractionSolid( "VacuumCryostat3_Solid", fSolidVacuumCryostat_2, fSolidTrunk, 0, myTrunk_3v[2] );
  fSolidVacuumCryostat_4 = new G4SubtractionSolid( "VacuumCryostat4_Solid", fSolidVacuumCryostat_3, fSolidTrunk, 0, myTrunk_3v[3] );
  fSolidVacuumCryostat_5 = new G4SubtractionSolid( "VacuumCryostat5_Solid", fSolidVacuumCryostat_4, fSolidTrunk, 0, myTrunk_3v[4] );
  fSolidVacuumCryostat_6 = new G4SubtractionSolid( "VacuumCryostat6_Solid", fSolidVacuumCryostat_5, fSolidTrunk, 0, myTrunk_3v[5] );
  fSolidVacuumCryostat   = new G4SubtractionSolid( "VacuumCryostat_Solid",  fSolidVacuumCryostat_6, fSolidTrunk, 0, myTrunk_3v[6] );
  fLogicVacuumCryostat   = new G4LogicalVolume( fSolidVacuumCryostat, DSMaterial::Get()->GetVacuum(), "VacuumCryostat_Logic" );
  fPhysicVacuumCryostat  = new G4PVPlacement( 0, myZeros, "VacuumCryostat", fLogicVacuumCryostat, fPhysicDS50, false, 0, myCheckOverlap );

  // Inner Cryostat
  fSolidInnerCryostat  = new G4Polycone( "InnerCryostat_Solid", 0, myTwoPi, myNumPointsInnerCryo, myInnerCryostatZ, myRmin, myInnerCryostatRout );
  fLogicInnerCryostat  = new G4LogicalVolume( fSolidInnerCryostat, DSMaterial::Get()->GetStainlessSteel(), "InnerCryostat_Logic" );
  fPhysicInnerCryostat = new G4PVPlacement( 0, myZeros, "InnerCryostat", fLogicInnerCryostat, fPhysicVacuumCryostat, false, 0, myCheckOverlap );

   // Gaseous Argon
  fSolidGasArgon  = new G4Polycone( "GaseousArgon_Solid", 0, myTwoPi, myNumPointsGasArgon, myGasArgonZ, myRmin,  myGasArgonRout);
  fLogicGasArgon  = new G4LogicalVolume( fSolidGasArgon, DSMaterial::Get()->GetGaseousArgon(), "GaseousArgon_Logic" );
  fPhysicGasArgon = new G4PVPlacement( 0, myZeros, "GaseousArgon", fLogicGasArgon, fPhysicInnerCryostat, false, 0, myCheckOverlap );

  // Liquid Argon
  fSolidOuterLiqArgon  = new G4Polycone( "OuterLiquidArgon_Solid", 0, myTwoPi, myNumPointsLiqArgon, myLiqArgonZ , myRmin, myLiqArgonRout );
  if(DSStorage::Get()->GetIsExternalLArScintillating())
    fLogicOuterLiqArgon  = new G4LogicalVolume( fSolidOuterLiqArgon, DSMaterial::Get()->GetLiquidArgon(), "OuterLiquidArgon_Logic" );
  else 
    fLogicOuterLiqArgon  = new G4LogicalVolume( fSolidOuterLiqArgon, DSMaterial::Get()->GetNSLiquidArgon(), "OuterLiquidArgon_Logic" );
    
  fPhysicOuterLiqArgon = new G4PVPlacement( 0, myZeros, "OuterLiquidArgon", fLogicOuterLiqArgon, fPhysicInnerCryostat, false, 0, myCheckOverlap );





  // --------------------- //
  // ------   TPC   ------ //
  // --------------------- //
   

  G4double myTPCshiftZ                = DSParameters::Get()->GetTPCShiftZ()*mm;               
  G4double myTeflonSupport_d          = DSParameters::Get()->GetTeflonSupportDiameter()*mm;  
  G4double myTeflonSupport_thickness  = DSParameters::Get()->GetTeflonSupportThickness()*mm; 
  G4double myTeflonSupport_h          = DSParameters::Get()->GetTeflonSupportHeight()*mm;

	  //  DSLog(routine) << "myTPCshiftZ= " <<myTPCshiftZ <<endlog;


  G4ThreeVector myTeflonSupportPos( 0, 0, myTPCshiftZ );

 
  G4RotationMatrix* rotZ120 = new G4RotationMatrix;
  rotZ120->rotateZ( 2.*M_PI/3.*rad );
  G4RotationMatrix* rotZ240 = new G4RotationMatrix;
  rotZ240->rotateZ( 4.*M_PI/3.*rad );

  fSolidTeflonSupport   = new G4Tubs( "TeflonSupport_Solid", 0, myTeflonSupport_d/2., myTeflonSupport_h/2., 0, myTwoPi);
  fLogicTeflonSupport   = new G4LogicalVolume( fSolidTeflonSupport, DSMaterial::Get()->GetTeflon(), "TeflonSupport_Logic" );
  fPhysicTeflonSupport  = new G4PVPlacement( 0, myTeflonSupportPos, "TeflonSupport", fLogicTeflonSupport, fPhysicOuterLiqArgon, false, 0, myCheckOverlap ); 

  fSolidInnerLiqArgon   = new G4Tubs( "InnerLiqArgon_Solid", 0, myTeflonSupport_d/2. - myTeflonSupport_thickness, myTeflonSupport_h/2., 0, myTwoPi); 
  fLogicInnerLiqArgon   = new G4LogicalVolume( fSolidInnerLiqArgon, DSMaterial::Get()->GetLiquidArgon(), "InnerLiqArgon_Logic" ); 
  fPhysicInnerLiqArgon  = new G4PVPlacement( 0, myZeros, "InnerLiqArgon", fLogicInnerLiqArgon, fPhysicTeflonSupport, false, 0, myCheckOverlap );

  fSolidInnerLiqArgonBetween   = new G4Tubs( "InnerLiqArgonBetween_Solid", myTeflonSupport_d/2. - myTeflonSupport_thickness, myTeflonSupport_d/2., myTeflonSupport_h/2., 0, 80.3*deg ); 
  fLogicInnerLiqArgonBetween   = new G4LogicalVolume( fSolidInnerLiqArgonBetween, DSMaterial::Get()->GetLiquidArgon(), "InnerLiqArgonBetween_Logic" );
  fPhysicInnerLiqArgonBetween1 = new G4PVPlacement( 0,       myZeros, "InnerLiqArgonBetween1", fLogicInnerLiqArgonBetween, fPhysicTeflonSupport, true, 0, myCheckOverlap );
  fPhysicInnerLiqArgonBetween2 = new G4PVPlacement( rotZ120, myZeros, "InnerLiqArgonBetween2", fLogicInnerLiqArgonBetween, fPhysicTeflonSupport, true, 1, myCheckOverlap );
  fPhysicInnerLiqArgonBetween3 = new G4PVPlacement( rotZ240, myZeros, "InnerLiqArgonBetween3", fLogicInnerLiqArgonBetween, fPhysicTeflonSupport, true, 2, myCheckOverlap );
  



  // PMT Assembly Mother Volume & Teflon Cap Disk
  G4double myPMTAssembly_h = DSParameters::Get()->GetPMTAssemblyHeight()*mm;       
  G4double myTeflonCap_d   = DSParameters::Get()->GetTeflonCapDiameter()*mm;
  G4double myTeflonCap_h   = DSParameters::Get()->GetTeflonCapHeight()*mm;
          
  G4ThreeVector myPMTAssemblyTop( 0, 0,    myTPCshiftZ + myTeflonSupport_h/2. + myPMTAssembly_h/2. );
  G4ThreeVector myPMTAssemblyBottom( 0, 0, myTPCshiftZ - myTeflonSupport_h/2. - myPMTAssembly_h/2. ); 
  G4ThreeVector myTeflonCapPos( 0, 0, -myPMTAssembly_h/2. + myTeflonCap_h/2. );
  DSLog(routine) << "myPMTAssemblyTop= " << myTPCshiftZ+myTeflonSupport_h/2.+myPMTAssembly_h/2. <<endlog;

  G4RotationMatrix* rotX180 = new G4RotationMatrix;
  rotX180->rotateX( M_PI*rad );


  fSolidPMTAssembly         = new G4Tubs( "PMTAssembly_Solid", 0, myTeflonCap_d/2., myPMTAssembly_h/2., 0, myTwoPi );

  fLogicPMTAssemblyTop      = new G4LogicalVolume( fSolidPMTAssembly, DSMaterial::Get()->GetLiquidArgon(), "PMTAssemblyTop_Logic" ); 
  fPhysicPMTAssemblyTop     = new G4PVPlacement( 0, myPMTAssemblyTop,    "PMTAssemblyTop",    fLogicPMTAssemblyTop, fPhysicOuterLiqArgon, false, 0, myCheckOverlap );

  fLogicPMTAssemblyBottom   = new G4LogicalVolume( fSolidPMTAssembly, DSMaterial::Get()->GetLiquidArgon(), "PMTAssemblyBottom_Logic" ); 
  fPhysicPMTAssemblyBottom  = new G4PVPlacement( rotX180, myPMTAssemblyBottom, "PMTAssemblyBottom", fLogicPMTAssemblyBottom, fPhysicOuterLiqArgon, false, 0, myCheckOverlap );



  fSolidTeflonCapDisk   = new G4Tubs( "TeflonCapDisk_Solid", 0, myTeflonCap_d/2., myTeflonCap_h/2., 0, myTwoPi );

  fLogicTeflonCapDiskTop   = new G4LogicalVolume( fSolidTeflonCapDisk, DSMaterial::Get()->GetTeflon(), "TeflonCapDiskTop_Logic" );
  fPhysicTeflonCapDiskTop  = new G4PVPlacement( 0, myTeflonCapPos, "TeflonCapDiskTop", fLogicTeflonCapDiskTop, fPhysicPMTAssemblyTop, false, 0, myCheckOverlap );

  fLogicTeflonCapDiskBottom   = new G4LogicalVolume( fSolidTeflonCapDisk, DSMaterial::Get()->GetTeflon(), "TeflonCapDiskBottom_Logic" );
  fPhysicTeflonCapDiskBottom  = new G4PVPlacement( 0, myTeflonCapPos, "TeflonCapDiskBottom", fLogicTeflonCapDiskBottom, fPhysicPMTAssemblyBottom, false, 0, myCheckOverlap );
  
  // PMTs
  // --- Note: fPhysicInnerLiqArgon is passed as argument just for the definition of the surface between the photocathode and the LAr
  DSDetectorPMTDS50( (G4VPhysicalVolume*) fPhysicPMTAssemblyTop,    (G4VPhysicalVolume*) fPhysicTeflonCapDiskTop );
  DSDetectorPMTDS50( (G4VPhysicalVolume*) fPhysicPMTAssemblyBottom, (G4VPhysicalVolume*) fPhysicTeflonCapDiskBottom );
   
  // Cathode Window + ITO
  G4double myLArBottomLayerThickness = DSParameters::Get()->GetLArBottomLayerThickness()*mm;
  G4double myCathodeWindow_d         = DSParameters::Get()->GetCathodeWindowDiameter()*mm;
  G4double myCathodeWindow_h         = DSParameters::Get()->GetCathodeWindowHeight()*mm; 
  G4double myITOThickness            = DSParameters::Get()->GetITOThickness()*mm;

  //G4ThreeVector myCathodeWinPos( 0, 0, - myTeflonSupport_h/2. + myLArBottomLayerThickness + myITOThickness + myCathodeWindow_h/2. );
  G4ThreeVector myCathodeWinPos( 0, 0, - myTeflonSupport_h/2. + myLArBottomLayerThickness + myCathodeWindow_h/2. );
  G4ThreeVector myITOBottomPos ( 0, 0, - myTeflonSupport_h/2. + myLArBottomLayerThickness + myITOThickness/2. ); 
  G4ThreeVector myITOTopPos    ( 0, 0, - myTeflonSupport_h/2. + myLArBottomLayerThickness + myITOThickness + myCathodeWindow_h + myITOThickness/2. ); 

   
  /// set in z = -189.065
  fSolidCathodeWindow     = new G4Tubs( "CathodeWindow_Solid", 0, myCathodeWindow_d/2., myCathodeWindow_h/2., 0, myTwoPi );
  fLogicCathodeWindow     = new G4LogicalVolume( fSolidCathodeWindow, DSMaterial::Get()->GetFusedSilica(), "CathodeWindow_Logic" );
  fPhysicCathodeWindow    = new G4PVPlacement( 0, myCathodeWinPos, "CathodeWindow", fLogicCathodeWindow, fPhysicInnerLiqArgon, false, 0, myCheckOverlap );

  fSolidITO               = new G4Tubs( "ITO_Solid", 0, myCathodeWindow_d/2., myITOThickness/2., 0, myTwoPi );
  fLogicITO               = new G4LogicalVolume( fSolidITO, DSMaterial::Get()->GetITO(), "ITO_Logic" );
  //fPhysicCathodeITOBottom = new G4PVPlacement( 0, myITOBottomPos, "CathodeWindowITOBottom", fLogicITO, fPhysicInnerLiqArgon, false, 0, myCheckOverlap );
  //fPhysicCathodeITOTop    = new G4PVPlacement( 0, myITOTopPos,    "CathodeWindowITOTop",    fLogicITO, fPhysicInnerLiqArgon, false, 0, myCheckOverlap );


  // Reflector
  G4double myReflector_id = DSParameters::Get()->GetReflectorInnerDiameter()*mm;
  G4double myReflector_od = DSParameters::Get()->GetReflectorOuterDiameter()*mm;
  G4double myReflector_h  = DSParameters::Get()->GetReflectorHeight()*mm;
  G4double myAboveGrid_id = DSParameters::Get()->GetAboveGridInnerDiameter()*mm; 
  G4double myAboveGrid_od = DSParameters::Get()->GetAboveGridOuterDiameter()*mm;
  G4double myAboveGrid_h  = DSParameters::Get()->GetAboveGridHeight()*mm; 
  G4double myGasGap       = DSParameters::Get()->GetGasPocketThickness()*mm;     
  G4double myTPBThickness = DSParameters::Get()->GetTPBThickness()*mm;
  G4double myLArThick     = 0.1 * mm;

  //G4ThreeVector myReflectorPos = myCathodeWinPos + G4ThreeVector( 0, 0, 2*myITOThickness + myCathodeWindow_h/2.);
  G4ThreeVector myReflectorPos = myCathodeWinPos + G4ThreeVector( 0, 0, myCathodeWindow_h/2.);
  G4ThreeVector myGasPocketPos( 0, 0, myReflector_h + myAboveGrid_h - myTPBThickness - myGasGap/2.);
  G4ThreeVector myLArLayerPos( 0, 0, myGasGap/2. - myLArThick/2.);

  G4double myReflectorProf_Z[10];
  G4double myReflectorProf_R[10];
  G4double myReflectorTPBProf_Z[10];
  G4double myReflectorTPBProf_R[10];
  G4double myReflectorArProf_Z[10];
  G4double myReflectorArProf_R[10];

  G4int myNZR = 0;
  myReflectorProf_Z[ myNZR ] = 0;                                     myReflectorProf_R[ myNZR++ ] = 0; 
  myReflectorProf_Z[ myNZR ] = 0;                                     myReflectorProf_R[ myNZR++ ] = myReflector_od/2.; 
  myReflectorProf_Z[ myNZR ] = myReflector_h;                         myReflectorProf_R[ myNZR++ ] = myReflector_od/2.; 
  myReflectorProf_Z[ myNZR ] = myReflector_h;                         myReflectorProf_R[ myNZR++ ] = myAboveGrid_od/2.;
  myReflectorProf_Z[ myNZR ] = myReflector_h + myAboveGrid_h;         myReflectorProf_R[ myNZR++ ] = myAboveGrid_od/2.;
  myReflectorProf_Z[ myNZR ] = myReflector_h + myAboveGrid_h;         myReflectorProf_R[ myNZR++ ] = 0; 

  G4int myNZR_TPB = 0;
  myReflectorTPBProf_Z[ myNZR_TPB ] = 0;                                     myReflectorTPBProf_R[ myNZR_TPB++ ] = 0; 
  myReflectorTPBProf_Z[ myNZR_TPB ] = 0;                                     myReflectorTPBProf_R[ myNZR_TPB++ ] = myReflector_id/2.; 
  myReflectorTPBProf_Z[ myNZR_TPB ] = myReflector_h;                         myReflectorTPBProf_R[ myNZR_TPB++ ] = myReflector_id/2.; 
  myReflectorTPBProf_Z[ myNZR_TPB ] = myReflector_h;                         myReflectorTPBProf_R[ myNZR_TPB++ ] = myAboveGrid_id/2.;
  myReflectorTPBProf_Z[ myNZR_TPB ] = myReflector_h + myAboveGrid_h;         myReflectorTPBProf_R[ myNZR_TPB++ ] = myAboveGrid_id/2.;
  myReflectorTPBProf_Z[ myNZR_TPB ] = myReflector_h + myAboveGrid_h;         myReflectorTPBProf_R[ myNZR_TPB++ ] = 0; 

  G4int myNZR_Ar = 0;
  myReflectorArProf_Z[ myNZR_Ar ] = myTPBThickness;                                      myReflectorArProf_R[ myNZR_Ar++ ] = 0; 
  myReflectorArProf_Z[ myNZR_Ar ] = myTPBThickness;                                      myReflectorArProf_R[ myNZR_Ar++ ] = myReflector_id/2. - myTPBThickness; 
  myReflectorArProf_Z[ myNZR_Ar ] = myReflector_h + myTPBThickness;                      myReflectorArProf_R[ myNZR_Ar++ ] = myReflector_id/2. - myTPBThickness;
  myReflectorArProf_Z[ myNZR_Ar ] = myReflector_h + myTPBThickness;                      myReflectorArProf_R[ myNZR_Ar++ ] = myAboveGrid_id/2. - myTPBThickness;
  myReflectorArProf_Z[ myNZR_Ar ] = myReflector_h + myAboveGrid_h - myTPBThickness;      myReflectorArProf_R[ myNZR_Ar++ ] = myAboveGrid_id/2. - myTPBThickness;
  myReflectorArProf_Z[ myNZR_Ar ] = myReflector_h + myAboveGrid_h - myTPBThickness;      myReflectorArProf_R[ myNZR_Ar++ ] = 0; 


  fSolidReflector   = new G4Polycone( "Reflector_Solid", 0, myTwoPi, myNZR, myReflectorProf_Z, myRmin, myReflectorProf_R);
  fLogicReflector   = new G4LogicalVolume( fSolidReflector, DSMaterial::Get()->GetTeflon(), "Reflector_Logic" );
  fPhysicReflector  = new G4PVPlacement( 0, myReflectorPos, "Reflector", fLogicReflector, fPhysicInnerLiqArgon, false, 0, myCheckOverlap );

  fSolidTPB         = new G4Polycone( "TPB_Solid", 0, myTwoPi, myNZR_TPB,myReflectorTPBProf_Z , myRmin, myReflectorTPBProf_R );
  fLogicTPB         = new G4LogicalVolume( fSolidTPB, DSMaterial::Get()->GetTPB(), "TPB_Logic" );
  fPhysicTPB        = new G4PVPlacement( 0, myZeros, "TPB", fLogicTPB, fPhysicReflector, false, 0, myCheckOverlap );

  fSolidActiveLAr   = new G4Polycone( "ActiveLAr_Solid", 0, myTwoPi, myNZR_Ar, myReflectorArProf_Z, myRmin, myReflectorArProf_R );
  fLogicActiveLAr   = new G4LogicalVolume( fSolidActiveLAr, DSMaterial::Get()->GetLiquidArgon(), "ActiveLAr_Logic" );
  fPhysicActiveLAr  = new G4PVPlacement( 0, myZeros, "ActiveLAr", fLogicActiveLAr, fPhysicTPB, false, 111, myCheckOverlap );

  fSolidGasPocket   = new G4Tubs( "GasPocket_Solid", 0, myAboveGrid_id/2. - myTPBThickness, myGasGap/2., 0, myTwoPi );
  if ( DSParameters::Get()->GetWithGasPocket() ) 
    fLogicGasPocket   = new G4LogicalVolume( fSolidGasPocket, DSMaterial::Get()->GetGaseousArgon(), "GasPocket_Logic" );
  else fLogicGasPocket= new G4LogicalVolume( fSolidGasPocket, DSMaterial::Get()->GetLiquidArgon(), "GasPocket_Logic" );
  fPhysicGasPocket  = new G4PVPlacement( 0, myGasPocketPos, "GasPocket", fLogicGasPocket, fPhysicActiveLAr, false, 0, myCheckOverlap );

  fSolidLArLayer   = new G4Tubs( "LArLayer_Solid", 0, myAboveGrid_id/2. - myTPBThickness, myLArThick/2., 0, myTwoPi );
  fLogicLArLayer   = new G4LogicalVolume( fSolidLArLayer, DSMaterial::Get()->GetPseudoArgon(), "LArLayer_Logic" );
  fPhysicLArLayer  = new G4PVPlacement( 0, myLArLayerPos, "LArLayer", fLogicLArLayer, fPhysicGasPocket, false, 0, myCheckOverlap );

  // Set the z coordinate of the LAr - GAr interface, necessary for S2 generation in DSLightX
  G4double myLArGArBoundaryPosZ = (myTeflonSupportPos + myReflectorPos + myGasPocketPos).z() - myGasGap/2.;
  DSStorage::Get()->SetLArGArBoundaryPosZ( myLArGArBoundaryPosZ );


  // The definition of the active LAr region is needed to set the range cuts for this volume 
  // to a smaller value with respect to the rest of the detector ( see DSPhysicsList::SetCuts() )
  G4Region* fLArRegion = new G4Region("LAr_Logic");
  fLogicActiveLAr->SetRegion( fLArRegion );
  fLArRegion->AddRootLogicalVolume(fLogicActiveLAr);


  // Grid
  G4double myGrid_h = 0.01 * cm;
  G4ThreeVector myGridPosition( 0, 0, myReflectorArProf_Z[ 2 ] + myGrid_h/2. );

  fSolidGrid  = new G4Tubs( "Grid_Solid", 0, myReflectorArProf_R[ 2 ], myGrid_h/2., 0, myTwoPi ); 
  fLogicGrid  = new G4LogicalVolume( fSolidGrid, DSMaterial::Get()->GetGridSteel(), "Grid_Logic" );
  //fLogicGrid  = new G4LogicalVolume( fSolidGrid, DSMaterial::Get()->GetLiquidArgon(), "Grid_Logic" );
  fPhysicGrid = new G4PVPlacement( 0, myGridPosition, "Grid", fLogicGrid, fPhysicActiveLAr, false, 0, myCheckOverlap );



  // Diving Bell
  G4double myDivingBell_h    = DSParameters::Get()->GetDivingBellHeight()*mm;
  G4double myDivingBellTop_h = DSParameters::Get()->GetDivingBellTopHeight()*mm;  
  G4double myDivingBell_od   = DSParameters::Get()->GetDivingBellOuterDiameter()*mm;
  G4double myDivingBell_id   = DSParameters::Get()->GetDivingBellInnerDiameter()*mm;  

  //G4ThreeVector myBellTopPos = myReflectorPos + G4ThreeVector( 0, 0, myReflectorProf_Z[ myNZR - 1 ] +  myITOThickness + myDivingBellTop_h/2. );
  G4ThreeVector myBellTopPos = myReflectorPos + G4ThreeVector( 0, 0, myReflectorProf_Z[ myNZR - 1 ] +  myDivingBellTop_h/2. );
  G4ThreeVector myBellRimPos = myReflectorPos + G4ThreeVector( 0, 0, myReflectorProf_Z[ myNZR - 1 ] - (myDivingBell_h - myDivingBellTop_h)/2. );
 
  G4ThreeVector myBellTopITOTop    = myBellTopPos + G4ThreeVector( 0, 0,   myDivingBellTop_h/2. + myITOThickness/2. ); 
  G4ThreeVector myBellTopITOBottom = myBellTopPos + G4ThreeVector( 0, 0, - myDivingBellTop_h/2. - myITOThickness/2. );
 

  fSolidBellTop   = new G4Tubs( "BellTop_Solid", 0,                  myDivingBell_od/2.,  myDivingBellTop_h/2.,                   0, myTwoPi );
  fSolidBellRim   = new G4Tubs( "BellRim_Solid", myDivingBell_id/2., myDivingBell_od/2., (myDivingBell_h - myDivingBellTop_h)/2., 0, myTwoPi );
  fLogicBellTop   = new G4LogicalVolume( fSolidBellTop, DSMaterial::Get()->GetFusedSilica(), "BellTop_Logic" );
  fLogicBellRim   = new G4LogicalVolume( fSolidBellRim, DSMaterial::Get()->GetFusedSilica(), "BellRim_Logic" );
  fPhysicBellTop  = new G4PVPlacement( 0, myBellTopPos, "BellTop", fLogicBellTop, fPhysicInnerLiqArgon, false, 0, myCheckOverlap );
  fPhysicBellRim  = new G4PVPlacement( 0, myBellRimPos, "BellRim", fLogicBellRim, fPhysicInnerLiqArgon, false, 0, myCheckOverlap );

  //fPhysicBellITOTop    = new G4PVPlacement( 0, myBellTopITOTop,    "BellITOTop",    fLogicITO, fPhysicInnerLiqArgon, false, 0, myCheckOverlap );
  //fPhysicBellITOBottom = new G4PVPlacement( 0, myBellTopITOBottom, "BellITOBottom", fLogicITO, fPhysicInnerLiqArgon, false, 0, myCheckOverlap );



  // Field Rings
  G4double myFieldRings_h         = DSParameters::Get()->GetFieldRingsHeight()*mm;
  G4double myFieldRings_thickness = DSParameters::Get()->GetFieldRingsThickness()*mm;

  G4ThreeVector myFieldRingsPos = myReflectorPos + G4ThreeVector( 0, 0, myFieldRings_h/2.);
 
  fSolidFieldRings   = new G4Tubs( "FieldRings_Solid", myReflector_od/2., myReflector_od/2. + myFieldRings_thickness, myFieldRings_h/2., 0, myTwoPi );  
  fLogicFieldRings   = new G4LogicalVolume( fSolidFieldRings, DSMaterial::Get()->GetMetalCopper(), "FieldRings_Logic" );
  fPhysicFieldRings  = new G4PVPlacement( 0, myFieldRingsPos, "FieldRings", fLogicFieldRings, fPhysicInnerLiqArgon, false, 0, myCheckOverlap );


  DefineSurfaces(); 

}

DSDetectorDS50::~DSDetectorDS50(){
  ; //delete fMessenger;
}

void  DSDetectorDS50::DefineSurfaces() {

  ////////////////////////////////////////
  // GAr -> LAr
  // 4-Jun-2014 P.Meyers: We only want to do the default Fresnel/Snell stuff
  //   and this only defines GAr -> LAr (and not LAr -> GAr), so comment it out.
  ////////////////////////////////////////
  
  //  if (DSParameters::Get()->GetWithGasPocket()) {  
  //    G4OpticalSurface *fOpGArLArSurface = new G4OpticalSurface("OpGArLArSurface");
  //    new G4LogicalBorderSurface("GArLArSurface", fPhysicGasPocket, fPhysicActiveLAr, fOpGArLArSurface); 
  //    fOpGArLArSurface->SetType( dielectric_dielectric );
  //    fOpGArLArSurface->SetModel( unified );
  //    fOpGArLArSurface->SetFinish( polished );
  //    G4MaterialPropertiesTable *fGArLArSurfProp = new G4MaterialPropertiesTable();
  //    fOpGArLArSurface->SetMaterialPropertiesTable( fGArLArSurfProp );
  //  }

  ////////////////////////////////////////
  // LAr -> Grid
  //  Note: in this model, all optical action takes place
  //  on entering the grid.
  ////////////////////////////////////////

  G4OpticalSurface *fOpGridLArSurface = new G4OpticalSurface("OpGridLArSurface");
  new G4LogicalBorderSurface("GridLArSurface", fPhysicActiveLAr, fPhysicGrid, fOpGridLArSurface); 
  fOpGridLArSurface->SetType( dielectric_dielectric );
  fOpGridLArSurface->SetModel( glisur );
  //fOpGridLArSurface->SetModel( unified );
  fOpGridLArSurface->SetFinish( polished );
  //fOpGridLArSurface -> SetPolish(0.1);
  //fOpGridLArSurface->SetFinish( ground );
  //fOpGridLArSurface->SetSigmaAlpha(0.1);
  G4double GridLArENE[4] = {0.1*eV, 8.0*eV, 8.3*eV, 20.0*eV};
  G4double GRUV = DSParameters::Get()->GetLArGridUVRef();
  G4double GRVIS = DSParameters::Get()->GetLArGridVisRef();
  G4double GridLArREF[4] = {GRVIS, GRVIS, GRUV, GRUV};
  G4MaterialPropertiesTable *fGridLArSurfProp = new G4MaterialPropertiesTable();
  if (DSParameters::Get()->GetWithNewGridModel()) 
    // the grid model is described in DSStorage.cc, the surface is treated in G4OpBoundaryProcess.cc
    fGridLArSurfProp->AddConstProperty("DOGRID",1);
  // Now use the following in old and new models.  By G4 convention, "reflectivity" is actually 1-absorption.
  fGridLArSurfProp->AddProperty("REFLECTIVITY", GridLArENE, GridLArREF, 4); 
  //  else { 
  //    fGridLArSurfProp->AddProperty("REFLECTIVITY", GridLArENE, GridLArREF, 4); 
  //    }
  fOpGridLArSurface->SetMaterialPropertiesTable( fGridLArSurfProp );

  ////////////////////////////////////////
  // Grid -> LAr (keeping backward labeling convention)
  //  Note: in this model, all optical action takes place
  //  on entering the grid.  Exit action is to just continue
  //  in straight line.
  ////////////////////////////////////////

  if (DSParameters::Get()->GetWithNewGridModel()) {
    G4OpticalSurface *fOpLArGridSurface = new G4OpticalSurface("OpLArGridSurface");
    new G4LogicalBorderSurface("LArGridSurface", fPhysicGrid, fPhysicActiveLAr, fOpLArGridSurface);

    fOpLArGridSurface->SetType( dielectric_dielectric );
  //  fOpLArGridSurface->SetModel( glisur );
  //  fOpLArGridSurface->SetFinish( polished );
    G4MaterialPropertiesTable *fLArGridSurfProp = new G4MaterialPropertiesTable();
    // the grid model is described in DSStorage.cc, the surface is treated in G4OpBoundaryProcess.cc
    fLArGridSurfProp->AddConstProperty("DOGRIDEXIT",1);
    fOpLArGridSurface->SetMaterialPropertiesTable( fLArGridSurfProp );
  }

  ////////////////////////////////////////
  // TPB <--> GAr and TPB <--> LAr
  // This surface will carry all the diffuse properties of the TPB
  // for both GAr and LAr.
  // Make this bi-directional
  ////////////////////////////////////////
  G4OpticalSurface *fOpTPBGArSurface = new G4OpticalSurface("OpTPBGArSurface");
  // Note: the following two work even when the "gas pocket" is LAr in no-pocket runs.
  new G4LogicalBorderSurface("GArTPBSurface", fPhysicGasPocket, fPhysicTPB, fOpTPBGArSurface );
  new G4LogicalBorderSurface("TPBGArSurface", fPhysicTPB, fPhysicGasPocket, fOpTPBGArSurface );
  new G4LogicalBorderSurface("LArTPBSurface", fPhysicActiveLAr, fPhysicTPB, fOpTPBGArSurface );
  new G4LogicalBorderSurface("TPBLArSurface", fPhysicTPB, fPhysicActiveLAr, fOpTPBGArSurface );
  new G4LogicalBorderSurface("LArLayerTPBSurface", fPhysicLArLayer, fPhysicTPB, fOpTPBGArSurface );
  new G4LogicalBorderSurface("TPBLArLayerSurface", fPhysicTPB, fPhysicLArLayer, fOpTPBGArSurface );
  fOpTPBGArSurface->SetType( dielectric_dielectric );
  fOpTPBGArSurface->SetModel( unified );
  fOpTPBGArSurface->SetFinish( ground );
  fOpTPBGArSurface->SetSigmaAlpha(0.3);

  G4double VISTRAN = DSParameters::Get()->GetArTPBVisTran(); 
  const G4int NUM = 4;
  G4double pp[NUM] = {0.1*eV, 8.0*eV, 8.3*eV, 20.0*eV};  // vis, vis, UV, UV
  G4double specularlobe[NUM] = {0., 0., 0., 0.};         //--
  G4double specularspike[NUM] = {0., 0., 0., 0.};        //----  gives all reflection to Lambertian lobe
  G4double backscatter[NUM] = {0., 0., 0., 0.};          //-- 
  G4double reflectivity[NUM] = {1.0, 1.0, 1.0, 1.0};     //  To set 1-absorption
  G4double transmitivity[NUM] = {VISTRAN, VISTRAN, 1.0, 1.0};    //  To set reflection vs. transmission, overridding Fresnel
                                                         //  For now, no angle dependence.

  G4MaterialPropertiesTable *fTPBGArSurfProp = new G4MaterialPropertiesTable();

  fTPBGArSurfProp -> AddProperty("SPECULARLOBECONSTANT",pp,specularlobe,NUM);
  fTPBGArSurfProp -> AddProperty("SPECULARSPIKECONSTANT",pp,specularspike,NUM);
  fTPBGArSurfProp -> AddProperty("BACKSCATTERCONSTANT",pp,backscatter,NUM);
  fTPBGArSurfProp -> AddProperty("REFLECTIVITY",pp,reflectivity,NUM);
  fTPBGArSurfProp -> AddProperty("TRANSMITTANCE",pp,transmitivity,NUM);

  fTPBGArSurfProp -> AddConstProperty("DOArTPB",1);

  fOpTPBGArSurface->SetMaterialPropertiesTable( fTPBGArSurfProp );

  ////////////////////////////////////////
  // ITO /////
  ////////////////////////////////////////
			 
  G4MaterialPropertiesTable *fITOSurfProp = new G4MaterialPropertiesTable();
  if ( DSParameters::Get()->GetWithITO()  ) fITOSurfProp->AddConstProperty("DOITO",1);
  
  ////////////////////////////////////////
  // BellTop (fusedsilica) <--> TPB and CathodeWindow <--> TPB (both with ITO)
  // In the current model, the diffuse nature of the TPB is handled 
  // entirely at the TPB-GAr/LAr surface, not here.
  // Make this bi-directional.
  ////////////////////////////////////////
  G4OpticalSurface *fOpWindowTPBSurface     = new G4OpticalSurface("OpWindowTPBSurface");
  new G4LogicalBorderSurface("TopWindowTPBSurface",    fPhysicBellTop,  fPhysicTPB, fOpWindowTPBSurface );
  new G4LogicalBorderSurface("BottomWindowTPBSurface", fPhysicCathodeWindow, fPhysicTPB, fOpWindowTPBSurface );
  new G4LogicalBorderSurface("TPBTopWindowSurface",    fPhysicTPB, fPhysicBellTop,  fOpWindowTPBSurface );
  new G4LogicalBorderSurface("TPBBottomWindowSurface", fPhysicTPB, fPhysicCathodeWindow, fOpWindowTPBSurface );
  fOpWindowTPBSurface->SetType( dielectric_dielectric );
  fOpWindowTPBSurface->SetModel( unified );
  //  fOpWindowTPBSurface->SetFinish( ground );
  fOpWindowTPBSurface->SetFinish( polished );
  //G4MaterialPropertiesTable *fWindowTPBSurfProp = new G4MaterialPropertiesTable();
  fOpWindowTPBSurface->SetMaterialPropertiesTable( fITOSurfProp );
  
  
  ////////////////////////////////////////
  // BellTop <--> LAr and CathodeWindow <--> LAr (both with ITO)
  // Make this bi-directional 
  ////////////////////////////////////////
  G4OpticalSurface *fOpWindowLArSurface     = new G4OpticalSurface("OpWindowLArSurface");
  new G4LogicalBorderSurface("TopWindowLArSurface",    fPhysicBellTop,       fPhysicInnerLiqArgon, fOpWindowLArSurface); 
  new G4LogicalBorderSurface("BottomWindowLArSurface", fPhysicCathodeWindow, fPhysicInnerLiqArgon, fOpWindowLArSurface); 
  new G4LogicalBorderSurface("LArTopWindowSurface",    fPhysicInnerLiqArgon, fPhysicBellTop,       fOpWindowLArSurface); 
  new G4LogicalBorderSurface("LArBottomWindowSurface", fPhysicInnerLiqArgon, fPhysicCathodeWindow, fOpWindowLArSurface); 
  fOpWindowLArSurface->SetType( dielectric_dielectric );
  fOpWindowLArSurface->SetModel( unified );
  //  fOpWindowLArSurface->SetFinish( ground );
  fOpWindowLArSurface->SetFinish( polished );
  //G4MaterialPropertiesTable *fWLArSurfProp = new G4MaterialPropertiesTable();
  fOpWindowLArSurface->SetMaterialPropertiesTable( fITOSurfProp );


  ////////////////////////////////////////
  // InnerLAr - LAr ???
  ////////////////////////////////////////

  //////////////////////////// if(DSLogger::GetSeverity() == DSLogger::development
  //    && theStep->GetPostStepPoint()) {
  
  
  ///////////////
  // TPB --> Teflon (Reflector)
  //  Should be no Teflon --> TPB as surface is defined as dielectric_metal.
  ////////////////////////////////////////
  G4double TeflonTPBENE[4] = {0.1*eV, 8.0*eV, 8.3*eV, 20.0*eV};
  G4double TREFUV  = DSParameters::Get()->GetTeflonTPBUVRef();
  G4double TREFVIS = DSParameters::Get()->GetTeflonTPBVisRef();
  G4double TeflonTPBREF[4] = {TREFVIS, TREFVIS ,TREFUV , TREFUV };
    
  G4OpticalSurface *fOpTPBTeflonSurface = new G4OpticalSurface("OpTBPTeflonSurface");
  new G4LogicalBorderSurface("TPBTeflonSurface", fPhysicTPB, fPhysicReflector, fOpTPBTeflonSurface );
  fOpTPBTeflonSurface->SetType( dielectric_metal );
  fOpTPBTeflonSurface->SetModel(unified);

  // PDM: though I can't see how, the following settings are giving the desired Lambertian reflection 
  fOpTPBTeflonSurface->SetFinish(groundfrontpainted);                           
  //fOpTPBTeflonSurface->SetFinish(ground);                           
  fOpTPBTeflonSurface->SetSigmaAlpha(0.1);
                          
  G4MaterialPropertiesTable *fTPBTeflonSurfProp = new G4MaterialPropertiesTable();
  fTPBTeflonSurfProp->AddProperty("REFLECTIVITY", TeflonTPBENE, TeflonTPBREF, 4);			 
  fOpTPBTeflonSurface->SetMaterialPropertiesTable( fTPBTeflonSurfProp );

  ////////////////////////////////////////
  // LiquidArgon --> Teflon (fPhysicTeflonCapDiskBottom and fPhysicTeflonCapDiskTop)
  //  Note: "InnerLiqArgon" is the LAr OUTSIDE the TPC active volume.
  //  Should be no Teflon --> LAr as surface is defined as dielectric_metal.
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArCapDiskSurface = new G4OpticalSurface("OpLArCapDiskSurface");
  new G4LogicalBorderSurface("InnerLArCapDiskTopSurface",fPhysicInnerLiqArgon ,fPhysicTeflonCapDiskBottom , fOpLArCapDiskSurface );
  new G4LogicalBorderSurface("InnerLArCapDiskBottomSurface",fPhysicInnerLiqArgon ,fPhysicTeflonCapDiskTop , fOpLArCapDiskSurface );
  fOpLArCapDiskSurface->SetType( dielectric_metal );
  fOpLArCapDiskSurface->SetModel(unified);

  // PDM: though I can't see how, the following settings are giving the desired Lambertian reflection 
  fOpLArCapDiskSurface->SetFinish(groundfrontpainted);                           
  fOpLArCapDiskSurface->SetSigmaAlpha(0.1);                          
  G4MaterialPropertiesTable *fLArCapDiskSurfProp = new G4MaterialPropertiesTable();

  G4double TeflonLArENE[4] = {0.1*eV, 8.0*eV, 8.3*eV, 20.0*eV};
  G4double TAREFUV  = DSParameters::Get()->GetTeflonLArUVRef();
  G4double TAREFVIS = DSParameters::Get()->GetTeflonLArVisRef();
  G4double TeflonLArREF[4] = {TAREFVIS, TAREFVIS ,TAREFUV , TAREFUV };
  fLArCapDiskSurfProp->AddProperty("REFLECTIVITY", TeflonLArENE, TeflonLArREF, 4);			 
  fOpLArCapDiskSurface->SetMaterialPropertiesTable( fLArCapDiskSurfProp );


  ////////////////////////////////////////
  // TPB - LAr
  //   Now done together with TPB <--> GAr
  ////////////////////////////////////////
  // G4OpticalSurface *fOpTPBLArSurface = new G4OpticalSurface("OpTBPLArSurface");
  // new G4LogicalBorderSurface("TPBLArSurface", fPhysicActiveLAr, fPhysicTPB, fOpTPBLArSurface );
  // fOpTPBLArSurface->SetType( dielectric_dielectric );
  // fOpTPBLArSurface->SetModel( unified );
  //  fOpTPBLArSurface->SetFinish( ground );                           
  //  fOpTPBLArSurface->SetSigmaAlpha(0.3);                          
  //  G4MaterialPropertiesTable *fTPBLArSurfProp = new G4MaterialPropertiesTable();
  //  fOpTPBLArSurface->SetMaterialPropertiesTable( fTPBLArSurfProp );

  
  ////////////////////////////////////////
  // Teflon - LAr
  //  PDM: These refer to EXTERNAL surfaces of the TPC.  I don't think they will be executed
  //  and they haven't been debugged.  (Some are defined only for teflon --> LAr, which shouldn't happen.)
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArTeflonSurface          = new G4OpticalSurface("OpLArTeflonSurface");
  new G4LogicalBorderSurface("LArTeflonSurface", fPhysicReflector, fPhysicInnerLiqArgon, fOpLArTeflonSurface );
  new G4LogicalBorderSurface("LArPMTTopAssemblySurface", fPhysicPMTAssemblyTop, fPhysicInnerLiqArgon, fOpLArTeflonSurface);
  new G4LogicalBorderSurface("LArPMTBottomAssemblySurface", fPhysicPMTAssemblyBottom, fPhysicInnerLiqArgon, fOpLArTeflonSurface);
  new G4LogicalBorderSurface("LArTeflonSupportSurface", fPhysicInnerLiqArgon, fPhysicTeflonSupport, fOpLArTeflonSurface); 
  fOpLArTeflonSurface->SetType( dielectric_metal );
  fOpLArTeflonSurface->SetModel( glisur );
  fOpLArTeflonSurface->SetFinish( ground );
  G4MaterialPropertiesTable *fLArTeflonSurfProp = new G4MaterialPropertiesTable();
  // Now defined above
  //  G4double TeflonLArENE[4] = {0.1*eV, 8.0*eV, 8.3*eV, 20.0*eV};
  //  G4double TAREFUV  = DSParameters::Get()->GetTeflonLArUVRef();
  //  G4double TAREFVIS = DSParameters::Get()->GetTeflonLArVisRef();
  //  G4double TeflonLArREF[4] = {TAREFVIS, TAREFVIS ,TAREFUV , TAREFUV };
  fLArTeflonSurfProp->AddProperty("REFLECTIVITY", TeflonLArENE, TeflonLArREF, 4);			 
  fOpLArTeflonSurface->SetMaterialPropertiesTable( fLArTeflonSurfProp );
  

  
  
  
  
  
  ////////////////////////////////////////
  // TPC - BScint
  ////////////////////////////////////////
  G4OpticalSurface *fOpElectropolishedStainlessSteelSurface = new G4OpticalSurface("OpElectropolishedStainlessSteelSurface");
  fOpElectropolishedStainlessSteelSurface->SetType(dielectric_metal);
  fOpElectropolishedStainlessSteelSurface->SetModel(unified);
  fOpElectropolishedStainlessSteelSurface->SetFinish(groundbackpainted);
  fOpElectropolishedStainlessSteelSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetElectropolishedStainlessSteelMPT()); 
  new G4LogicalBorderSurface("DS50OuterSurface",
			     fMotherVolume,
			     fPhysicDS50,
			     fOpElectropolishedStainlessSteelSurface);


  ////////////////////////////////////////
  // Trunks - BScint
  ////////////////////////////////////////
  G4OpticalSurface *fOpUntreatedStainlessSteelSurface = new G4OpticalSurface("OpUntreatedStainlessSteelSurface");
  fOpUntreatedStainlessSteelSurface->SetType(dielectric_metal);
  fOpUntreatedStainlessSteelSurface->SetModel(unified);
  fOpUntreatedStainlessSteelSurface->SetFinish(groundbackpainted);
  fOpUntreatedStainlessSteelSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetUntreatedStainlessSteelMPT());  
  new G4LogicalBorderSurface("Trunk0Surface",
   			 fMotherVolume,
   			 fPhysicTrunk0,
   			 fOpUntreatedStainlessSteelSurface);
  
  new G4LogicalBorderSurface("Trunk1Surface",
   			 fMotherVolume,
   			 fPhysicTrunk1,
   			 fOpUntreatedStainlessSteelSurface);
  
  new G4LogicalBorderSurface("Trunk2Surface",
   			 fMotherVolume,
   			 fPhysicTrunk2,
   			 fOpUntreatedStainlessSteelSurface);
  
  new G4LogicalBorderSurface("Trunk3Surface",
   			 fMotherVolume,
   			 fPhysicTrunk3,
   			 fOpUntreatedStainlessSteelSurface);
  
  new G4LogicalBorderSurface("Trunk4Surface",
   			 fMotherVolume,
   			 fPhysicTrunk4,
   			 fOpUntreatedStainlessSteelSurface);
  
  new G4LogicalBorderSurface("Trunk5Surface",
   			 fMotherVolume,
   			 fPhysicTrunk5,
   			 fOpUntreatedStainlessSteelSurface);
  
  new G4LogicalBorderSurface("Trunk6Surface",
   			 fMotherVolume,
   			 fPhysicTrunk6,
   			 fOpUntreatedStainlessSteelSurface);
  
  new G4LogicalBorderSurface("TrunkArSurface",
   			 fMotherVolume,
		         fPhysicTrunkAr,
                         fOpUntreatedStainlessSteelSurface);

}  



/*
 * $Log: DSDetectorDS50.cc,v $
 * Revision 1.13  2015/01/17 11:31:52  pagnes
 * PAr model added form optical tuning
 *
 * Revision 1.12  2014/11/20 15:32:06  dfranco
 * added a command to remove scintillation process from liquid argon between TPC and cryostat
 *
 * Revision 1.11  2014/07/25 14:10:09  perassos
 * Improved handling of the LArGArBoundaryZ variable
 *
 * Revision 1.10  2014/07/02 17:07:21  meyers
 * Update model for LAr-TPB and GAr-TPB interface.  Now all Lambertian (reflection and transmission).
 *
 * Revision 1.9  2014/06/27 19:54:24  meyers
 * More updates of GAr-TPB and LAr-TPB surfaces
 *
 * Revision 1.8  2014/06/10 13:33:45  meyers
 * Make fused silica-LAr optical surface bi-directional
 *
 * Revision 1.7  2014/06/09 21:13:19  meyers
 * Minor adjustments to Fused silica --> TPB optical model and make it bi-directional
 *
 * Revision 1.6  2014/06/09 14:08:18  meyers
 * Make GAr-TPB optical surface bi-directional
 *
 * Revision 1.5  2014/06/09 13:45:03  meyers
 * Mods to GAr-->TPB optical surface to give more realistic reflection.  Still doesnt execute for TPB-->GAr.
 *
 * Revision 1.4  2014/06/04 13:39:56  meyers
 * Change fused silica -> LAr (with ITO) from ground to polished
 *
 * Revision 1.3  2014/06/04 13:27:48  meyers
 * Commented out GAr -> LAr surface action -- we want default.  Makes no change.
 *
 * Revision 1.2  2014/06/03 13:31:35  meyers
 * Migrate TPC grid and ITO optics updates to g4ds10
 *
 * Revision 1.1  2014/05/07 12:21:01  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.24  2014/03/19 16:37:27  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.23  2014/03/11 16:49:56  dfranco
 * read optical parameters from external file (DSOptics.dat)
 *
 * Revision 1.22  2014/03/11 09:51:23  dfranco
 * preliminary optical optimization
 *
 * Revision 1.21  2014/01/21 10:50:49  perassos
 * Range cuts set to 1um in LAr and to 1 mm elsewhere
 *
 * Revision 1.20  2013/08/06 13:59:47  dfranco
 * Added 3 variables to store the deposited energy in LAr, scintillator, and water. The last two are not yet implemented. g4rooter has been updated with 3 new variables: tpcene, vetoene, and muene
 *
 * Revision 1.19  2013/08/05 10:56:14  perassos
 * Grid added to DS50; GridSteel Rindex defined; GridSteel set as Grid Material
 *
 * Revision 1.18  2013/06/11 22:48:33  dfranco
 * Added ITO optical boundary. The Sernelius function is defined in DSStorage, and called by G4OpBoundaryProcess. New ITO bool variable added to DSDetectorDS50.cc as surface property (G4MaterialPropertyTable)
 *
 * Revision 1.17  2013/06/10 13:55:37  dfranco
 * preliminary (but not tuned) working optics for DS50
 *
 * Revision 1.16  2013/06/05 23:03:32  swesterd
 * moved optical boundary MPTs to DSMaterial and gave the trunks optical boundary properties consistent with untreated stainless steel
 *
 * Revision 1.15  2013/05/31 15:40:06  dfranco
 * fixed some optical properties of the TPC
 *
 * Revision 1.14  2013/05/30 12:34:53  dfranco
 * Fixing the optical properties of the TPC. Not yet concluded
 *
 * Revision 1.13  2013/05/29 16:40:01  dfranco
 * changed optical properties of the materials in order to have consistent array lengths within each material
 *
 * Revision 1.12  2013/05/29 10:50:56  dfranco
 * Restore the correct TPB material
 *
 * Revision 1.11  2013/05/26 03:23:02  swesterd
 * added the data file for veto PMTs and optical boundary properties for the outside of the cryostat
 *
 * Revision 1.10  2013/05/06 14:59:51  perassos
 * Updates on the TPC surface properties and geometry
 *
 * Revision 1.9  2013/05/02 12:41:21  perassos
 * A Few Changes to the TPC Optical Processes
 *
 * Revision 1.8  2013/04/30 15:01:54  perassos
 * Debug DSDetectorDS50.cc
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
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
