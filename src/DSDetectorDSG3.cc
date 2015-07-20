#include <iostream>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "DSDetectorDSG3.hh"
#include "DSDetectorPMTDSG3.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSStorage.hh"
#include "DSParameters.hh"
//#include "DSDetectorPMTDSG3.hh"


using namespace std;



///////////////////////////////////////////////////////////////////
//                                                               //
//   For a description of the G3 detector, refer to docdb #944   //
//                                                               //
///////////////////////////////////////////////////////////////////



DSDetectorDSG3::DSDetectorDSG3(G4VPhysicalVolume *myMotherVolume) {
  
  fMotherVolume = myMotherVolume;

  const double myTwoPi = 2*M_PI*rad;
  G4bool   myCheckOverlap   = DSStorage::Get()->GetCheckOverlap();

  DSLog(routine) << " Constructing DSG3 Geometry" << endlog ;

  G4double myCryostatShiftZ = 172.*cm; 
  G4ThreeVector myZeros( 0., 0., 0.);
  

  // ----------------------- //
  // ------   Trunks  ------ //
  // ----------------------- //

  // To be added...

  // --------------------------- //
  // ------   Cryostats   ------ //
  // --------------------------- //

  
  //  Load cryostat profiles                   
  //  
  //  Note: profiles in DSG3CryostatProfilesFromTechDrawings.dat are extracted from technical drawings
  //  assuming as reference z (z = 0) the plane containing the junction between the outer cryostat top and bottom.
  //  From tech drawings, the inner cryostat top/bottom junction is located at about z = -7.5 cm (this offset is already
  //  included in the cryostats profiles).
  //  Profiles are expressed in cm.

  G4double myOuterCryostatZ[ 60 ];
  G4double myOuterCryostatRout[ 60 ];

  G4double myVacuumCryostatZ[ 60 ];
  G4double myVacuumCryostatRout[ 60 ];

  G4double myInnerCryostatZ[ 60 ];
  G4double myInnerCryostatRout[ 60 ];

  G4double myGasArgonZ[ 60 ];
  G4double myGasArgonRout[ 60 ];

  G4double myLiqArgonZ[ 60 ];
  G4double myLiqArgonRout[ 60 ];

  G4double myRmin[60] = {0.};
  
  

  G4int    myNumPointsOuterCryo = 0;
  G4int    myNumPointsVacCryo   = 0;
  G4int    myNumPointsInnerCryo = 0;
  G4int    myNumPointsGasArgon  = 0;
  G4int    myNumPointsLiqArgon  = 0;

  int index = 0;


  // Outer Cryostat - Outer Surface
  DSIO::Get()->GetStreamDSG3CryostatProfile() >> myNumPointsOuterCryo;
  for( index = 0; index < myNumPointsOuterCryo; index++ ) { 
    DSIO::Get()->GetStreamDSG3CryostatProfile() >> myOuterCryostatZ[index] >> myOuterCryostatRout[index];
    myOuterCryostatRout[index] *= cm;
    myOuterCryostatZ[index] = myOuterCryostatZ[index]*cm + myCryostatShiftZ;
  }

  // Outer Cryostat - Inner Surface
  DSIO::Get()->GetStreamDSG3CryostatProfile() >> myNumPointsVacCryo;
  for( index = 0; index < myNumPointsVacCryo;   index++ ) {  
    DSIO::Get()->GetStreamDSG3CryostatProfile() >> myVacuumCryostatZ[index] >> myVacuumCryostatRout[index];
    myVacuumCryostatRout[index] *= cm;
    myVacuumCryostatZ[index] = myVacuumCryostatZ[index]*cm + myCryostatShiftZ;
  }

  // Inner Cryostat - Outer Surface
  DSIO::Get()->GetStreamDSG3CryostatProfile() >> myNumPointsInnerCryo;
  for( index = 0; index < myNumPointsInnerCryo; index++ ) {
    DSIO::Get()->GetStreamDSG3CryostatProfile() >> myInnerCryostatZ[index]  >> myInnerCryostatRout[index];
    myInnerCryostatRout[index] *= cm;
    myInnerCryostatZ[index] = myInnerCryostatZ[index]*cm + myCryostatShiftZ;
  }

  // GAr
  DSIO::Get()->GetStreamDSG3CryostatProfile() >> myNumPointsGasArgon;
  for( index = 0; index < myNumPointsGasArgon;  index++ ) {   
    DSIO::Get()->GetStreamDSG3CryostatProfile() >> myGasArgonZ[index]  >> myGasArgonRout[index];
    myGasArgonRout[index] *= cm;
    myGasArgonZ[index] = myGasArgonZ[index]*cm + myCryostatShiftZ;
  }

  // LAr
  DSIO::Get()->GetStreamDSG3CryostatProfile() >> myNumPointsLiqArgon;
  for( index = 0; index < myNumPointsLiqArgon;  index++ ) {
    DSIO::Get()->GetStreamDSG3CryostatProfile() >> myLiqArgonZ[index]  >> myLiqArgonRout[index];
    myLiqArgonRout[index] *= cm;
    myLiqArgonZ[index] = myLiqArgonZ[index]*cm + myCryostatShiftZ;
  }
  
  DSIO::Get()->CloseStreamDSG3CryostatProfile();


  // Outer Cryostat
  fSolidDSOuterCryostat  = new G4Polycone( "OuterCryostat0_Solid", 0., myTwoPi, myNumPointsOuterCryo, myOuterCryostatZ, myRmin, myOuterCryostatRout  );
  fLogicDSOuterCryostat  = new G4LogicalVolume( fSolidDSOuterCryostat, DSMaterial::Get()->GetStainlessSteel(), "OuterCryostat_Logic" );
  fPhysicDSOuterCryostat = new G4PVPlacement( 0, myZeros, "OuterCryostat", fLogicDSOuterCryostat, fMotherVolume, false, 0, myCheckOverlap );
                                                                                                          
  // Vacuum Region
  fSolidVacuumCryostat  = new G4Polycone( "VacuumCryostat0_Solid", 0, myTwoPi, myNumPointsVacCryo,  myVacuumCryostatZ, myRmin, myVacuumCryostatRout);
  fLogicVacuumCryostat  = new G4LogicalVolume( fSolidVacuumCryostat, DSMaterial::Get()->GetVacuum(), "VacuumCryostat_Logic" );
  fPhysicVacuumCryostat = new G4PVPlacement( 0, myZeros, "VacuumCryostat", fLogicVacuumCryostat, fPhysicDSOuterCryostat, false, 0, myCheckOverlap );

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
  fLogicOuterLiqArgon  = new G4LogicalVolume( fSolidOuterLiqArgon, DSMaterial::Get()->GetLiquidArgon(), "OuterLiquidArgon_Logic" );
  fPhysicOuterLiqArgon = new G4PVPlacement( 0, myZeros, "OuterLiquidArgon", fLogicOuterLiqArgon, fPhysicInnerCryostat, false, 0, myCheckOverlap );





  // TPC is composed by two coaxial cylinders, the field rings and the (copper) and the field cage (PTFE).
  // The field cage contains the gas pocket, the active LAr and the grid (immersed in LAr).
  // The anode and cathode windows limit the active volume.
  // The two volumes, PMTAssemblyTop and Bottom, placed above and below the TPC, host the PMTs.


  // ------------------- //
  // ---  Dimensions --- //
  // ------------------- //
                                                                                // Difference in z (from technical drawings) between:
  G4double myZGapCryostatsMatingSurfaces = 7.54125*cm;                          //  - the outer cryostat and the inner cryostat mating surfaces


  G4double myZGapInMatSurfPMTAssemblyTop = 14.2*cm;                             //  - the inner cryostat mating surface and the top surface of the top pmt assembly
  G4double myZGapPMTAssemblyTopAnodeWin = 1*cm;                                 //  - the bottom surface of the top pmt assembly and the top surface of the  anode window  
  G4double myZGapCathodeWinPMTAssemblyBot = myZGapPMTAssemblyTopAnodeWin;       //  - the bottom surface of the cathode window and the top surface of the bottom pmt assembly
  G4double myZGapAnodeWinGrid  = 1.3*cm;                                        //  - the bottom surface of the anode window and the grid. Hanguo: the value in the drawing
                                                                                //     is wrong. The correct one is ~1.3 cm, the grid being a few cm below LAr surface and the 
                                                                                //     gas pocket being 1 cm thick          
  

  G4double myPMTAssembly_d = 313.6*cm;
  G4double myPMTAssembly_h = 14.00*cm;                              
  G4double myCapDisk_h = 7.62*cm;

  G4double myAnodeWindow_d = 313.6*cm;
  G4double myAnodeWindow_h = 1.27*cm;  

  G4double myFieldCage_d     = 300*cm;
  G4double myFieldCage_h     = 300*cm;
  G4double myFieldCageWall_t = 3.81*cm;

  G4double myTPB_t = 0.1*mm;                // Same as DS50
  G4double myGasPocket_t = 1.0*cm ;         // From Hanguo 

  G4double myFieldRings_ext_d = myAnodeWindow_d; 
  G4double myFieldRings_h = myFieldCage_h;

  G4double myGrid_h = 0.01*cm;

  G4double myCathodeWindow_d = 313.6*cm;
  G4double myCathodeWindow_h = 1.27*cm;

  G4double myTPCShiftZ = myZGapCryostatsMatingSurfaces
                       + myZGapInMatSurfPMTAssemblyTop 
                       + myPMTAssembly_h 
                       + myZGapPMTAssemblyTopAnodeWin 
                       + myAnodeWindow_h 
                       + (myZGapAnodeWinGrid + myFieldCage_h)/2.;

  // Volumes Positions

  G4ThreeVector myTPCCenterPos( 0., 0., -myTPCShiftZ + myCryostatShiftZ );

  G4ThreeVector myGasPocketPos( 0., 0., myFieldCage_h/2.- myTPB_t - myGasPocket_t/2.);
  G4ThreeVector myActiveLArPos( 0., 0., -myGasPocket_t/2. );

  G4ThreeVector myGridPos( 0., 0., myFieldCage_h/2. - myZGapAnodeWinGrid + myGasPocket_t/2. );

  G4ThreeVector myAnodeWindowPos = myTPCCenterPos + G4ThreeVector( 0., 0., myFieldCage_h/2. + myAnodeWindow_h/2.);
  G4ThreeVector myCathodeWindowPos = myTPCCenterPos - G4ThreeVector( 0., 0., myFieldCage_h/2. + myCathodeWindow_h/2.);

  G4ThreeVector myPMTAssemblyTopPos = myAnodeWindowPos + G4ThreeVector( 0., 0.,  myAnodeWindow_h/2. + myZGapPMTAssemblyTopAnodeWin + myPMTAssembly_h/2.);
  G4ThreeVector myPMTAssemblyBotPos = myCathodeWindowPos - G4ThreeVector( 0., 0., myCathodeWindow_h/2. + myZGapCathodeWinPMTAssemblyBot + myPMTAssembly_h/2.);

  G4ThreeVector myCapDiskPos( 0., 0.,  myPMTAssembly_h/2. - myCapDisk_h/2.);

  G4RotationMatrix* rotX180 = new G4RotationMatrix;
  rotX180->rotateX( M_PI*rad );



  // ----------------------------------------------------- //
  // ---  Definition of the volumes from top to bottom --- //
  // ----------------------------------------------------- //

  // Top PMTs hosting volume
  fSolidPMTAssembly     = new G4Tubs( "PMTAssembly_Solid", 0, myPMTAssembly_d/2., myPMTAssembly_h/2., 0, myTwoPi );
  fLogicPMTAssemblyTop  = new G4LogicalVolume( fSolidPMTAssembly, DSMaterial::Get()->GetLiquidArgon(), "PMTAssemblyTop_Logic" ); 
  fPhysicPMTAssemblyTop = new G4PVPlacement( 0, myPMTAssemblyTopPos, "PMTAssemblyTop", fLogicPMTAssemblyTop, fPhysicOuterLiqArgon, false, 0, myCheckOverlap );

  fSolidCapDisk  = new G4Tubs( "CapDisk_Solid", 0, myPMTAssembly_d/2., myCapDisk_h/2., 0, myTwoPi );
  fLogicCapDiskTop  = new G4LogicalVolume( fSolidCapDisk, DSMaterial::Get()->GetMetalCopper(), "CapDiskTop_Logic" );
  fPhysicCapDiskTop = new G4PVPlacement( 0, myCapDiskPos, "CapDiskTop", fLogicCapDiskTop, fPhysicPMTAssemblyTop, false, 0, myCheckOverlap );

  fSolidAnodeWindow  = new G4Tubs( "AnodeWidow_Solid", 0, myAnodeWindow_d/2., myAnodeWindow_h/2., 0, myTwoPi );
  fLogicAnodeWindow  = new G4LogicalVolume( fSolidAnodeWindow, DSMaterial::Get()->GetFusedSilica(), "AnodeWindow_Logic");
  fPhysicAnodeWindow = new G4PVPlacement( 0, myAnodeWindowPos, "AnodeWindow", fLogicAnodeWindow, fPhysicOuterLiqArgon, false, 0, myCheckOverlap ); 


  fSolidFieldRings  = new G4Tubs( "FieldRings_Solid", 0, myFieldRings_ext_d/2., myFieldRings_h/2., 0, myTwoPi );
  fLogicFieldRings  = new G4LogicalVolume( fSolidFieldRings, DSMaterial::Get()->GetMetalCopper(), "FieldRings_Logic" );
  fPhysicFieldRings = new G4PVPlacement( 0, myTPCCenterPos, "FieldRings", fLogicFieldRings, fPhysicOuterLiqArgon, false, 0, myCheckOverlap );


  // PTFE cylindrical container 
  fSolidFieldCage  = new G4Tubs( "FieldCage_Solid", 0, myFieldCage_d/2. + myFieldCageWall_t, myFieldCage_h/2., 0, myTwoPi );
  fLogicFieldCage  = new G4LogicalVolume( fSolidFieldCage, DSMaterial::Get()->GetTeflon(), "FieldCage_Logic" );
  fPhysicFieldCage = new G4PVPlacement( 0, myZeros, "FieldCage", fLogicFieldCage, fPhysicFieldRings, false, 0, myCheckOverlap );


  fSolidTPB  = new G4Tubs( "TPB_Solid", 0, myFieldCage_d/2., myFieldCage_h/2., 0, myTwoPi );
  fLogicTPB  = new G4LogicalVolume( fSolidTPB, DSMaterial::Get()->GetTPB(), "TPB_Logic" );
  fPhysicTPB = new G4PVPlacement( 0, myZeros, "TPB", fLogicTPB, fPhysicFieldCage, false, 0, myCheckOverlap );


  fSolidGasPocket  = new G4Tubs( "GasPocket_Solid", 0, myFieldCage_d/2. - myTPB_t, myGasPocket_t/2., 0, myTwoPi );
  if ( DSParameters::Get()->GetWithGasPocket() ) 
    fLogicGasPocket  = new G4LogicalVolume( fSolidGasPocket, DSMaterial::Get()->GetGaseousArgon(), "GasPocket_Logic" );
  else
    fLogicGasPocket  = new G4LogicalVolume( fSolidGasPocket, DSMaterial::Get()->GetLiquidArgon(), "GasPocket_Logic" );
  fPhysicGasPocket = new G4PVPlacement( 0, myGasPocketPos, "GasPocket", fLogicGasPocket, fPhysicTPB, false, 0, myCheckOverlap );


  // Set the z coordinate of the LAr - GAr interface, necessary for S2 generation in DSLightX
  G4double myLArGArBoundaryPosZ = (myGasPocketPos + myTPCCenterPos).z() - myGasPocket_t/2.;
  DSStorage::Get()->SetLArGArBoundaryPosZ( myLArGArBoundaryPosZ );


  fSolidActiveLAr  = new G4Tubs( "ActiveLAr_Solid", 0, myFieldCage_d/2. - myTPB_t, (myFieldCage_h - myGasPocket_t)/2. - myTPB_t, 0, myTwoPi );
  fLogicActiveLAr  = new G4LogicalVolume( fSolidActiveLAr, DSMaterial::Get()->GetLiquidArgon(), "ActiveLAr_Logic" );
  fPhysicActiveLAr = new G4PVPlacement( 0, myActiveLArPos, "ActiveLAr", fLogicActiveLAr, fPhysicTPB, false, 0, myCheckOverlap );


  fSolidGrid  = new G4Tubs( "Grid_Solid", 0, myFieldCage_d/2. - myTPB_t, myGrid_h/2., 0, myTwoPi );
  fLogicGrid  = new G4LogicalVolume( fSolidGrid, DSMaterial::Get()->GetGridSteel(), "Grid_Logic" ); 
  fPhysicGrid = new G4PVPlacement( 0, myGridPos, "Grid", fLogicGrid, fPhysicActiveLAr, false, 0, myCheckOverlap );


  fSolidCathodeWindow  = new G4Tubs( "CathodeWindow_Solid", 0, myCathodeWindow_d/2., myCathodeWindow_h/2., 0, myTwoPi );
  fLogicCathodeWindow  = new G4LogicalVolume( fSolidCathodeWindow, DSMaterial::Get()->GetFusedSilica(), "CathodeWindow_Logic" );
  fPhysicCathodeWindow = new G4PVPlacement( 0, myCathodeWindowPos, "CathodeWindow", fLogicCathodeWindow, fPhysicOuterLiqArgon, false, 0, myCheckOverlap );


  // Bottom PMTs hosting volume
  fLogicPMTAssemblyBot  = new G4LogicalVolume( fSolidPMTAssembly, DSMaterial::Get()->GetLiquidArgon(), "PMTAssemblyBottom_Logic" ); 
  fPhysicPMTAssemblyBot = new G4PVPlacement( rotX180, myPMTAssemblyBotPos, "PMTAssemblyBottom", fLogicPMTAssemblyBot, fPhysicOuterLiqArgon, false, 0, myCheckOverlap );

  fLogicCapDiskBot  = new G4LogicalVolume( fSolidCapDisk, DSMaterial::Get()->GetMetalCopper(), "CapDiskBot_Logic" );
  fPhysicCapDiskBot = new G4PVPlacement( 0, myCapDiskPos, "CapDiskBot", fLogicCapDiskBot, fPhysicPMTAssemblyBot, false, 0, myCheckOverlap );

  // PMTs?
  //DSDetectorPMTDSG3( (G4VPhysicalVolume*) fPhysicPMTAssemblyTop, (G4VPhysicalVolume*) fPhysicCapDiskTop ); 
  DSDetectorPMTDSG3( fPhysicPMTAssemblyTop, fPhysicCapDiskTop ); 
  DSDetectorPMTDSG3( fPhysicPMTAssemblyBot, fPhysicCapDiskBot ); 


  // The definition of the active LAr region is needed to set the range cuts for this volume 
  // to a smaller value with respect to the rest of the detector ( see DSPhysicsList::SetCuts() )
  G4Region* fLArRegion = new G4Region("LAr_Logic");
  fLogicActiveLAr->SetRegion( fLArRegion );
  fLArRegion->AddRootLogicalVolume(fLogicActiveLAr);


  DefineSurfaces(); 

}



DSDetectorDSG3::~DSDetectorDSG3(){
  ; //delete fMessenger;
}



void  DSDetectorDSG3::DefineSurfaces() {

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
  fOpTPBGArSurface->SetType( dielectric_dielectric );
  fOpTPBGArSurface->SetModel( unified );
  fOpTPBGArSurface->SetFinish( ground );
  fOpTPBGArSurface->SetSigmaAlpha(0.3);
 
  const G4int NUM = 4;
  G4double pp[NUM] = {0.1*eV, 8.0*eV, 8.3*eV, 20.0*eV};  // vis, vis, UV, UV
  G4double specularlobe[NUM] = {0., 0., 0., 0.};         //--
  G4double specularspike[NUM] = {0., 0., 0., 0.};        //----  gives all reflection to Lambertian lobe
  G4double backscatter[NUM] = {0., 0., 0., 0.};          //-- 
  G4double reflectivity[NUM] = {1.0, 1.0, 1.0, 1.0};     //  To set 1-absorption
  G4double transmitivity[NUM] = {0.6, 0.6, 1.0, 1.0};    //  To set reflection vs. transmission, overridding Fresnel
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
  new G4LogicalBorderSurface("TopWindowTPBSurface",    fPhysicAnodeWindow,  fPhysicTPB, fOpWindowTPBSurface );
  new G4LogicalBorderSurface("BottomWindowTPBSurface", fPhysicCathodeWindow, fPhysicTPB, fOpWindowTPBSurface );
  new G4LogicalBorderSurface("TPBTopWindowSurface",    fPhysicTPB, fPhysicAnodeWindow,  fOpWindowTPBSurface );
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
  new G4LogicalBorderSurface("TopWindowLArSurface",    fPhysicAnodeWindow,   fPhysicOuterLiqArgon, fOpWindowLArSurface); 
  new G4LogicalBorderSurface("BottomWindowLArSurface", fPhysicCathodeWindow, fPhysicOuterLiqArgon, fOpWindowLArSurface); 
  new G4LogicalBorderSurface("LArTopWindowSurface",    fPhysicOuterLiqArgon, fPhysicAnodeWindow,   fOpWindowLArSurface); 
  new G4LogicalBorderSurface("LArBottomWindowSurface", fPhysicOuterLiqArgon, fPhysicCathodeWindow, fOpWindowLArSurface); 
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
  new G4LogicalBorderSurface("TPBTeflonSurface", fPhysicTPB, fPhysicFieldCage, fOpTPBTeflonSurface );
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
  new G4LogicalBorderSurface("InnerLArCapDiskTopSurface", fPhysicOuterLiqArgon, fPhysicCapDiskBot , fOpLArCapDiskSurface );
  new G4LogicalBorderSurface("InnerLArCapDiskBottomSurface", fPhysicOuterLiqArgon, fPhysicCapDiskTop , fOpLArCapDiskSurface );
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
  new G4LogicalBorderSurface("LArTeflonSurface", fPhysicFieldCage, fPhysicOuterLiqArgon, fOpLArTeflonSurface );
  new G4LogicalBorderSurface("LArPMTTopAssemblySurface", fPhysicPMTAssemblyTop, fPhysicOuterLiqArgon, fOpLArTeflonSurface);
  new G4LogicalBorderSurface("LArPMTBottomAssemblySurface", fPhysicPMTAssemblyBot, fPhysicOuterLiqArgon, fOpLArTeflonSurface);
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
  new G4LogicalBorderSurface("DSG3OuterSurface",
			     fMotherVolume,
			     fPhysicDSOuterCryostat,
			     fOpElectropolishedStainlessSteelSurface);



}  
