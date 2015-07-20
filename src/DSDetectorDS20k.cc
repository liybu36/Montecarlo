#include <iostream>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Orb.hh"

#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"

#include "DSDetectorDS20k.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSStorage.hh"
#include "DSParameters.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

using namespace std;

////////////////////////////////////////////////
//////        Detector Description      ////////
////////////////////////////////////////////////

// 24th April 2015
// schematics of the 20k detector, starting from scaling the 5K geometry (the cryostats in particular) and from some drawings from Hanguo
// OPTICS is directly derived from DS50
// 
// Placements are realized with no shifts. The correct position is calculated when making the solids.
//
// GENERAL SCHEME:
// - nested single volumes, from the outside to the center: Vacuum Jacket (steel), Vacuum, Inner Cryo (ssteel). 
// - the fill of the inner Volume is divided in two region at z = 1200, which correspond to the LAr-GAr interface (Inactive LAr and InactiveGar) 
//   - the TPC lays inside the InactiveLAr region (Copper Rings, Reflector, TPBSide , Active LAr)
//   - the Gas Pocket lays in the upper part (reflector, Active GAr)
// - other volumes: 
//   - the grid is placed 0.5 cm below the LAr surface
//   - Fused silica windows with SiPm on the outfacing surface and TPB on the innerfacing one
//   - Pseudo Argon Laryer (condensation on the top window)

// addition note on TPBSide and TPBTop
// TPBSide is a volume sorrunding the LAr Active volume on the side and bottom. TPBTop is only between
// LArLAyer (Pseudo Ar) and the top window

DSDetectorDS20k::DSDetectorDS20k(G4VPhysicalVolume *myMotherVolume) {

  G4Colour  myWhite   (1.0, 1.0, 1.0) ;  // white
  G4Colour  myGray    (0.5, 0.5, 0.5) ;  // gray
  G4Colour  myBlack   (0.0, 0.0, 0.0) ;  // black
  G4Colour  myRed     (1.0, 0.0, 0.0) ;  // red
  G4Colour  myGreen   (0.0, 1.0, 0.0) ;  // green
  G4Colour  myBlue    (0.0, 0.0, 1.0) ;  // blue
  G4Colour  myCyan    (0.0, 1.0, 1.0) ;  // cyan
  G4Colour  myMagenta (1.0, 0.0, 1.0) ;  // magenta 
  G4Colour  myYellow  (1.0, 1.0, 0.0) ;  // yellow
  
  fMotherVolume = myMotherVolume;
  DSLog(routine) << " Constructing DS20k Geometry" << endlog ;

  const double myTwoPi = 2*M_PI*rad;
  G4bool   myCheckOverlap   = DSStorage::Get()->GetCheckOverlap();
  G4ThreeVector myZeros( 0., 0., 0.);

  // Parameters from drawings  
  double TPCEdge             = 120*cm ;     //cold
  double myTeflonThickness   = 1.5*2.54*cm;   // 1.5 inches
  double myRingThickness     = 0.5*cm;
  double myTPCOuterRadius    = GetOctagonInnerRadius(TPCEdge);    // CHANGED from dav to pao
  double myTeflonOuterRadius = myTPCOuterRadius + myTeflonThickness;  
  double myRingOuterRadius   = myTeflonOuterRadius + myRingThickness ;
  double myWindowsThickness    = 2.54/4 * cm ;
  double myTPCHeight          = 2400*mm;
  // place the center of the TPC in the center of the ref system
  double myLArGArBoundaryPosZ = 1200*mm;
  DSStorage::Get()->SetLArGArBoundaryPosZ( myLArGArBoundaryPosZ );
 
  // Other parameters  
  double myTPBThickness        = 0.1 *mm ; 
  double myPArThickness        = 0.1 * mm ; 
  double myGasPocketThickness  = 1*cm ; 
  double mySiPmThickness       = 1. * mm ;
  double mySiPmBoardThickness  = 0.0*mm ; 
  double myActiveGasLogicThickness = myGasPocketThickness+myWindowsThickness+mySiPmThickness+mySiPmBoardThickness;
  double myTPCHeightLAr = myTPCHeight - myTPCHeight/2. + myLArGArBoundaryPosZ;
  double myTPCHeightGAr = myTPCHeight - myTPCHeight/2. - myLArGArBoundaryPosZ;
  
  
   // --------------------------- //
  // ------   Cryostats   ------ //
  // --------------------------- //
  
  //Load cryostat profiles                   
  // NOTE:  Each cryostat profile is evaluated in the ref frame with z = 0 at the CENTER of THE TPC
  G4double myCryostatShiftZ = 0.; 
  G4double LArGarIntefaceZ = 1200 ; 
  
  G4double myOuterCryostatZ[ 300 ];
  G4double myOuterCryostatRout[ 300 ];

  G4double myVacuumCryostatZ[ 300 ];
  G4double myVacuumCryostatRout[ 300 ];

  G4double myInnerCryostatZ[ 300 ];
  G4double myInnerCryostatRout[ 300 ];

  G4double myGasArgonZ[ 300 ];
  G4double myGasArgonRout[ 300 ];

  G4double myLiqArgonZ[ 300 ];
  G4double myLiqArgonRout[ 300 ];

  G4double myRmin[300];
  for(int ii = 0; ii<300;++ii) myRmin[ii] = 0. ;
  
  G4int    myNumPointsOuterCryo = 0;
  G4int    myNumPointsVacCryo   = 0;
  G4int    myNumPointsInnerCryo = 0;
  G4int    myNumPointsLarGAs_tot = 0;
  G4int    myNumPointsGasArgon  = 0;
  G4int    myNumPointsLiqArgon  = 0;
  int i =0 ; 
  // Outer Cryostat - Outer Surface
  ifstream profileOuterVacCryo("../data/detector/DS20kVacuumOuterProfile.dat") ; 
  profileOuterVacCryo >> myNumPointsOuterCryo;
  for( i = 0; i < myNumPointsOuterCryo; i++ ) { 
    profileOuterVacCryo >> myOuterCryostatZ[i]  >> myOuterCryostatRout[i];
    myOuterCryostatZ[i] += myCryostatShiftZ;
  }
  profileOuterVacCryo.close(); 

  // Outer Cryostat - Inner Surface 
  ifstream profileInnerVacCryo("../data/detector/DS20kVacuumInnerProfile.dat") ; 
  profileInnerVacCryo  >> myNumPointsVacCryo;
  for( i = 0; i < myNumPointsVacCryo;   i++ ) {  
    profileInnerVacCryo >> myVacuumCryostatZ[i] >> myVacuumCryostatRout[i];
    myVacuumCryostatZ[i] += myCryostatShiftZ;
  }
  profileInnerVacCryo.close(); 

  // Inner Cryostat - Outer Surface 
  ifstream profileOuterCryo("../data/detector/DS20kCryoOuterProfile.dat") ; 
  profileOuterCryo  >>myNumPointsInnerCryo ;
  for( i = 0; i < myNumPointsInnerCryo;   i++ ) {  
    profileOuterCryo >> myInnerCryostatZ[i] >> myInnerCryostatRout[i];
    myInnerCryostatZ[i] += myCryostatShiftZ;
  }
  profileInnerVacCryo.close(); 
  // Inactive LAr and Gar - Outer Surface   - to be divided in 2 parts
  ifstream profileInnerCryo("../data/detector/DS20kCryoInnerProfile.dat") ; 
  profileInnerCryo >> myNumPointsLarGAs_tot;
  for( i = 0; i < myNumPointsLarGAs_tot; i++ ) {
    G4double tmpreadz, tmpreadr ; 
    profileInnerCryo>> tmpreadz  >> tmpreadr ; 
    if (tmpreadz>LArGarIntefaceZ) {
       myGasArgonRout[i-myNumPointsLiqArgon+1] =tmpreadr ;
       myGasArgonZ[i-myNumPointsLiqArgon+1]= tmpreadz ; 
       ++myNumPointsGasArgon; 
    } else {
       myLiqArgonRout[i]=tmpreadr ;
       myLiqArgonZ[i]= tmpreadz ; 
       ++myNumPointsLiqArgon; 
    }
    myLiqArgonZ[i] += myCryostatShiftZ;
    myGasArgonZ[i] += myCryostatShiftZ;
  }
  myGasArgonZ[0] =LArGarIntefaceZ  ;
  myGasArgonRout[0] = myLiqArgonRout[myNumPointsLiqArgon-1];
  myLiqArgonZ[myNumPointsLiqArgon]    =LArGarIntefaceZ ;
  myLiqArgonRout[myNumPointsLiqArgon] = myLiqArgonRout[myNumPointsLiqArgon-1];
  profileInnerCryo.close(); 
  
  //-------------------------//
  //      Outer VacCryo      //
  //-------------------------//
    
  G4Polycone      *fSolidDS20k   = new G4Polycone( "OuterCryostat_Solid", 0, myTwoPi, myNumPointsOuterCryo,myOuterCryostatZ, myRmin, myOuterCryostatRout  );
  G4LogicalVolume *fLogicDS20k   = new G4LogicalVolume( fSolidDS20k, DSMaterial::Get()->GetStainlessSteel(), "OuterCryostat_Logic" );
  G4PVPlacement   *fPhysicDS20k  = new G4PVPlacement( 0, G4ThreeVector(0.,0.,0.), "OuterCryostat", fLogicDS20k,fMotherVolume, false, 0, myCheckOverlap ); 
  //fLogicDS20k->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  G4Polycone      *fSolidDS20kVac   = new G4Polycone( "OuterCryostat_SolidVac", 0, myTwoPi, myNumPointsVacCryo,myVacuumCryostatZ, myRmin, myVacuumCryostatRout  );
  G4LogicalVolume *fLogicDS20kVac   = new G4LogicalVolume( fSolidDS20kVac, DSMaterial::Get()->GetVacuum(), "OuterCryostat_LogicVac" );
  G4PVPlacement   *fPhysicDS20kVac  = new G4PVPlacement( 0, myZeros, "OuterCryostatVac", fLogicDS20kVac, fPhysicDS20k, false, 0, myCheckOverlap );  
  G4VisAttributes *myVacuumAttributes = new G4VisAttributes(myWhite);

  //-------------------------//
  //      Inner Cryo         //
  //-------------------------//
  
  G4Polycone      *fSolidInnerCryo   = new G4Polycone( "InnerCryo_Solid", 0, myTwoPi, myNumPointsInnerCryo,myInnerCryostatZ, myRmin,  myInnerCryostatRout );
  G4LogicalVolume *fLogicInnerCryo   = new G4LogicalVolume( fSolidInnerCryo, DSMaterial::Get()->GetStainlessSteel(), "SolidInnerCryo_Logic" );
  G4PVPlacement   *fPhysicInnerCryo  = new G4PVPlacement( 0, myZeros, "InnerCryo", fLogicInnerCryo,fPhysicDS20kVac , false, 0, myCheckOverlap );  
fLogicInnerCryo->SetVisAttributes(G4Colour  (1.0, 1.0, 1.0, 0.01));
  //  fLogicInnerCryo->SetVisAttributes(G4VisAttributes::GetInvisible());
 
  // evrything must be divided in 2 pieces in Z for LAr and GAr regions  
  G4Polycone      *fSolidInactiveLar   = new G4Polycone( "SolidInactiveLar_Solid", 0, myTwoPi, myNumPointsLiqArgon+1,myLiqArgonZ, myRmin, myLiqArgonRout  );
  G4LogicalVolume *fLogicInactiveLar   = new G4LogicalVolume( fSolidInactiveLar, DSMaterial::Get()->GetNSLiquidArgon(), "SolidInactiveLar_Logic" );
  G4PVPlacement   *fPhysicInactiveLar  = new G4PVPlacement( 0, myZeros, "SolidInactiveLar", fLogicInactiveLar, fPhysicInnerCryo, false, 0, myCheckOverlap );  
  G4VisAttributes *myLarAttributes = new G4VisAttributes(G4Colour  (0.0, .1, .9, 0.9));
  fLogicInactiveLar->SetVisAttributes(G4VisAttributes::GetInvisible());
//  fLogicInnerCryo->SetVisAttributes(G4Colour  (0.0, 0.0, 1.0));
  
  G4Polycone      *fSolidInactiveGar   = new G4Polycone( "SolidInactiveGar_Solid", 0, myTwoPi,myNumPointsGasArgon ,myGasArgonZ, myRmin, myGasArgonRout  );
  G4LogicalVolume *fLogicInactiveGar   = new G4LogicalVolume( fSolidInactiveGar, DSMaterial::Get()->GetGaseousArgon(), "SolidInactiveGar_Logic" );
  G4PVPlacement   *fPhysicInactiveGar  = new G4PVPlacement( 0, myZeros, "SolidInactiveGar", fLogicInactiveGar,fPhysicInnerCryo , false, 0, myCheckOverlap );  
  G4VisAttributes *myGarAttributes = new G4VisAttributes(G4Colour  (0.0, 0.0, 1.0, 0.9));
  fLogicInactiveGar->SetVisAttributes(myGarAttributes);

  //-------------------------//
  //    Shaping Rings        //
  //-------------------------//
  // field shaping rings for now assume a constant wall of 5mm thick outside the PTFE in 
  // average for MC first order estimation should be enough for now. 
  //(6mm thick 2.2cm high ring spaced 2.5cm and 95 rings total).

  
  double myZRingBottom[2] 	 = {+LArGarIntefaceZ-myTPCHeight,LArGarIntefaceZ};
  
  double myZRingTop[2]		 = {-myTPCHeightGAr/2., +myTPCHeightGAr/2.};
  double myInnerRadiusRing[2] = {0.,0.};
  double myOuterRadiusRing[2] = {myRingOuterRadius,myRingOuterRadius};
  
//ring top (necessary?) -  up to now, the tpc is totally immersed in Inactive Liquid Argon, no upper part (in GAr) is required
/*  G4Polyhedra *fSolidRingTop    = new G4Polyhedra( "RingTop_Solid",
        					    0,
        					    myTwoPi,
        					    8,
        					    2,
        					    myZRingTop,
        					    myInnerRadiusRing,
        					    myOuterRadiusRing );

  G4LogicalVolume *fLogicRingTop    = new G4LogicalVolume(fSolidRingTop,
				                         DSMaterial::Get()->GetMetalCopper(),
				                         "RingTop_Logic");
  
  G4VPhysicalVolume *fPhysicRingTop     = new G4PVPlacement(0,
							   G4ThreeVector(0,0,0 ),
							   //G4ThreeVector(0,0, 0.),
							   "FieldCageRingTop",
							   fLogicRingTop,
							   fPhysicInactiveGar, //fMotherVolume,
							   false,
							   0,
							   myCheckOverlap);
*/  //ring bottom
  G4Polyhedra *fSolidRingBottom = new G4Polyhedra( "RingBottom_Solid",
        					    0,
        					    myTwoPi,
        					    8,
        					    2,
        					    myZRingBottom,
        					    myInnerRadiusRing,
        					    myOuterRadiusRing );
  
  G4LogicalVolume *fLogicRingBottom = new G4LogicalVolume(fSolidRingBottom,
				                         DSMaterial::Get()->GetMetalCopper(),
				                         "RingBottom_Logic");

  G4VPhysicalVolume *fPhysicRingBottom  = new G4PVPlacement(0,
							   G4ThreeVector(0,0, 0),
							   "FieldCageRingBottom",
							   fLogicRingBottom,
							   fPhysicInactiveLar, //fMotherVolume,
							   false,
							   0,
							   myCheckOverlap);
    
  
  G4VisAttributes *myRingBottomAttributes = new G4VisAttributes(myRed);
  fLogicRingBottom->SetVisAttributes(myRingBottomAttributes);
  
  //G4VisAttributes *myRingTopAttributes = new G4VisAttributes(myBlue);
  //fLogicRingTop->SetVisAttributes(myRingTopAttributes);

  //-------------------------//
  //    Teflon Reflector       //
  //-------------------------//
  
  // outwards from the inner active wall, a 1-inch thick PTFE reflector (octagonal shape).
  double myZTeflonBottom[2]          = { LArGarIntefaceZ-myTPCHeight,LArGarIntefaceZ};
  double myZTeflonTop[2]             = {LArGarIntefaceZ,LArGarIntefaceZ + myGasPocketThickness};
  double myInnerRadiusTeflon[2] = {0,0};
  double myOuterRadiusTeflon[2] = {myTeflonOuterRadius,myTeflonOuterRadius};
  
  //reflector bottom 
  G4Polyhedra *fSolidTeflonBottom = new G4Polyhedra( "TeflonBottom_Solid",
        					    0,
        					    myTwoPi,
        					    8,
        					    2,
        					    myZTeflonBottom,
        					    myInnerRadiusTeflon,
        					    myOuterRadiusTeflon );
  
  
  G4LogicalVolume *fLogicTeflonBottom = new G4LogicalVolume(fSolidTeflonBottom,
				                         DSMaterial::Get()->GetTeflon(),
				                         "TeflonBottom_Logic");


  fPhysicTeflonBottom  = new G4PVPlacement(0,
							   G4ThreeVector(0,0,0),
							   "TeflonBottom",
							   fLogicTeflonBottom,
							   fPhysicRingBottom,
							   false,
							   0,
							   myCheckOverlap);
  //-------------------------//
  //          TPC            //
  //-------------------------//
  
  // TPC active volume section: Octagonal shape with edge 64.2cm when cold. (or Inscribed circle radius = 77.5cm)
  double myZTPCBottom[2]          = { LArGarIntefaceZ-myTPCHeight+myTPBThickness,LArGarIntefaceZ};
  double myZTPBSide[2]          = { LArGarIntefaceZ-myTPCHeight,LArGarIntefaceZ};
  double myInnerRadiusTPC[2] = {0,0};
  double myOuterRadiusTPC[2] = {myTPCOuterRadius,myTPCOuterRadius};
  double myOuterRadiusTPBSide[2] = {myTPCOuterRadius+myTPBThickness,myTPCOuterRadius+myTPBThickness};


  
  // TPB layer (side and bottom) 
  G4Polyhedra *fSolidTPBSide = new G4Polyhedra( "TPBSide_Solid",
        					    0,
        					    myTwoPi,
        					    8,
        					    2,
        					    myZTPBSide,
        					    myInnerRadiusTPC,
        					    myOuterRadiusTPBSide );
  
  
  
  G4LogicalVolume *fLogicTPBSide= new G4LogicalVolume(fSolidTPBSide,
				                         DSMaterial::Get()->GetTPB(),
				                         "TPBSide_Logic");
  fPhysicTPBSide  = new G4PVPlacement(0,
	             				   G4ThreeVector(0,0,0),
							   "TPBSide",
							   fLogicTPBSide,
							   fPhysicTeflonBottom,
							   false,
							   0,
							   myCheckOverlap);
  // TPC Active Volume
  G4Polyhedra *fSolidActiveLAr = new G4Polyhedra( "ActiveLAr_Solid",
        					    0,
        					    myTwoPi,
        					    8,
        					    2,
        					    myZTPCBottom,
        					    myInnerRadiusTPC,
        					    myOuterRadiusTPC );
  
  
  
  G4LogicalVolume *fLogicActiveLAr= new G4LogicalVolume(fSolidActiveLAr,
				                         DSMaterial::Get()->GetLiquidArgon(),
				                         "LAr_Logic");
  fPhysicActiveLAr  = new G4PVPlacement(0,
	             				   G4ThreeVector(0,0,0),
							   "ActiveLAr",
							   fLogicActiveLAr,
							   fPhysicTPBSide, //fPhysicTPBSide,//fPhysicRingBottom,
							   false,
							   0,
							   myCheckOverlap);
  
  fLogicActiveLAr->SetVisAttributes(myBlue);
  G4Region* fLArRegion = new G4Region("LAr_Logic");
  fLogicActiveLAr->SetRegion( fLArRegion );
  fLArRegion->AddRootLogicalVolume(fLogicActiveLAr);



  //-------------------------------------------//
  //  TopReflector and  Active GAr Region      //
  //-------------------------------------------//
    
  //reflector top 
  G4Polyhedra *fSolidTeflonTop    = new G4Polyhedra( "TeflonTop_Solid",
        					    0,
        					    myTwoPi,
        					    8,
        					    2,
        					    myZTeflonTop,
        					    myOuterRadiusTPC,
        					    myOuterRadiusRing );
  G4LogicalVolume *fLogicTeflonTop    = new G4LogicalVolume(fSolidTeflonTop,
				                         DSMaterial::Get()->GetTeflon(),
				                         "TeflonTop_Logic");

  G4VPhysicalVolume *fPhysicTeflonTop     = new G4PVPlacement(0,
							   G4ThreeVector(0,0,0),
							   "TeflonTop",
							   fLogicTeflonBottom,
							   fPhysicInactiveGar , 
							   false,
							   0,
							   myCheckOverlap);
  double myZTPCTop[2]          = { LArGarIntefaceZ,LArGarIntefaceZ+myActiveGasLogicThickness};
  
  G4Polyhedra *fSolidGasPocket    = new G4Polyhedra( "GasPocket_Solid",
        					    0,
        					    myTwoPi,
        					    8,
        					    2,
        					    myZTPCTop,
        					    myInnerRadiusTPC,
        					    myOuterRadiusTPC );

  G4LogicalVolume *fLogicGasPocket    = new G4LogicalVolume(fSolidGasPocket,
				                         DSMaterial::Get()->GetGaseousArgon(),
				                         "GasPocket_Logic");

  fPhysicGasPocket     = new G4PVPlacement(0,
							   G4ThreeVector(0,0,0),
							   "GasPocket",
							   fLogicGasPocket,
							   fPhysicInactiveGar,//fPhysicRingTop,
							   false,
							   0,
							   myCheckOverlap);
 //-------------------------------//
  //          Top and Bottom caps  //
  //-------------------------------//
  
  // Top/Bottom of the active volume section: Octagonal shape with edge 120cm when cold. (or Inscribed circle radius = ~150cm)

  double myZTopWindow[2]          = {LArGarIntefaceZ+myGasPocketThickness,LArGarIntefaceZ+myGasPocketThickness+myWindowsThickness};
  double myZBotWindow[2]          = { LArGarIntefaceZ-myTPCHeight-myWindowsThickness, LArGarIntefaceZ-myTPCHeight};
  double myZTPBTop[2]          = {LArGarIntefaceZ+myGasPocketThickness-myTPBThickness,LArGarIntefaceZ+myGasPocketThickness};
  double myZPArLayerTop[2]          = {LArGarIntefaceZ+myGasPocketThickness -myTPBThickness -myPArThickness,LArGarIntefaceZ+myGasPocketThickness-myTPBThickness};
  double myZTPBBottom[2]          = { LArGarIntefaceZ-myTPCHeight, LArGarIntefaceZ-myTPCHeight+myTPBThickness};
  double z1sipm = myZTopWindow[1]; 
  double z2sipm = myZTopWindow[1]+mySiPmBoardThickness+mySiPmThickness ; 
  double myZSiPmTop[2]             = { z1sipm, z2sipm};
  z1sipm =  myZBotWindow[0]-mySiPmBoardThickness-mySiPmThickness;
  z2sipm =  myZBotWindow[0];
  double myZSiPmBottom[2]          = { z1sipm, z2sipm};
  
  
  
  // top window
  G4Polyhedra *fSolidTopWindow = new G4Polyhedra( "TopWindow_Solid",  0,
        					    myTwoPi,   8,
        					    2,  myZTopWindow,
        					    myInnerRadiusTPC,
        					    myOuterRadiusTPC );
 G4LogicalVolume *fLogicTopWindow   = new G4LogicalVolume(fSolidTopWindow,
				                         DSMaterial::Get()->GetFusedSilica(),
				                         "TopWindow_Logic");

  fPhysicTopWindow  = new G4PVPlacement(0,
							   G4ThreeVector(0,0, 0),
							   "TopWindow",  fLogicTopWindow,
							   fPhysicGasPocket,   false,
							   0,   myCheckOverlap);
  // TPB in between gar and window (top) 
  G4Polyhedra * fSolidTPBTop = new G4Polyhedra( "TPBTop_Solid",  0,
        					    myTwoPi,   8,
        					    2,  myZTPBTop,
        					    myInnerRadiusTPC,
        					    myOuterRadiusTPC );
  G4LogicalVolume *fLogicTPBTop   = new G4LogicalVolume(fSolidTPBTop,
				                         DSMaterial::Get()->GetTPB(),
				                         "TPBTop_Logic");

  fPhysicTPBTop  = new G4PVPlacement(0,
					  G4ThreeVector(0,0,0),
					  "TPBTop",  fLogicTPBTop,
					  fPhysicGasPocket,   false,
					  0,   myCheckOverlap);
  // Pseudo Ar layer 
  G4Polyhedra * fSolidLArLayer= new G4Polyhedra( "LArLayer_Solid",  0,
        					    myTwoPi,   8,
        					    2, myZPArLayerTop ,
        					    myInnerRadiusTPC,
        					    myOuterRadiusTPC );
  G4LogicalVolume *fLogicLArLayer   = new G4LogicalVolume(fSolidLArLayer,
				                         DSMaterial::Get()->GetPseudoArgon(),
				                         "LArLayer_Logic");

  fPhysicTPBTop  = new G4PVPlacement(0,
					  G4ThreeVector(0,0, 0.),
					  "LArLayer",  fLogicLArLayer,
					  fPhysicGasPocket,   false,
					  0,   myCheckOverlap);
  // bottom window
  G4Polyhedra * fSolidBotWindow   = new G4Polyhedra( "BotWindow_Solid",
        					    0,
        					    myTwoPi,
        					    8,
        					    2,
        					    myZBotWindow,
        					    myInnerRadiusTPC,
        					    myOuterRadiusTPC );
  
  
   G4LogicalVolume *fLogicBotWindow    = new G4LogicalVolume(fSolidBotWindow,
				                         DSMaterial::Get()->GetFusedSilica(),
				                         "BotWindow_Logic");

  fPhysicBotWinwow    = new G4PVPlacement(0,
							   G4ThreeVector(0,0,0),
							   "BotWindow",
							   fLogicBotWindow,
							   fPhysicInactiveLar,//fPhysicRingTop,
							   false,
							   0,
							   myCheckOverlap);
fLogicBotWindow->SetVisAttributes(G4Color(1.0, 1.0, 0.0, 0.6));
fLogicTopWindow->SetVisAttributes(G4Color(1.0, 1.0, 0.0, 0.6));


  //---------------------------//
  //      SiPM arrays          //
  //--------------------------//

  // Top Array 
  G4Polyhedra *fSolidSiPMTop = new G4Polyhedra( "SiPMTop_Solid",
        					    0,
        					    myTwoPi,
        					    8,
        					    2,
        					    myZSiPmTop,
        					    myInnerRadiusTPC,
        					    myOuterRadiusTPC );
  
  G4LogicalVolume *fLogicSiPMTop    = new G4LogicalVolume(fSolidSiPMTop,
				                         DSMaterial::Get()->GetMetalSilicon(),
				                         "SiPMTop_Logic");

  fPhysicSiPmTop  = new G4PVPlacement(0,
							   G4ThreeVector(0,0,0),
							   "SiPMTop",
							   fLogicSiPMTop,
							   fPhysicGasPocket,//fPhysicRingBottom,
							   false,
							   0,
							   myCheckOverlap);
  
  // Bottom Array 
  G4Polyhedra *fSolidSiPMBottom = new G4Polyhedra( "SiPMBottom_Solid",
        					    0,
        					    myTwoPi,
        					    8,
        					    2,
        					    myZSiPmBottom,
        					    myInnerRadiusTPC,
        					    myOuterRadiusTPC );
  
  G4LogicalVolume *fLogicSiPMBottom    = new G4LogicalVolume(fSolidSiPMBottom,
				                         DSMaterial::Get()->GetMetalSilicon(),
				                         "SiPMBottom_Logic");

  fPhysicSiPmBottom  = new G4PVPlacement(0,
							   G4ThreeVector(0,0,0),
							   "SiPMBottom",
							   fLogicSiPMBottom,
							   fPhysicInactiveLar,//fPhysicRingBottom,
							   false,
							   0,
							   myCheckOverlap);
  
  // Grid
  G4double myGrid_h = 0.01 * cm;
  double myZGrid[2] = {LArGarIntefaceZ-.5*cm,  LArGarIntefaceZ-.5*cm+myGrid_h};
  G4ThreeVector myGridPosition( 0, 0, myLArGArBoundaryPosZ - 0.5*cm );

  fSolidGrid  =  new G4Polyhedra( "Grid_Solid",
        					    0,
        					    myTwoPi,
        					    8,
        					    2,
        					    myZGrid,
        					    myInnerRadiusTPC,
        					    myOuterRadiusTPC );
						    
  fLogicGrid  = new G4LogicalVolume( fSolidGrid, DSMaterial::Get()->GetGridSteel(), "Grid_Logic" );
  fPhysicGrid = new G4PVPlacement( 0, myZeros, "Grid", fLogicGrid,fPhysicActiveLAr , false, 0, myCheckOverlap );

  DefineSurfaces(); 

  // make SiPM as pe storing material 
  DSStorage::Get()->SetPMTMaterialIndex(fLogicSiPMBottom->GetMaterial()->GetIndex());
  
  
}

DSDetectorDS20k::~DSDetectorDS20k(){
  ; //delete fMessenger;
}

void  DSDetectorDS20k::DefineSurfaces() {
  ////////////////////////////////////////
  // Copied from DS DetectorDS50 - 24th april 2015
  ////////////////////////////////////////
    
  // LAR-Grid  
  G4OpticalSurface *fOpGridLArSurface = new G4OpticalSurface("OpGridLArSurface");
  new G4LogicalBorderSurface("GridLArSurface", fPhysicActiveLAr, fPhysicGrid, fOpGridLArSurface); 
  fOpGridLArSurface->SetType( dielectric_dielectric );
  fOpGridLArSurface->SetModel( glisur );
  fOpGridLArSurface->SetFinish( polished );
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
  fOpGridLArSurface->SetMaterialPropertiesTable( fGridLArSurfProp );

  // Grid->LAR 
  if (DSParameters::Get()->GetWithNewGridModel()) {
    G4OpticalSurface *fOpLArGridSurface = new G4OpticalSurface("OpLArGridSurface");
    new G4LogicalBorderSurface("LArGridSurface", fPhysicGrid, fPhysicActiveLAr, fOpLArGridSurface);
  cout << " With DS50 new grid model " << endl ; 
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
  G4OpticalSurface *fOpTPBGArSurface = new G4OpticalSurface("OpTPBGArSurface");
  new G4LogicalBorderSurface("GArTPBSurface", fPhysicGasPocket, fPhysicTPBTop, fOpTPBGArSurface );
  new G4LogicalBorderSurface("TPBGArSurface", fPhysicTPBTop, fPhysicGasPocket, fOpTPBGArSurface );
  new G4LogicalBorderSurface("LArTPBSurfaceSide", fPhysicActiveLAr, fPhysicTPBSide, fOpTPBGArSurface );
  new G4LogicalBorderSurface("TPBLArSurfaceSide", fPhysicTPBSide, fPhysicActiveLAr, fOpTPBGArSurface );
//  new G4LogicalBorderSurface("LArTPBSurfaceBot", fPhysicActiveLAr, fPhysicTPBBottom, fOpTPBGArSurface );
//  new G4LogicalBorderSurface("TPBLArSurfaceBot", fPhysicTPBBottom, fPhysicActiveLAr, fOpTPBGArSurface );
  new G4LogicalBorderSurface("LArLayerTPBSurface", fPhysicLArLayer, fPhysicTPBTop, fOpTPBGArSurface );
  new G4LogicalBorderSurface("TPBLArLayerSurface", fPhysicTPBTop, fPhysicLArLayer, fOpTPBGArSurface );
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
  new G4LogicalBorderSurface("TopWindowTPBSurface",   fPhysicTopWindow ,  fPhysicTPBTop, fOpWindowTPBSurface );
  new G4LogicalBorderSurface("BottomWindowTPBSurface", fPhysicBotWinwow, fPhysicTPBSide, fOpWindowTPBSurface );
  new G4LogicalBorderSurface("TPBTopWindowSurface",    fPhysicTPBTop,fPhysicTopWindow ,  fOpWindowTPBSurface );
  new G4LogicalBorderSurface("TPBBottomWindowSurface", fPhysicTPBSide, fPhysicBotWinwow, fOpWindowTPBSurface );
  fOpWindowTPBSurface->SetType( dielectric_dielectric );
  fOpWindowTPBSurface->SetModel( unified );
  //  fOpWindowTPBSurface->SetFinish( ground );
  fOpWindowTPBSurface->SetFinish( polished );
  //G4MaterialPropertiesTable *fWindowTPBSurfProp = new G4MaterialPropertiesTable();
  fOpWindowTPBSurface->SetMaterialPropertiesTable( fITOSurfProp );

  ////////////////////////////////////////
  // BellTop (fusedsilica) <--> SiPM and CathodeWindow <--> SiPM (both with ITO)
  // Make this bi-directional.
  ////////////////////////////////////////
  G4OpticalSurface *fOpWindowSiPmSurface     = new G4OpticalSurface("OpWindowSiPmSurface");
  new G4LogicalBorderSurface("TopWindowSiPmSurface",   fPhysicTopWindow ,  fPhysicSiPmTop, fOpWindowSiPmSurface );
  new G4LogicalBorderSurface("BottomWindowSiPmSurface", fPhysicBotWinwow, fPhysicSiPmBottom, fOpWindowSiPmSurface );
  //new G4LogicalBorderSurface("SiPmTopWindowSurface",    fPhysicSiPmTop,fPhysicTopWindow ,  fOpWindowSiPmSurface );
  //new G4LogicalBorderSurface("SiPmBottomWindowSurface", fPhysicSiPmBottom, fPhysicBotWinwow, fOpWindowSiPmSurface );
  fOpWindowSiPmSurface->SetType( dielectric_dielectric );
  fOpWindowSiPmSurface->SetModel( unified );
  //  fOpWindowSiPmSurface->SetFinish( ground );
  fOpWindowSiPmSurface->SetFinish( polished );
  //G4MaterialPropertiesTable *fWindowSiPmSurfProp = new G4MaterialPropertiesTable();
  fOpWindowSiPmSurface->SetMaterialPropertiesTable( fITOSurfProp );

  ///////////////
  // TPB --> Teflon (Reflector)
  //  Should be no Teflon --> TPB as surface is defined as dielectric_metal.
  ////////////////////////////////////////
  G4double TeflonTPBENE[4] = {0.1*eV, 8.0*eV, 8.3*eV, 20.0*eV};
  G4double TREFUV  = DSParameters::Get()->GetTeflonTPBUVRef();
  G4double TREFVIS = DSParameters::Get()->GetTeflonTPBVisRef();
  G4double TeflonTPBREF[4] = {TREFVIS, TREFVIS ,TREFUV , TREFUV };
  G4OpticalSurface *fOpTPBTeflonSurface = new G4OpticalSurface("OpTBPTeflonSurface");
  new G4LogicalBorderSurface("TPBTeflonSurface", fPhysicTPBSide,fPhysicTeflonBottom , fOpTPBTeflonSurface );
  fOpTPBTeflonSurface->SetType( dielectric_metal );
  fOpTPBTeflonSurface->SetModel(unified);

  // PDM: though I can't see how, the following settings are giving the desired Lambertian reflection 
  fOpTPBTeflonSurface->SetFinish(groundfrontpainted);                           
  //fOpTPBTeflonSurface->SetFinish(ground);                           
  fOpTPBTeflonSurface->SetSigmaAlpha(0.1);
                          
  G4MaterialPropertiesTable *fTPBTeflonSurfProp = new G4MaterialPropertiesTable();
  fTPBTeflonSurfProp->AddProperty("REFLECTIVITY", TeflonTPBENE, TeflonTPBREF, 4);			 
  fOpTPBTeflonSurface->SetMaterialPropertiesTable( fTPBTeflonSurfProp );
  
}  

double DSDetectorDS20k::GetOctagonInnerRadius(double edge) {  return edge /2 *(1+sqrt(2)); }
double DSDetectorDS20k::GetOctagonOuterRadius(double edge) {  return edge /2. *sqrt(4  + 2 *sqrt(2)) ; }
