#include <iostream>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Orb.hh"

#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"

#include "DSDetectorDS5k.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSStorage.hh"
#include "DSParameters.hh"
//#include "DSDetectorPMTDS5k.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

using namespace std;

////////////////////////////////////////////////
//////        Detector Description      ////////
////////////////////////////////////////////////


DSDetectorDS5k::DSDetectorDS5k(G4VPhysicalVolume *myMotherVolume) {

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
  DSLog(routine) << " Constructing DS5k Geometry" << endlog ;

  const double myTwoPi = 2*M_PI*rad;
  G4bool   myCheckOverlap   = DSStorage::Get()->GetCheckOverlap();


  G4ThreeVector myZeros( 0., 0., 0.);
 
  double myTPCHeight          = 1800*mm;
  double myLArGArBoundaryPosZ = 1599*mm;

  DSStorage::Get()->SetLArGArBoundaryPosZ( myLArGArBoundaryPosZ );
    
  double TPCEdge = 64.2*cm ;
    
  double myTeflonThickness = 2.54*cm;
  double myRingThickness = 0.5*cm;
  

  double myTPCOuterRadius    = GetOctagonInnerRadius(TPCEdge);    // CHANGED from dav to pao
  double myTeflonOuterRadius = myTPCOuterRadius + myTeflonThickness;  
  double myRingOuterRadius   = myTeflonOuterRadius + myRingThickness ;

  //CHANGED from dav to pao
  //double myTPCHeightLAr = myTPCHeight - myTPCHeight/2. + myLArGArBoundaryPosZ;
  //double myTPCHeightGAr = myTPCHeight - myTPCHeight/2. - myLArGArBoundaryPosZ;
  double myGasPocketThickness = 1*cm ; 
  double myWindowsThickness = 2.4/4 * cm ;
  double mySiPmThickness = 1. * mm ;
  double mySiPmBoardThickness  = 0.0*mm ; 
  double myActiveGasLogicThickness = myGasPocketThickness+myWindowsThickness+mySiPmThickness+mySiPmBoardThickness;
  double myTPCHeightLAr = myTPCHeight - myTPCHeight/2. + myLArGArBoundaryPosZ;
  double myTPCHeightGAr = myTPCHeight - myTPCHeight/2. - myLArGArBoundaryPosZ;
  
  
  // --------------------------- //
  // --- Scintillator Sphere --- //
  // --------------------------- //
  G4Orb* fScintillatorSolid = new G4Orb("Scintillator_Solid", 2.5*m );  
  G4LogicalVolume* fScintillatorLogic= new G4LogicalVolume(fScintillatorSolid, DSMaterial::Get()->GetBoronScintillator(),"Scintillator__Logic");
  G4VPhysicalVolume *fScintillatorPhysic  = new G4PVPlacement(0,
							   G4ThreeVector(0,0, 0),
							   "Scintillator_Physic",
							   fScintillatorLogic,
							   fMotherVolume,
							   false,
							   0,
							   myCheckOverlap);
  fScintillatorLogic->SetVisAttributes(G4VisAttributes::GetInvisible());
   // --------------------------- //
  // ------   Cryostats   ------ //
  // --------------------------- //
  
  //Load cryostat profiles                   
  // NOTE:  Each cryostat profile is evaluated in the ref frame with z = 0 at the level of the internal cryostat rib
  G4double myCryostatShiftZ = 0.; 
  // RANDOMLY CHOOSEN:
//  G4double myTPCShiftZ = 870.; 
  G4double LArGarIntefaceZ = 1740 ; 
  
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
  ifstream profileOuterVacCryo("../data/detector/ds5KVacCryoOuterProfile.dat") ; 
  profileOuterVacCryo >> myNumPointsOuterCryo;
  for( i = 0; i < myNumPointsOuterCryo; i++ ) { 
    profileOuterVacCryo >> myOuterCryostatZ[i]  >> myOuterCryostatRout[i];
    myOuterCryostatZ[i] += myCryostatShiftZ;
  }
  profileOuterVacCryo.close(); 

  // Outer Cryostat - Inner Surface 
  ifstream profileInnerVacCryo("../data/detector/ds5KVacCryoInnerProfile.dat") ; 
  profileInnerVacCryo  >> myNumPointsVacCryo;
  for( i = 0; i < myNumPointsVacCryo;   i++ ) {  
    profileInnerVacCryo >> myVacuumCryostatZ[i] >> myVacuumCryostatRout[i];
    myVacuumCryostatZ[i] += myCryostatShiftZ;
  }
  profileInnerVacCryo.close(); 

  // Inner Cryostat - Outer Surface 
  ifstream profileOuterCryo("../data/detector/ds5KCryoOuterProfile.dat") ; 
  profileOuterCryo  >>myNumPointsInnerCryo ;
  for( i = 0; i < myNumPointsInnerCryo;   i++ ) {  
    profileOuterCryo >> myInnerCryostatZ[i] >> myInnerCryostatRout[i];
    myInnerCryostatZ[i] += myCryostatShiftZ;
   // cout << myInnerCryostatRout[i] <<endl;
  }
  profileInnerVacCryo.close(); 
  // Inactive LAr and Gar - Outer Surface   - to be divided by 2
  ifstream profileInnerCryo("../data/detector/ds5KCryoInnerProfile.dat") ; 
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
  //      Parentesys         //
  //       Chimneyes         //
  //-------------------------//
  
  G4Tubs * fMainChimneySolid  = new G4Tubs("MainChimneySolid",   95.*mm , 100*mm, 220*mm, 0., myTwoPi ) ; 
  G4Tubs * fRightChimneySolid = new G4Tubs("RightChimneySolid", 30.*mm , 35*mm, 40*cm, 0., myTwoPi ) ; 
  G4Tubs * fLetfChimneySolid  = new G4Tubs("LeftChimneySolid", 35.*mm , 40*mm, 35*cm, 0., myTwoPi ) ; 
  G4Polycone      *fSolidDS5k_noCh   = new G4Polycone( "OuterCryostat_Solid_noCh", 0, myTwoPi, myNumPointsOuterCryo,myOuterCryostatZ, myRmin, myOuterCryostatRout  );
  G4UnionSolid    *fSolidDS5k_oneCh  = new G4UnionSolid("OuterCryostat1ch_Solid",fSolidDS5k_noCh,fMainChimneySolid, 0 , G4ThreeVector(0,0, 2215.7*mm) ) ; 
  G4UnionSolid    *fSolidDS5k_twoCh  = new G4UnionSolid("OuterCryostat2ch_Solid",fSolidDS5k_oneCh,fLetfChimneySolid, 0 , G4ThreeVector(-805*mm,0, 1830*mm) ) ; 
  G4UnionSolid    *fSolidDS5k_thrCh  = new G4UnionSolid("OuterCryostat3ch_Solid",fSolidDS5k_twoCh,fRightChimneySolid, 0 , G4ThreeVector(0*mm,905*mm, 1730*mm) ) ; 
  G4UnionSolid    *fSolidDS5k        = new G4UnionSolid("OuterCryostat_Solid",fSolidDS5k_thrCh,fRightChimneySolid, 0 , G4ThreeVector(905*mm,0, 1730*mm) ) ; 
  
  //-------------------------//
  //      Outer VacCryo      //
  //-------------------------//
    
  G4LogicalVolume *fLogicDS5k   = new G4LogicalVolume( fSolidDS5k, DSMaterial::Get()->GetStainlessSteel(), "OuterCryostat_Logic" );
  // the following shift has been roghly tuned, in order to place the center of the TPC in the center of the scintillator sphere 
  G4PVPlacement   *fPhysicDS5k  = new G4PVPlacement( 0, G4ThreeVector(0.,0.,-84.*cm), "OuterCryostat", fLogicDS5k,fScintillatorPhysic, false, 0, myCheckOverlap ); 
  //fLogicDS5k->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  G4Polycone      *fSolidDS5kVac   = new G4Polycone( "OuterCryostat_SolidVac", 0, myTwoPi, myNumPointsVacCryo,myVacuumCryostatZ, myRmin, myVacuumCryostatRout  );
  G4LogicalVolume *fLogicDS5kVac   = new G4LogicalVolume( fSolidDS5kVac, DSMaterial::Get()->GetVacuum(), "OuterCryostat_LogicVac" );
  G4PVPlacement   *fPhysicDS5kVac  = new G4PVPlacement( 0, myZeros, "OuterCryostatVac", fLogicDS5kVac, fPhysicDS5k, false, 0, myCheckOverlap );  
  G4VisAttributes *myVacuumAttributes = new G4VisAttributes(myWhite);

  //put lower part of chimneys in Vacuumm
/*  G4Tubs * fMainChimneySolidVac = new G4Tubs("MainChimneySolidVac", 95.*mm , 100*mm, 70*mm, 0., myTwoPi ) ; 
  G4Tubs * fRightChimneySolidVac = new G4Tubs("RightChimneySolidVac", 30.*mm , 35*mm, 60*mm, 0., myTwoPi ) ; 
  G4Tubs * fLetfChimneySolidVac = new G4Tubs("LeftChimneySolidVac", 35.*mm , 40*mm, 60*mm, 0., myTwoPi ) ; 
  G4LogicalVolume * fMainChimneyLogicVac = new G4LogicalVolume( fMainChimneySolidVac, DSMaterial::Get()->GetStainlessSteel(), "fMainChimneyLogicVac" );
  G4LogicalVolume * fRightChimneyLogicVac = new G4LogicalVolume( fRightChimneySolidVac, DSMaterial::Get()->GetStainlessSteel(), "fRightChimneyLogicVac" );
  G4LogicalVolume * fLetfChimneyLogicVac = new G4LogicalVolume( fLetfChimneySolidVac, DSMaterial::Get()->GetStainlessSteel(), "fLetfChimneyLogicVac" );
  G4PVPlacement * fMainChimneyPhysicVac  = new G4PVPlacement( 0,G4ThreeVector(0,0,(2115+70)*mm), "fMainChimneyPhysicVac", fMainChimneyLogicVac, fPhysicDS5kVac, false, 0, myCheckOverlap ); 
  G4PVPlacement * fRightChimneyPhysVac1 = new G4PVPlacement( 0,G4ThreeVector(905*mm,0,(1790)*mm), "fRightChimneyPhysicVac1", fRightChimneyLogicVac, fPhysicDS5kVac, false, 0, myCheckOverlap ); 
  G4PVPlacement * fRightChimneyPhysVac2 = new G4PVPlacement( 0,G4ThreeVector(0,905*mm,1790*mm), "fRightChimneyPhysicVac2", fRightChimneyLogicVac, fPhysicDS5kVac, false, 0, myCheckOverlap ); 
  G4PVPlacement * fLetfChimneyPhysVac = new G4PVPlacement( 0,G4ThreeVector(-805*mm,0,1890*mm), "fLeftChimneyPhysicVac", fLetfChimneyLogicVac, fPhysicDS5kVac, false, 0, myCheckOverlap ); 
*/  
  //-------------------------//
  //      Inner Cryo         //
  //-------------------------//

  G4Polycone      *fSolidInnerCryo   = new G4Polycone( "InnerCryo_Solid", 0, myTwoPi, myNumPointsInnerCryo,myInnerCryostatZ, myRmin,  myInnerCryostatRout );
  G4LogicalVolume *fLogicInnerCryo   = new G4LogicalVolume( fSolidInnerCryo, DSMaterial::Get()->GetStainlessSteel(), "SolidInnerCryo_Logic" );
  G4PVPlacement   *fPhysicInnerCryo  = new G4PVPlacement( 0, myZeros, "InnerCryo", fLogicInnerCryo,fPhysicDS5kVac , false, 0, myCheckOverlap );  
fLogicInnerCryo->SetVisAttributes(G4Colour  (1.0, 1.0, 1.0, 0.01));
  //  fLogicInnerCryo->SetVisAttributes(G4VisAttributes::GetInvisible());
 
  // evrything must be divided in 2 pieces in Z for LAr and GAr regions  
  G4Polycone      *fSolidInactiveLar   = new G4Polycone( "SolidInactiveLar_Solid", 0, myTwoPi, myNumPointsLiqArgon+1,myLiqArgonZ, myRmin, myLiqArgonRout  );
  G4LogicalVolume *fLogicInactiveLar   = new G4LogicalVolume( fSolidInactiveLar, DSMaterial::Get()->GetNSLiquidArgon(), "SolidInactiveLar_Logic" );
  G4PVPlacement   *fPhysicInactiveLar  = new G4PVPlacement( 0, myZeros, "SolidInactiveLar", fLogicInactiveLar, fPhysicInnerCryo, false, 0, myCheckOverlap );  
  G4VisAttributes *myLarAttributes = new G4VisAttributes(G4Colour  (0.0, .1, .9, 0.9));
  fLogicInactiveLar->SetVisAttributes(G4VisAttributes::GetInvisible());
  
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
  //(6mm thick 2.2cm high ring spaced 2.5cm and 71 rings total).

  //changed 
  //double myZRingBottom[2]          = {-myTPCHeightLAr/2., myTPCHeightLAr/2.-myTPCHeightGAr/2.};
  //double myZRingTop[2]             = {myTPCHeightLAr/2.-myTPCHeightGAr/2., myTPCHeightLAr/2.+myTPCHeightGAr/2.};
  
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
							   //G4ThreeVector(0,0,0.),
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

  //cout <<  -(myTPCHeight - myTPCHeightLAr)/2. << " " << myTPCHeight - (myTPCHeight - myTPCHeightLAr)/2. <<endl ;
  //cout << myTPCHeightLAr << " " << myTPCHeightGAr   <<endl ;

  //-------------------------//
  //    Teflon Reflector       //
  //-------------------------//
  
  // outwards from the inner active wall, a 1-inch thick PTFE reflector (octagonal shape).
  //double myZTeflonBottom[2]          = {0, myTPCHeightLAr};
  //double myZTeflonTop[2]             = {0, myTPCHeightGAr};
  double myZTeflonBottom[2]          = { LArGarIntefaceZ-myTPCHeight,LArGarIntefaceZ};
//  double myZTeflonBottom[2]          = {-myTPCHeightLAr/2., myTPCHeightLAr/2.-myTPCHeightGAr/2.};
  double myZTeflonTop[2]             = {myTPCHeightLAr/2.-myTPCHeightGAr/2., myTPCHeightLAr/2.+myTPCHeightGAr/2.};
  double myInnerRadiusTeflon[2] = {0,0};
  double myOuterRadiusTeflon[2] = {myTeflonOuterRadius,myTeflonOuterRadius};
  
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


  G4VPhysicalVolume *fPhysicTeflonBottom  = new G4PVPlacement(0,
							   G4ThreeVector(0,0,0),
							   "TeflonBottom",
							   fLogicTeflonBottom,
							   fPhysicRingBottom,
							   false,
							   0,
							   myCheckOverlap);

/*  
  G4Polyhedra *fSolidTeflonTop    = new G4Polyhedra( "TeflonTop_Solid",
        					    0,
        					    myTwoPi,
        					    8,
        					    2,
        					    myZTeflonTop,
        					    myInnerRadiusTeflon,
        					    myOuterRadiusTeflon );
  G4LogicalVolume *fLogicTeflonTop    = new G4LogicalVolume(fSolidTeflonTop,
				                         DSMaterial::Get()->GetTeflon(),
				                         "TeflonTop_Logic");

  G4VPhysicalVolume *fPhysicTeflonTop     = new G4PVPlacement(0,
							   G4ThreeVector(0,0,0),
							   "TeflonTop",
							   fLogicTeflonBottom,
							   fPhysicRingTop,
							   false,
							   0,
							   myCheckOverlap);
  fLogicTeflonBottom->SetVisAttributes(myRingTopAttributes);
 */ 
  //-------------------------//
  //          TPC            //
  //-------------------------//
  
  // TPC active volume section: Octagonal shape with edge 64.2cm when cold. (or Inscribed circle radius = 77.5cm)
  
  double myZTPCBottom[2]          = { LArGarIntefaceZ-myTPCHeight,LArGarIntefaceZ};
  double myInnerRadiusTPC[2] = {0,0};
  double myOuterRadiusTPC[2] = {myTPCOuterRadius,myTPCOuterRadius};
  
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
  G4VPhysicalVolume *fPhysicTPCBottom  = new G4PVPlacement(0,
							   G4ThreeVector(0,0,0),
							   "ActiveLAr",
							   fLogicActiveLAr,
							   fPhysicTeflonBottom,//fPhysicRingBottom,
							   false,
							   0,
							   myCheckOverlap);
  
  fLogicActiveLAr->SetVisAttributes(myBlue);
  G4Region* fLArRegion = new G4Region("LAr_Logic");
  fLogicActiveLAr->SetRegion( fLArRegion );
  fLArRegion->AddRootLogicalVolume(fLogicActiveLAr);


  //-------------------------//
  //  Active GAr Region      //
  //-------------------------//
    
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

  G4VPhysicalVolume *fPhysicGasPocket     = new G4PVPlacement(0,
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
  
  // Top/Bottom of the active volume section: Octagonal shape with edge 64.2cm when cold. (or Inscribed circle radius = 77.5cm)

  double myZTopWindow[2]          = {LArGarIntefaceZ+myGasPocketThickness,LArGarIntefaceZ+myGasPocketThickness+myWindowsThickness};
  //double myZTopWindow[2]          = {-myWindowsThickness/2. , myWindowsThickness/2.};
  double myZBotWindow[2]          = { LArGarIntefaceZ-myTPCHeight-myWindowsThickness, LArGarIntefaceZ-myTPCHeight};
//  double myZBotWindow[2]          = {
  double z1sipm = myZTopWindow[1]; 
  double z2sipm = myZTopWindow[1]+mySiPmBoardThickness+mySiPmThickness ; 
  double myZSiPmTop[2]             = { z1sipm, z2sipm};
  z1sipm =  myZBotWindow[0]-mySiPmBoardThickness-mySiPmThickness;
  z2sipm =  myZBotWindow[0];
  double myZSiPmBottom[2]          = { z1sipm, z2sipm};
  
  
  
  G4Polyhedra *fSolidTopWindow = new G4Polyhedra( "TopWindow_Solid",  0,
        					    myTwoPi,   8,
        					    2,  myZTopWindow,
        					    myInnerRadiusTPC,
        					    myOuterRadiusTPC );
 G4LogicalVolume *fLogicTopWindow   = new G4LogicalVolume(fSolidTopWindow,
				                         DSMaterial::Get()->GetFusedSilica(),
				                         "TopWindow_Logic");

  G4VPhysicalVolume *fPhysicTopWindow  = new G4PVPlacement(0,
							   G4ThreeVector(0,0, 0),
							   "TopWindow",  fLogicTopWindow,
							   fPhysicGasPocket,   false,
							   0,   myCheckOverlap);

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

  G4VPhysicalVolume *fPhysicBotWinwow    = new G4PVPlacement(0,
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

  G4VPhysicalVolume *fPhysicSiPmTop  = new G4PVPlacement(0,
							   G4ThreeVector(0,0,0),
							   "SiPMTop",
							   fLogicSiPMTop,
							   fPhysicGasPocket,//fPhysicRingBottom,
							   false,
							   0,
							   myCheckOverlap);
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

  G4VPhysicalVolume *fPhysicSiPmBottom  = new G4PVPlacement(0,
							   G4ThreeVector(0,0,0),
							   "SiPMBottom",
							   fLogicSiPMBottom,
							   fPhysicInactiveLar,//fPhysicRingBottom,
							   false,
							   0,
							   myCheckOverlap);
  

  DefineSurfaces(); 

}

DSDetectorDS5k::~DSDetectorDS5k(){
  ; //delete fMessenger;
}

void  DSDetectorDS5k::DefineSurfaces() {
 ;
}  

double DSDetectorDS5k::GetOctagonInnerRadius(double edge) {  return edge /2 *(1+sqrt(2)); }
double DSDetectorDS5k::GetOctagonOuterRadius(double edge) {  return edge /2. *sqrt(4  + 2 *sqrt(2)) ; }
