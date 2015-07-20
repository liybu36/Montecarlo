#include <iostream>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"

#include "DSDetectorDSG2.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSStorage.hh"
#include "DSParameters.hh"
#include "DSDetectorPMTDSG2.hh"
#include "G4VisAttributes.hh"

using namespace std;

////////////////////////////////////////////////
//////        Detector Description      ////////
////////////////////////////////////////////////
/* THESE NUMBERS ARE ALL WRONG

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


DSDetectorDSG2::DSDetectorDSG2(G4VPhysicalVolume *myMotherVolume) {
  
  fMotherVolume = myMotherVolume;

  const double myTwoPi = 2*M_PI*rad;
  G4bool   myCheckOverlap   = DSStorage::Get()->GetCheckOverlap();

  DSLog(routine) << " Constructing DSG2 Geometry" << endlog ;


  //G4double myCryostatShiftZ = DSParameters::Get()->GetCryostatShiftZ();
  G4ThreeVector myZeros( 0., 0., 0.);
  
  G4double epsilon = 0.1*mm;

  // --------------------------- //
  // ------   Cryostats   ------ //
  // --------------------------- //
  G4double cryoOuterR = 910.*mm;
  G4double cryoInnerR = 831.5*mm;
  G4double cryoOuterH = (1990.-303.3+240+229.5)*mm; //2235.56*mm;// 2000.*mm; // Cylindrical + top dome
  G4double cryoInnerH = (1858.-277.2+240+229.5)*mm; //2030.46*mm;//1843.1*mm; // Cylindrical + top dome
  G4double cryoOuterDomeH = 303.3*mm; //303.6*mm;
  G4double cryoInnerDomeH = 277.2*mm; //220.1*mm; //202.1*mm;

  G4double cryoOuterBodyH = cryoOuterH-cryoOuterDomeH;
  G4double cryoInnerBodyH = cryoInnerH-cryoInnerDomeH;
  G4double cryoOuterThickness = 5*mm;
  G4double cryoOuterDomeThickness = 9.5*mm;
  G4double cryoInnerThickness = 5*mm;
  G4double cryoInnerDomeThickness = 9.5*mm;

  G4ThreeVector outerCryoShift(0,0,cryoOuterBodyH/2.);
  G4ThreeVector innerCryoShift(0,0,cryoInnerBodyH/2.);
  G4ThreeVector vacShift(0,0,(cryoOuterBodyH-2*cryoOuterDomeThickness)/2.);
  G4ThreeVector larShift(0,0,(cryoInnerBodyH-2*cryoInnerDomeThickness)/2.);

  G4bool curvedBottom = true;

  //Outer Cryostat Layer
  fSolidCryoOuterTop = new G4Ellipsoid("CryoOuterTop_Solid", 
				       cryoOuterR,
				       cryoOuterR,
				       cryoOuterDomeH,
				       -epsilon,
				       cryoOuterDomeH);
  fSolidCryoOuterBody = new G4Tubs("CryOuterBody_Solid",
				  0.,
				  cryoOuterR,
				  cryoOuterBodyH/2.,
				  0,
				  2*M_PI);
  //Cryostat geometry with flat bottom
  fSolidCryoOuter = new G4UnionSolid("CryoOuter_Solid",
				     fSolidCryoOuterBody,
				     fSolidCryoOuterTop,
				     0,
				     outerCryoShift);

  fSolidCryoOuterBot = new G4Ellipsoid("CryoOuterBot_Solid",
				       cryoOuterR,
				       cryoOuterR,
				       cryoOuterDomeH,
				       -cryoOuterDomeH,
				       epsilon);
  //Cryostat geometry with curved bottom
  fSolidCryoOuterFull = new G4UnionSolid("CryoOuterFull_Solid",
					 fSolidCryoOuter,
					 fSolidCryoOuterBot,
					 0,
					 -outerCryoShift);

  if(curvedBottom)
    {
      fLogicCryoOuter = new G4LogicalVolume(fSolidCryoOuterFull,
					    DSMaterial::Get()->GetStainlessSteel(),
					    "CryoOuter_Logic");
    }
  else
    {
      fLogicCryoOuter = new G4LogicalVolume(fSolidCryoOuter,
					    DSMaterial::Get()->GetStainlessSteel(),
					    "CryoOuter_Logic");
    }
  fLogicCryoOuter->SetVisAttributes(new G4VisAttributes(G4Colour(0.2,0.2,0.2)));
  fPhysicCryoOuter = new G4PVPlacement(0,
				       myZeros,
				       "CryoOuter",
				       fLogicCryoOuter,
				       myMotherVolume,
				       false,
				       0,
				       myCheckOverlap);

  //Vacuum Layer
  G4double vacR = cryoOuterR-cryoOuterThickness;
  G4double vacDomeR = cryoOuterR-cryoOuterDomeThickness;
  G4double vacDomeH = cryoOuterDomeH-cryoOuterDomeThickness;
  G4double vacBodyH = cryoOuterBodyH-cryoOuterDomeThickness;
  fSolidVacTop = new G4Ellipsoid("VacTop_Solid", 
				 vacDomeR,
				 vacDomeR,
				 vacDomeH,
				 -epsilon,
				 vacDomeH);
  fSolidVacBody = new G4Tubs("VacBody_Solid",
			     0.,
			     vacR,
			     vacBodyH/2.,
			     0,
			     2*M_PI);
  //Vacuum layer for flat bottom
  fSolidVac = new G4UnionSolid("Vac_Solid",
			       fSolidVacBody,
			       fSolidVacTop,
			       0,
			       vacShift);

  fSolidVacBot = new G4Ellipsoid("VacBot_Solid",
				 vacDomeR,
				 vacDomeR,
				 vacDomeH,
				 -vacDomeH,
				 epsilon);
  //Vacuum layer for curved bottom
  fSolidVacFull = new G4UnionSolid("VacFull_Solid",
				   fSolidVac,
				   fSolidVacBot,
				   0,
				   -vacShift);

  if(curvedBottom)
    {
      fLogicVac = new G4LogicalVolume(fSolidVacFull,
				      DSMaterial::Get()->GetVacuum(),
				      "Vac_Logic");
    }
  else
    {
      fLogicVac = new G4LogicalVolume(fSolidVac,
				      DSMaterial::Get()->GetVacuum(),
				      "Vac_Logic");
    }
  fLogicVac->SetVisAttributes(new G4VisAttributes(G4Colour(1,1,1)));
  fPhysicVac = new G4PVPlacement(0,
				 myZeros,
				 "Vac",
				 fLogicVac,
				 fPhysicCryoOuter,
				 false,
				 0,
				 myCheckOverlap);

  //Inner Cryostat Layer
  fSolidCryoInnerTop = new G4Ellipsoid("CryoInnerTop_Solid", 
				       cryoInnerR,
				       cryoInnerR,
				       cryoInnerDomeH,
				       -epsilon,
				       cryoInnerDomeH);
  fSolidCryoInnerBody = new G4Tubs("CryoInnercBody_Solid",
				   0.,
				   cryoInnerR,
				   cryoInnerBodyH/2.,
				   0,
				   2*M_PI);
  //Inner cryostat for flat bottom
  fSolidCryoInner = new G4UnionSolid("CryoInner_Solid",
				     fSolidCryoInnerBody,
				     fSolidCryoInnerTop,
				     0,
				     innerCryoShift);

  fSolidCryoInnerBot = new G4Ellipsoid("CryoInnerBot_Solid",
				       cryoInnerR,
				       cryoInnerR,
				       cryoInnerDomeH,
				       -cryoInnerDomeH,
				       epsilon);
  //Inenr cryostat for curved bottom
  fSolidCryoInnerFull = new G4UnionSolid("CryoInnerFull_Solid",
					 fSolidCryoInner,
					 fSolidCryoInnerBot,
					 0,
					 -innerCryoShift);

  if(curvedBottom)
    {
      fLogicCryoInner = new G4LogicalVolume(fSolidCryoInnerFull,
					    DSMaterial::Get()->GetStainlessSteel(),
					    "CryoInner_Logic");
    }
  else
    {
      fLogicCryoInner = new G4LogicalVolume(fSolidCryoInner,
					    DSMaterial::Get()->GetStainlessSteel(),
					    "CryoInner_Logic");
    }
  fLogicCryoInner->SetVisAttributes(new G4VisAttributes(G4Colour(.2,.2,.2)));
  fPhysicCryoInner = new G4PVPlacement(0,
				       myZeros,
				       "CryoInner",
				       fLogicCryoInner,
				       fPhysicVac,
				       false,
				       0,
				       myCheckOverlap);

  // Fill cryostat with LAr
  G4double LArR = cryoInnerR-cryoInnerThickness;
  G4double LArDomeR = cryoInnerR-cryoInnerDomeThickness;
  G4double LArDomeH = cryoInnerDomeH-cryoInnerDomeThickness;
  G4double LArBodyH = cryoInnerBodyH-cryoInnerDomeThickness;
  fSolidLArTop = new G4Ellipsoid("LArTop_Solid", 
				 LArDomeR,
				 LArDomeR,
				 LArDomeH,
				 -epsilon,
				 LArDomeH);
  fSolidLArBody = new G4Tubs("LArBody_Solid",
			     0.,
			     LArR,
			     LArBodyH/2.,
			     0,
			     2*M_PI);
  //Liquid argon for flat bottom
  fSolidLAr = new G4UnionSolid("LAr_Solid",
			       fSolidLArBody,
			       fSolidLArTop,
			       0,
			       larShift);

  fSolidLArBot = new G4Ellipsoid("LArBot_Solid",
				 LArDomeR,
				 LArDomeR,
				 LArDomeH,
				 -LArDomeH,
				 epsilon);
  //Liquid argon for curved bottom
  fSolidLArFull = new G4UnionSolid("LArFull_Solid",
				   fSolidLAr,
				   fSolidLArBot,
				   0,
				   -larShift);

  if(curvedBottom)
    {
      fLogicLAr = new G4LogicalVolume(fSolidLArFull,
				      DSMaterial::Get()->GetLiquidArgon(),
				      "LAr_Logic");
    }
  else
    {
      fLogicLAr = new G4LogicalVolume(fSolidLAr,
				      DSMaterial::Get()->GetLiquidArgon(),
				      "LAr_Logic");
    }
  fLogicLAr->SetVisAttributes(new G4VisAttributes(G4Colour(0,.1,.9)));
  fPhysicLAr = new G4PVPlacement(0,
				 myZeros,
				 "LAr",
				 fLogicLAr,
				 fPhysicCryoInner,
				 false,
				 0,
				 myCheckOverlap);

  // The definition of the LAr region is needed to set the range cuts for this volume 
  // to a smaller value with respect to the rest of the detector ( see DSPhysicsList::SetCuts() )
  G4Region* fLArRegion = new G4Region("LAr_Logic");
  fLogicLAr->SetRegion( fLArRegion );
  fLArRegion->AddRootLogicalVolume(fLogicLAr);


  /*
  // GAr layer
  G4double GArH = 150.*mm;
  G4double GArShiftZ = LArBodyH/2.;//LArDomeH+LArBodyH/2-GArH-cryoInnerThickness;
  G4ThreeVector GArShift(0, 0, GArShiftZ);
  fSolidGAr = new G4Ellipsoid("GAr_Solid",
			      LArDomeR,
			      LArDomeR,
			      LArDomeH,
			      LArDomeH-GArH,
			      LArDomeH);
  fLogicGAr = new G4LogicalVolume(fSolidGAr,
				  DSMaterial::Get()->GetGaseousArgon(),
				  "GAr_Logic");
  fLogicGAr->SetVisAttributes(new G4VisAttributes(G4Colour(.4,.4,0)));
  fPhysicGAr = new G4PVPlacement(0,
				 larShift,
				 "GAr",
				 fLogicGAr,
				 fPhysicLAr,
				 false,
				 0,
				 myCheckOverlap);
  */
  //Teflon Insulator and reflector
  G4double teflonThickness = 30.*mm;
  G4double teflonIR = 1483.*mm/2.;
  G4double teflonOR = teflonIR + teflonThickness;
  G4double teflonH = 1350.*mm + 229.5*mm; //1110*mm; //1350.*mm;
  G4double fsThickness = 12.7*mm;
  G4double lowerDeadSpace = 150.*mm/2;
  G4double teflonShiftZ = teflonH/2.+lowerDeadSpace+fsThickness-LArBodyH/2.;
  G4ThreeVector teflonShift(0, 0, teflonShiftZ);
  fSolidTeflon = new G4Tubs("Teflon_Solid",
			    teflonIR,
			    teflonOR,
			    teflonH/2.,
			    0,
			    2*M_PI);
  fLogicTeflon = new G4LogicalVolume(fSolidTeflon,
				     DSMaterial::Get()->GetTeflon(),
				     "Teflon_Logic");
  fLogicTeflon->SetVisAttributes(new G4VisAttributes(G4Colour(1,1,1)));
  fPhysicTeflon = new G4PVPlacement(0,
				    teflonShift,
				    "Teflon",
				    fLogicTeflon,
				    fPhysicLAr,
				    false,
				    0,
				    myCheckOverlap);

  //Field Cage Rings
  G4double ringThickness = 5.*mm;
  G4double ringR = teflonIR + 15*mm;
  fSolidRings = new G4Tubs("Rings_Solid",
			   ringR,
			   ringR + ringThickness,
			   teflonH/2.,
			   0,
			   2*M_PI);
  fLogicRings = new G4LogicalVolume(fSolidRings,
				    DSMaterial::Get()->GetMetalCopper(),
				    "Rings_Logic");
  fLogicRings->SetVisAttributes(new G4VisAttributes(G4Colour(1,0,0)));
  fPhysicRings = new G4PVPlacement(0,
				   myZeros,
				   "FiledCageRings",
				   fLogicRings,
				   fPhysicTeflon,
				   false,
				   0,
				   myCheckOverlap);

  //Lower TPC Window
  G4double lowerWindowShiftZ = teflonShiftZ-teflonH/2.-fsThickness/2.;
  G4ThreeVector lowerWindowShift(0,0,lowerWindowShiftZ);
  fSolidLowerWindow = new G4Tubs("LowerWindow_Solid",
				 0,
				 teflonOR,
				 fsThickness/2.,
				 0,
				 2*M_PI);
  fLogicLowerWindow = new G4LogicalVolume(fSolidLowerWindow,
					  DSMaterial::Get()->GetFusedSilica(),
					  "LowerWindow_Logic");
  fPhysicLowerWindow = new G4PVPlacement(0,
					 lowerWindowShift,
					 "LowerFSWindow",
					 fLogicLowerWindow,
					 fPhysicLAr,
					 false,
					 0,
					 myCheckOverlap);

  //Upper TPC Window
  G4double fsTopThickness = 25.4*mm;
  G4double upperWindowShiftZ = teflonShiftZ+teflonH/2.+fsTopThickness/2.;
  G4ThreeVector upperWindowShift(0,0,upperWindowShiftZ);
  fSolidUpperWindow = new G4Tubs("UpperWindow_Solid",
				 0,
				 teflonOR,
				 fsTopThickness/2.,
				 0,
				 2*M_PI);
  fLogicUpperWindow = new G4LogicalVolume(fSolidUpperWindow,
					  DSMaterial::Get()->GetFusedSilica(),
					  "UpperWindow_Logic");
  fPhysicUpperWindow = new G4PVPlacement(0,
					 upperWindowShift,
					 "UpperFSWindow",
					 fLogicUpperWindow,
					 fPhysicLAr,
					 false,
					 0,
					 myCheckOverlap);
  

  //Gas Pocket
  G4double gasThickness = 12.7*mm;
  G4ThreeVector gasPocketShift(0,0,(gasThickness-fsTopThickness)/2.);

  fSolidGasPocket = new G4Tubs("GasPocket_Solid",
			       0,
			       teflonIR,
			       gasThickness/2.,
			       0,
			       2*M_PI);
  fLogicGasPocket = new G4LogicalVolume(fSolidGasPocket,
					DSMaterial::Get()->GetGaseousArgon(),
					"GasPocket_Logic");
  fLogicGasPocket->SetVisAttributes(new G4VisAttributes(G4Colour(.4,.4,0)));
  fPhysicGasPocket = new G4PVPlacement(0,
				       gasPocketShift,
				       "GasPocket",
				       fLogicGasPocket,
				       fPhysicUpperWindow,
				       false,
				       0,
				       myCheckOverlap);


  // Set the z coordinate of the LAr - GAr interface, necessary for S2 generation in DSLightX
  G4double myLArGArBoundaryPosZ = ( gasPocketShift + upperWindowShift ).z() - gasThickness/2.;
  DSStorage::Get()->SetLArGArBoundaryPosZ( myLArGArBoundaryPosZ );

  

  // PMT Assembly Mother Volume & Teflon Cap Disk
  G4double myPMTAssembly_h = 130.*mm;//DSParameters::Get()->GetPMTAssemblyHeight()*mm;       
  G4double myTeflonCap_d   = 2*teflonOR; //DSParameters::Get()->GetTeflonCapDiameter()*mm;
  G4double myTeflonCap_h   = DSParameters::Get()->GetTeflonCapHeight()*mm;
          
  G4ThreeVector myPMTAssemblyTop( 0, 0,    upperWindowShiftZ + fsTopThickness/2. + myPMTAssembly_h/2. );
  G4ThreeVector myPMTAssemblyBottom( 0, 0, lowerWindowShiftZ - fsThickness/2. - myPMTAssembly_h/2. ); 
  G4ThreeVector myTeflonCapPos( 0, 0, -myPMTAssembly_h/2. + myTeflonCap_h/2. );

  G4RotationMatrix* rotX180 = new G4RotationMatrix;
  rotX180->rotateX( M_PI*rad );


  fSolidPMTAssemblyTub      = new G4Tubs( "PMTAssemblyTub_Solid", 0, myTeflonCap_d/2., myPMTAssembly_h/2., 0, myTwoPi );
  fSolidPMTAssemblyTop      = new G4IntersectionSolid("PMTAssemblyTop_Solid", fSolidPMTAssemblyTub, fSolidLArFull, 0, -1*myPMTAssemblyTop);
  fSolidPMTAssemblyBottom   = new G4IntersectionSolid("PMTAssemblyBottom_Solid",fSolidPMTAssemblyTub,fSolidLArFull,0, -1*myPMTAssemblyBottom);

  fLogicPMTAssemblyTop      = new G4LogicalVolume( fSolidPMTAssemblyTop, DSMaterial::Get()->GetLiquidArgon(), "PMTAssemblyTop_Logic" ); 
  fLogicPMTAssemblyTop->SetVisAttributes(new G4VisAttributes(G4Colour(0,.3,.8)));
  fPhysicPMTAssemblyTop     = new G4PVPlacement( 0, myPMTAssemblyTop,    "PMTAssemblyTop",    fLogicPMTAssemblyTop, fPhysicLAr, false, 0, myCheckOverlap );

  fLogicPMTAssemblyBottom   = new G4LogicalVolume( fSolidPMTAssemblyBottom, DSMaterial::Get()->GetLiquidArgon(), "PMTAssemblyBottom_Logic" ); 
  fLogicPMTAssemblyBottom->SetVisAttributes(new G4VisAttributes(G4Colour(0,.3,.8)));
  fPhysicPMTAssemblyBottom  = new G4PVPlacement( rotX180, myPMTAssemblyBottom, "PMTAssemblyBottom", fLogicPMTAssemblyBottom, fPhysicLAr, false, 0, myCheckOverlap );

  fSolidTeflonCapDisk   = new G4Tubs( "TeflonCapDisk_Solid", 0, myTeflonCap_d/2., myTeflonCap_h/2., 0, myTwoPi );

  fLogicTeflonCapDiskTop   = new G4LogicalVolume( fSolidTeflonCapDisk, DSMaterial::Get()->GetTeflon(), "TeflonCapDiskTop_Logic" );
  fPhysicTeflonCapDiskTop  = new G4PVPlacement( 0, myTeflonCapPos, "TeflonCapDiskTop", fLogicTeflonCapDiskTop, fPhysicPMTAssemblyTop, false, 0, myCheckOverlap );
  
  fLogicTeflonCapDiskBottom   = new G4LogicalVolume( fSolidTeflonCapDisk, DSMaterial::Get()->GetTeflon(), "TeflonCapDiskBottom_Logic" );
  fPhysicTeflonCapDiskBottom  = new G4PVPlacement(0, myTeflonCapPos, "TeflonCapDiskBottom", fLogicTeflonCapDiskBottom, fPhysicPMTAssemblyBottom, false, 0, myCheckOverlap );

  
  // PMTs
  // --- Note: fPhysicInnerLiqArgon is passed as argument just for the definition of the surface between the photocathode and the LAr
  DSDetectorPMTDSG2( (G4VPhysicalVolume*) fPhysicPMTAssemblyTop,    (G4VPhysicalVolume*) fPhysicTeflonCapDiskTop );
  DSDetectorPMTDSG2( (G4VPhysicalVolume*) fPhysicPMTAssemblyBottom, (G4VPhysicalVolume*) fPhysicTeflonCapDiskBottom );

  /*
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

  //G4ThreeVector myReflectorPos = myCathodeWinPos + G4ThreeVector( 0, 0, 2*myITOThickness + myCathodeWindow_h/2.);
  G4ThreeVector myReflectorPos = myCathodeWinPos + G4ThreeVector( 0, 0, myCathodeWindow_h/2.);
  G4ThreeVector myGasPocketPos( 0, 0, myReflector_h + myAboveGrid_h - myTPBThickness - myGasGap/2.);

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


  fSolidReflector   = new G4Polycone( "Reflector_Solid", 0, myTwoPi, myNZR, myReflectorProf_R, myReflectorProf_Z );
  fLogicReflector   = new G4LogicalVolume( fSolidReflector, DSMaterial::Get()->GetTeflon(), "Reflector_Logic" );
  fPhysicReflector  = new G4PVPlacement( 0, myReflectorPos, "Reflector", fLogicReflector, fPhysicInnerLiqArgon, false, 0, myCheckOverlap );

  fSolidTPB         = new G4Polycone( "TPB_Solid", 0, myTwoPi, myNZR_TPB, myReflectorTPBProf_R, myReflectorTPBProf_Z );
  fLogicTPB         = new G4LogicalVolume( fSolidTPB, DSMaterial::Get()->GetTPB(), "TPB_Logic" );
  fPhysicTPB        = new G4PVPlacement( 0, myZeros, "TPB", fLogicTPB, fPhysicReflector, false, 0, myCheckOverlap );

  fSolidActiveLAr   = new G4Polycone( "ActiveLAr_Solid", 0, myTwoPi, myNZR_Ar, myReflectorArProf_R, myReflectorArProf_Z );
  fLogicActiveLAr   = new G4LogicalVolume( fSolidActiveLAr, DSMaterial::Get()->GetLiquidArgon(), "ActiveLAr_Logic" );
  fPhysicActiveLAr  = new G4PVPlacement( 0, myZeros, "ActiveLAr", fLogicActiveLAr, fPhysicTPB, false, 0, myCheckOverlap );

  fSolidGasPocket   = new G4Tubs( "GasPocket_Solid", 0, myAboveGrid_id/2. - myTPBThickness, myGasGap/2., 0, myTwoPi );
  fLogicGasPocket   = new G4LogicalVolume( fSolidGasPocket, DSMaterial::Get()->GetGaseousArgon(), "GasPocket_Logic" );
  fPhysicGasPocket  = new G4PVPlacement( 0, myGasPocketPos, "GasPocket", fLogicGasPocket, fPhysicActiveLAr, false, 0, myCheckOverlap );



  // Grid
  G4double myGrid_h = 0.01 * cm;
  G4ThreeVector myGridPosition( 0, 0, myReflectorArProf_Z[ 2 ] + myGrid_h/2. );

  fSolidGrid  = new G4Tubs( "Grid_Solid", 0, myReflectorArProf_R[ 2 ], myGrid_h/2., 0, myTwoPi ); 
  fLogicGrid  = new G4LogicalVolume( fSolidGrid, DSMaterial::Get()->GetGridSteel(), "Grid_Logic" );
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
  */

  DefineSurfaces(); 

}

DSDetectorDSG2::~DSDetectorDSG2(){
  ; //delete fMessenger;
}

void  DSDetectorDSG2::DefineSurfaces() {

  ////////////////////////////////////////
  // TPC - BScint
  ////////////////////////////////////////
  G4OpticalSurface *fOpElectropolishedStainlessSteelSurface = new G4OpticalSurface("OpElectropolishedStainlessSteelSurface");
  fOpElectropolishedStainlessSteelSurface->SetType(dielectric_metal);
  fOpElectropolishedStainlessSteelSurface->SetModel(unified);
  fOpElectropolishedStainlessSteelSurface->SetFinish(groundbackpainted);
  fOpElectropolishedStainlessSteelSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetElectropolishedStainlessSteelMPT()); 
  //  fOpElectropolishedStainlessSteelSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetLumirrorMPT()); 
  fDSG2OuterSurface = new G4LogicalBorderSurface("DSG2OuterSurface",
						 fMotherVolume,
						 fPhysicCryoOuter,
						 fOpElectropolishedStainlessSteelSurface);

  fDSG2OuterSurface2 = new G4LogicalBorderSurface("DSG2OuterSurface2",
						 fPhysicCryoOuter,
						 fMotherVolume,
						 fOpElectropolishedStainlessSteelSurface);


  /*
  ////////////////////////////////////////
  // GAr - LAr
  ////////////////////////////////////////
  G4OpticalSurface *fOpGArLArSurface = new G4OpticalSurface("OpGArLArSurface");
  new G4LogicalBorderSurface("GArLArSurface", fPhysicGasPocket, fPhysicActiveLAr, fOpGArLArSurface); 
  fOpGArLArSurface->SetType( dielectric_dielectric );
  fOpGArLArSurface->SetModel( unified );
  fOpGArLArSurface->SetFinish( polished );
  G4MaterialPropertiesTable *fGArLArSurfProp = new G4MaterialPropertiesTable();
  fOpGArLArSurface->SetMaterialPropertiesTable( fGArLArSurfProp );


  ////////////////////////////////////////
  // TPB - GAr
  ////////////////////////////////////////
  G4OpticalSurface *fOpTPBGArSurface = new G4OpticalSurface("OpTPBGArSurface");
  new G4LogicalBorderSurface("TPBGArSurface", fPhysicGasPocket, fPhysicTPB, fOpTPBGArSurface );
  fOpTPBGArSurface->SetType( dielectric_dielectric );
  fOpTPBGArSurface->SetModel( unified );
  fOpTPBGArSurface->SetFinish( ground );
  G4MaterialPropertiesTable *fTPBGArSurfProp = new G4MaterialPropertiesTable();
  fOpTPBGArSurface->SetMaterialPropertiesTable( fTPBGArSurfProp );

  ////////////////////////////////////////
  // ITO /////
  ////////////////////////////////////////
			 
  G4MaterialPropertiesTable *fITOSurfProp = new G4MaterialPropertiesTable();
  fITOSurfProp->AddConstProperty("DOITO",1);
  
  ////////////////////////////////////////
  // BellTop (fusedsilica) - TPB and CathodeWindow - TPB
  ////////////////////////////////////////
  G4OpticalSurface *fOpWindowTPBSurface     = new G4OpticalSurface("OpWindowTPBSurface");
  new G4LogicalBorderSurface("TopWindowTPBSurface",    fPhysicBellTop,  fPhysicTPB, fOpWindowTPBSurface );
  new G4LogicalBorderSurface("BottomWindowTPBSurface", fPhysicCathodeWindow, fPhysicTPB, fOpWindowTPBSurface );
  fOpWindowTPBSurface->SetType( dielectric_dielectric );
  fOpWindowTPBSurface->SetModel( unified );
  fOpWindowTPBSurface->SetFinish( ground );
  G4MaterialPropertiesTable *faWindowTPBSurfProp = new G4MaterialPropertiesTable();
  fOpWindowTPBSurface->SetMaterialPropertiesTable( fITOSurfProp );
  
  
  ////////////////////////////////////////
  // BellTop - LAr and CathodeWindow - LAr 
  ////////////////////////////////////////
  G4OpticalSurface *fOpWindowLArSurface     = new G4OpticalSurface("OpWindowLArSurface");
  new G4LogicalBorderSurface("TopWindowLArSurface",    fPhysicBellTop,       fPhysicInnerLiqArgon, fOpWindowLArSurface); 
  new G4LogicalBorderSurface("BottomWindowLArSurface", fPhysicCathodeWindow, fPhysicInnerLiqArgon, fOpWindowLArSurface); 
  fOpWindowLArSurface->SetType( dielectric_dielectric );
  fOpWindowLArSurface->SetModel( unified );
  fOpWindowLArSurface->SetFinish( ground );
  G4MaterialPropertiesTable *fWLArSurfProp = new G4MaterialPropertiesTable();
  fOpWindowLArSurface->SetMaterialPropertiesTable( fITOSurfProp );


  ////////////////////////////////////////
  // InnerLAr - LAr ???
  ////////////////////////////////////////

  ////////////////////////////////////////
  // Teflon (Reflector) - TPB
  ////////////////////////////////////////
  G4OpticalSurface *fOpTPBTeflonSurface = new G4OpticalSurface("OpTBPTeflonSurface");
  new G4LogicalBorderSurface("TPBTeflonSurface", fPhysicTPB, fPhysicReflector, fOpTPBTeflonSurface );
  fOpTPBTeflonSurface->SetType( dielectric_metal );
  fOpTPBTeflonSurface->SetModel(unified);
  fOpTPBTeflonSurface->SetFinish(groundfrontpainted);                           
  fOpTPBTeflonSurface->SetSigmaAlpha(0.1);                          
  G4MaterialPropertiesTable *fTPBTeflonSurfProp = new G4MaterialPropertiesTable();
  G4double TeflonTPBENE[2] = {0.1*eV, 20.0*eV};
  G4double TeflonTPBREF[2] = {0.99, 0.99};
  G4double TeflonTPBEFF[2] = {0.00, 0.00};
  fTPBTeflonSurfProp->AddProperty("REFLECTIVITY", TeflonTPBENE, TeflonTPBREF, 2);			 
  fTPBTeflonSurfProp->AddProperty("EFFICIENCY",   TeflonTPBENE, TeflonTPBEFF, 2);		       
  fOpTPBTeflonSurface->SetMaterialPropertiesTable( fTPBTeflonSurfProp );

  ////////////////////////////////////////
  // Teflon (fPhysicTeflonCapDiskBottom and fPhysicTeflonCapDiskTop) - LiquidArgon
  ////////////////////////////////////////
  G4OpticalSurface *fOpLArCapDiskSurface = new G4OpticalSurface("OpLArCapDiskSurface");
  new G4LogicalBorderSurface("InnerLArCapDiskTopSurface",fPhysicInnerLiqArgon ,fPhysicTeflonCapDiskBottom , fOpLArCapDiskSurface );
  new G4LogicalBorderSurface("InnerLArCapDiskBottomSurface",fPhysicInnerLiqArgon ,fPhysicTeflonCapDiskTop , fOpLArCapDiskSurface );
  fOpLArCapDiskSurface->SetType( dielectric_metal );
  fOpLArCapDiskSurface->SetModel(unified);
  fOpLArCapDiskSurface->SetFinish(groundfrontpainted);                           
  fOpLArCapDiskSurface->SetSigmaAlpha(0.1);                          
  G4MaterialPropertiesTable *fLArCapDiskSurfProp = new G4MaterialPropertiesTable();
  G4double TeflonLArENE[2] = {0.1*eV, 20.0*eV};
  G4double TeflonLArREF[2] = {0.99, 0.99};
  G4double TeflonLArEFF[2] = {0.00, 0.00};
  fLArCapDiskSurfProp->AddProperty("REFLECTIVITY", TeflonLArENE, TeflonLArREF, 2);			 
  fLArCapDiskSurfProp->AddProperty("EFFICIENCY",   TeflonLArENE, TeflonLArEFF, 2);		       
  fOpLArCapDiskSurface->SetMaterialPropertiesTable( fLArCapDiskSurfProp );


  ////////////////////////////////////////
  // TPB - LAr
  ////////////////////////////////////////
  G4OpticalSurface *fOpTPBLArSurface = new G4OpticalSurface("OpTBPLArSurface");
  new G4LogicalBorderSurface("TPBLArSurface", fPhysicActiveLAr, fPhysicTPB, fOpTPBLArSurface );
  fOpTPBLArSurface->SetType( dielectric_dielectric );
  fOpTPBLArSurface->SetModel( unified );
  fOpTPBLArSurface->SetFinish( polished );                           
  fOpTPBLArSurface->SetSigmaAlpha(0.0);                          
  G4MaterialPropertiesTable *fTPBLArSurfProp = new G4MaterialPropertiesTable();
  fOpTPBLArSurface->SetMaterialPropertiesTable( fTPBLArSurfProp );

  
  ////////////////////////////////////////
  // Teflon - LAr
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
  fLArTeflonSurfProp->AddConstProperty("REFLECTIVITY", 0.2);
  fOpLArTeflonSurface->SetMaterialPropertiesTable( fLArTeflonSurfProp );
  
  
  
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
  */
}  
