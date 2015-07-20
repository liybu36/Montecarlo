#include <iostream>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4PVPlacement.hh"
#include "G4UIcommand.hh"

#include "DSDetectorPMTDSG3.hh"
#include "DSStorage.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"




using namespace std;


DSDetectorPMTDSG3::DSDetectorPMTDSG3( G4VPhysicalVolume* myMotherVolume, G4VPhysicalVolume* myCapDiskVolume ){

  fMotherVolume  = myMotherVolume;
  fCapDiskVolume = myCapDiskVolume;


  const double myTwoPi    = 2*M_PI*rad;
  G4ThreeVector myZeros(0,0,0);
  G4bool   myCheckOverlap = DSStorage::Get()->GetCheckOverlap();
  


  // Variables to distinguish bewteen top and bottom
  G4int  myOffset = 379;
  if( myMotherVolume->GetRotation() == 0 ) {
    myOffset = 0;
    DSLog(routine) << " Constructing Top TPC PMTs Geometry  " << endlog;
  }
  else DSLog(routine) << " Constructing Bottom TPC PMTs Geometry  " << endlog;


  G4ThreeVector myDir( 0, 0, 1 );
  G4double myPMTBody_hh = myCapDiskVolume->GetLogicalVolume()->GetSolid()->DistanceToOut( myZeros, myDir );    // Cap Disk Half Height
  G4double myPMTBody_d = 5*2.54*cm;
  G4double myPMTWallThickness = 0.55*mm;
  G4double myPhotocathodeThickness = 3.5866*mm;


  // The PMT model is assumed to be the SiGHT.
  // It is modeled as the sum of a cylinder and a hemisphere, as the vacuum volume inside it.


  // PMT Body
  fSolidPMTBody = new G4Tubs( "PMTBody_Solid", 0, myPMTBody_d/2., myPMTBody_hh, 0, myTwoPi);
  fLogicPMTBody = new G4LogicalVolume( fSolidPMTBody, DSMaterial::Get()->GetFusedSilica(), "PMTBody_Logic");


  fSolidPMTBodyVac = new G4Tubs( "PMTBodyVac_Solid", 0, myPMTBody_d/2. - myPMTWallThickness, myPMTBody_hh - myPMTWallThickness/2., 0, myTwoPi);
  fLogicPMTBodyVac = new G4LogicalVolume( fSolidPMTBodyVac, DSMaterial::Get()->GetVacuum(), "PMTBodyVac_Logic");


  // PMT Head
  fSolidPMTHead = new G4Sphere( "PMTHead_Solid", 0, myPMTBody_d/2., 0, myTwoPi, myTwoPi/4., myTwoPi/2 );
  fLogicPMTHead = new G4LogicalVolume( fSolidPMTHead, DSMaterial::Get()->GetBialkali(), "PMTHead_Logic");


  fSolidPMTHeadVac = new G4Sphere( "PMTHeadVac_Solid", 0, myPMTBody_d/2. - myPhotocathodeThickness, 0, myTwoPi, myTwoPi/4., myTwoPi/2 );
  fLogicPMTHeadVac = new G4LogicalVolume( fSolidPMTHeadVac, DSMaterial::Get()->GetVacuum(), "PMTHeadVac_Logic");
   

  // IMPORTANT: set here the index of the cathode material
  DSStorage::Get()->SetPMTMaterialIndex(fLogicPMTHead->GetMaterial()->GetIndex());



  // Positioning
  //
  //   The PMT array is composed by 23 lines, the i-th containing myNpmts[i] PMTs.
  //   The array is constructed starting from the y coordinate of the first PMT in each line.
  //
  G4double myCenterToCenterDist = 5.5 * 2.54 * cm;
  G4double myXaxisUnit = myCenterToCenterDist * pow( 3, 0.5) / 2.;
  G4double myYaxisUnit = myCenterToCenterDist / 2.;

  G4double myNpmts[12] = { 8, 11, 14, 15, 16, 17, 18, 19, 20, 21, 20, 21};                // The pattern is symmetric, so only half of it is listed here
  G4double myXpmt00 = -11;                                                                // x coord of the first PMT in the first line, in x axis units
  G4double myYpmt0[12] = { -7, -10, -13, -14, -15, -16, -17, -18, -19, -20, -19, -20 };   // y coord of the first PMT in each line, in y axis units

  G4int myLineIdx = 0;    // index for selecting the line's npmts and ypmt0. 0 ---> 11 ---> 0 
  G4int myLine = 0;       // Line number
  G4int myPMTIdx = 0;     // PMT#

  G4double myMotherVolume_hh = myMotherVolume->GetLogicalVolume()->GetSolid()->DistanceToOut( myZeros, myDir ); 
  G4ThreeVector myPMTBodyVacPos( 0, 0, -myPMTWallThickness/2. );
  G4ThreeVector myPMTHeadZ( 0, 0, myMotherVolume_hh - 2*myPMTBody_hh ); 



  // Actual placement
  for(int idx = 0; idx < 23; idx++){

    idx < 12 ? myLineIdx = idx : myLineIdx--;

    for(int pmt = 0; pmt < myNpmts[myLineIdx]; pmt++){
 
      G4ThreeVector myPMTBodyPos( (myXpmt00 + myLine)*myXaxisUnit, (myYpmt0[myLineIdx] + 2*pmt) * myYaxisUnit, 0 );
      G4ThreeVector myPMTHeadPos= myPMTBodyPos + myPMTHeadZ; 

      fPhysicPMTBody[myPMTIdx] = new G4PVPlacement( 0, myPMTBodyPos, "PMTBody", fLogicPMTBody, myCapDiskVolume, false, 0, myCheckOverlap );

      G4String myPMTHeadName = "TPMT_" + G4UIcommand::ConvertToString( myPMTIdx + myOffset );

      fPhysicPMTHead[myPMTIdx] = new G4PVPlacement( 0, myPMTHeadPos, myPMTHeadName, fLogicPMTHead, myMotherVolume, false, 0, myCheckOverlap );

      myPMTIdx++;
    }


    myLine++;
  }

  fPhysicPMTBodyVac = new G4PVPlacement( 0, myPMTBodyVacPos, "PMTBodyVac", fLogicPMTBodyVac, fPhysicPMTBody[0], false, 0, myCheckOverlap);
  fPhysicPMTHeadVac = new G4PVPlacement( 0, myZeros, "PMTHeadVac", fLogicPMTHeadVac, fPhysicPMTHead[0], false, 0, myCheckOverlap );

  

  DefineSurfaces();

}


DSDetectorPMTDSG3::~DSDetectorPMTDSG3(){
  ;
}


void  DSDetectorPMTDSG3::DefineSurfaces(){

  // Photocathode - Vacuum
  fOpPMTVacuumSurface = new G4OpticalSurface("OpPMTLArSurface");
  fPMTVacuumSurface = new G4LogicalBorderSurface("PMTVacuumSurface", fPhysicPMTHeadVac, fPhysicPMTHead[0], fOpPMTVacuumSurface ); 
  fOpPMTVacuumSurface->SetType( dielectric_metal );
  fOpPMTVacuumSurface->SetModel( glisur );
  fOpPMTVacuumSurface->SetFinish( polished );  
  fPMTVacuumSurfProp = new G4MaterialPropertiesTable();
  fPMTVacuumSurfProp->AddConstProperty("REFLECTIVITY", 1.0);
  fPMTVacuumSurfProp->AddConstProperty("EFFICIENCY",   0.0);
  fOpPMTVacuumSurface->SetMaterialPropertiesTable( fPMTVacuumSurfProp );
  

  // Photocathode - LAr
  fOpPMTLArSurface = new G4OpticalSurface("OpPMTLArSurface");
  for(int i = 0; i < 379; i++)  fPMTLArSurface[i] = new G4LogicalBorderSurface("PMTLArSurface", fMotherVolume, fPhysicPMTHead[i], fOpPMTLArSurface );  
  fOpPMTLArSurface->SetType( dielectric_dielectric );
  fOpPMTLArSurface->SetModel( unified );
  fOpPMTLArSurface->SetFinish( polished );  
  fPMTLArSurfProp = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable *fLArPMTSurfProp = new G4MaterialPropertiesTable();
  G4double PMTLArENE[4] = {0.1*eV, 8.0*eV, 8.3*eV, 20.0*eV};
  G4double TAREFUVK  = DSParameters::Get()->GetPMTLArUVRef();
  G4double TAREFVIK = DSParameters::Get()->GetPMTLArVisRef();
  G4double PMTLArREF[4] = {TAREFVIK, TAREFVIK ,TAREFUVK , TAREFUVK };
  fLArPMTSurfProp->AddProperty("REFLECTIVITY", PMTLArENE, PMTLArREF, 4);			 
  fOpPMTLArSurface->SetMaterialPropertiesTable( fPMTLArSurfProp );
  


  // Teflon - LAr
  fOpTeflonLArSurface = new G4OpticalSurface("OpTeflonLArSurface");
  fTeflonLArSurface = new G4LogicalBorderSurface("TeflonLArSurface", fCapDiskVolume, fMotherVolume, fOpTeflonLArSurface );
  fOpTeflonLArSurface->SetType(dielectric_dielectric);
  fOpTeflonLArSurface->SetModel(unified);
  fOpTeflonLArSurface->SetSigmaAlpha(0.1);
  //where sigma_alpha is in [rad]
  fOpTeflonLArSurface->SetFinish(polished);

  G4MaterialPropertiesTable *fLArTeflonSurfProp = new G4MaterialPropertiesTable();
  G4double TeflonLArENE[4] = {0.1*eV, 8.0*eV, 8.3*eV, 20.0*eV};
  G4double TAREFUVP  = DSParameters::Get()->GetTeflonLArUVRef();
  G4double TAREFVISP = DSParameters::Get()->GetTeflonLArVisRef();
  G4double TeflonLArREF[4] = {TAREFVISP, TAREFVISP ,TAREFUVP , TAREFUVP };
  fLArTeflonSurfProp->AddProperty("REFLECTIVITY", TeflonLArENE, TeflonLArREF, 4);			 

  fTeflonLArSurfProp = new G4MaterialPropertiesTable();
  fTeflonLArSurfProp->AddProperty("REFLECTIVITY", TeflonLArENE, TeflonLArREF, 4);
  fOpTeflonLArSurface->SetMaterialPropertiesTable( fTeflonLArSurfProp );
}


/*
 * $Log: DSDetectorPMTDSG3.cc,v $
 * Revision 1.1  2014/07/25 14:07:19  perassos
 * Implementation of the DSG3 (TPC + NV). Configuration #8
 *
 * Revision 1.2  2014/05/08 11:00:39  pagnes
 * Added generator in PMTs stem
 *
 * Revision 1.1  2014/05/07 12:21:01  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.11  2014/03/19 16:37:28  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.10  2014/03/11 16:49:56  dfranco
 * read optical parameters from external file (DSOptics.dat)
 *
 * Revision 1.9  2013/06/22 07:21:21  dfranco
 * Fixed a bug in the photoelectron absorption in DSG3. Back to the previous QE method
 *
 * Revision 1.8  2013/06/21 13:10:03  dfranco
 * Small change in optical surfaces
 *
 * Revision 1.7  2013/06/19 18:35:28  swesterd
 * added DSScintCelll and made tpc PMTs' QE and reflections work like veto PMTs
 *
 * Revision 1.6  2013/05/31 15:40:06  dfranco
 * fixed some optical properties of the TPC
 *
 * Revision 1.5  2013/05/30 12:34:53  dfranco
 * Fixing the optical properties of the TPC. Not yet concluded
 *
 * Revision 1.4  2013/05/29 16:41:09  dfranco
 * added logger
 *
 *
 */
