#include "DSDetectorPMTNeutronVetoG3.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSMaterial.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "DSParameters.hh"
#include <iostream>
#include <fstream>
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UIcommand.hh"
#include "G4MaterialPropertiesTable.hh"


using namespace std;



DSDetectorPMTNeutronVetoG3::DSDetectorPMTNeutronVetoG3( G4VPhysicalVolume *myMotherVolume ) {
  
  fMotherVolume = myMotherVolume;
    
  DSLog(routine) << " Constructing PMTNeutronVetoG3 Geometry" << endlog ;

  G4int volumeNumber_bulb         = 11000;
  G4int volumeNumber_base         = 12000;
  G4int volumeNumber_photocathode = 13000;

  //////////////////////////////////////////////
  //               dimensions                 //
  //////////////////////////////////////////////

  G4double PMT_height = 26.2231*mm; //32.75*mm;  //half the height of the PMT
  G4double PMT_cyl_h = 89.*mm;
  G4double PMT_cyl_d = 84.5*mm;  
  G4double PMT_Rmin = 0.*mm; //smallest radius of the PMT
  G4double PMT_Rmax = 101.*mm; //largest radius of the PMT 
  G4double PHOTOCATHODE_Rmax = 95.*mm; // largest radius of PC
  G4double PHOTOCATHODE_h = 23.20*mm; // full height of PC

  const G4ThreeVector PMT_slide( 0, 0, (0.5*PMT_cyl_h)*mm);
  const G4ThreeVector PC_slide( 0, 0, (3*PMT_height-PHOTOCATHODE_h));
  const G4ThreeVector PMTtop_slide( 0, 0, (2.*PMT_height)*mm);

  G4RotationMatrix* RotationPMTtop = new G4RotationMatrix;
  RotationPMTtop->rotateY(M_PI);
 

  //////////////////////////////////////////////
  //        PMT Volumes definition            //
  //////////////////////////////////////////////

  // Cathode
  fSolidPhotocathode= new G4Paraboloid( "Phothode_Solid", PHOTOCATHODE_h, 0, PHOTOCATHODE_Rmax);
  fLogicPhotocathode= new G4LogicalVolume( fSolidPhotocathode, DSMaterial::Get()->GetBialkali(), "Photocathode_Logic");
  fLogicPhotocathode->SetVisAttributes( new G4VisAttributes(G4Colour(1,0,0)) );

  // IMPORTANT: set here the index of the cathode material
  DSStorage::Get()->SetVetoPMTMaterialIndex( fLogicPhotocathode->GetMaterial()->GetIndex() );


  fSolidPMThalfbulb = new G4Paraboloid("PMTbottom_Solid", PMT_height, PMT_Rmin, PMT_Rmax);
  fSolidPMTbulb     = new G4UnionSolid("PMTbulb_Solid", fSolidPMThalfbulb, fSolidPMThalfbulb, RotationPMTtop, PMTtop_slide);

  // Base
  fSolidPMTcyl      = new G4Tubs("PMTcyl_Solid", 0, PMT_cyl_d/2., PMT_cyl_h/2., 0., twopi*rad);
  fSolidPMTbase     = new G4SubtractionSolid("PMTbase_Solid", fSolidPMTcyl, fSolidPMThalfbulb, 0, PMT_slide);
  fLogicPMTbase     = new G4LogicalVolume( fSolidPMTbase, DSMaterial::Get()->GetStainlessSteel(), "PMTbase_Logic");
  fLogicPMTbase->SetVisAttributes(new G4VisAttributes(G4Colour(0,0,1)));



  // PMT placement
  G4double mySSSphereDiameter = 736.6*cm;
  G4double myOffsetPMTSSSphere = 14*mm;
  G4double myRadiusPMTbase = mySSSphereDiameter/2. - myOffsetPMTSSSphere - PMT_cyl_h/2.;
  G4double myRadiusPMTbulb = myRadiusPMTbase - PMT_cyl_h/2.;

  G4double myTheta_min = 45*deg;  
  G4double myTheta_step = 15*deg;
  G4double myPhi_step = 30*deg;   // phi_min = 0;

  G4double myPMTnum = 0;
  G4double myNumPMTsInColumn = 8;  
  G4double myNumColumns = 12;

  G4bool myCheckOverlap = DSStorage::Get()->GetCheckOverlap();


  for( G4int col_idx = 0; col_idx < myNumColumns; col_idx++){
  
    for(G4int idx = 0; idx < myNumPMTsInColumn; idx++){

      G4RotationMatrix* myPMTRotation = new G4RotationMatrix;

      G4double myTheta = myTheta_min + idx * myTheta_step;
      G4double myPhi = col_idx * myPhi_step;

      myPMTRotation->rotateZ( -myPhi );
      myPMTRotation->rotateY( -(M_PI + myTheta) );
      
      G4ThreeVector myUnitaryVect( sin( myTheta) * cos( myPhi ), sin( myTheta ) * sin( myPhi ), cos( myTheta ) );
      G4ThreeVector myPMTbasePos    = myRadiusPMTbase * myUnitaryVect; 
      G4ThreeVector myPMTbulbPos    = myRadiusPMTbulb * myUnitaryVect;

      
      fPhysicPMTbase = new G4PVPlacement( myPMTRotation, 
                                  myPMTbasePos, 
                                  "VPMTMuMetal", 
                                  fLogicPMTbase, 
                                  fMotherVolume, 
                                  true, 
                                  volumeNumber_base + myPMTnum, 
                                  myCheckOverlap );

	    fLogicPMTbulb     = new G4LogicalVolume(fSolidPMTbulb, DSMaterial::Get()->GetStainlessSteel(), "PMTbulb_Logic");
      fPhysicPMTbulb = new G4PVPlacement( myPMTRotation,  
                                  myPMTbulbPos, 
                                  "VPMTbulb",   
                                  fLogicPMTbulb, 
                                  fMotherVolume, 
                                  true, 
                                  volumeNumber_bulb + myPMTnum, 
                                  myCheckOverlap );

      G4String myPMTcathodeName = "VPMT_" +  G4UIcommand::ConvertToString( myPMTnum );
      fPhysicPhotocathode = new G4PVPlacement( RotationPMTtop, 
                                  PC_slide, 
                                  myPMTcathodeName, 
                                  fLogicPhotocathode, 
                                  fPhysicPMTbulb, 
                                  true, 
                                  volumeNumber_photocathode + myPMTnum, 
                                  myCheckOverlap ); 

      myPMTnum++;
    }

  }

  DefineSurfaces();
}



DSDetectorPMTNeutronVetoG3::~DSDetectorPMTNeutronVetoG3(){
  ; //delete fMessenger;
}



void DSDetectorPMTNeutronVetoG3::DefineSurfaces(){

  fOpVPhotocathodeSurface = new G4OpticalSurface("OpVPhotocathodeSurface");
  fOpVPhotocathodeSurface->SetType(dielectric_metal);
  fOpVPhotocathodeSurface->SetModel(glisur);
  fOpVPhotocathodeSurface->SetFinish(polished);
  fOpVPhotocathodeSurface->
    SetMaterialPropertiesTable(DSMaterial::Get()->GetVPhotocathodeMPT());
  
  fOpPMTBackSurface = new G4OpticalSurface("OpPMTBackSurface");
  fOpPMTBackSurface->SetType(dielectric_metal);
  fOpPMTBackSurface->SetModel(unified);
  fOpPMTBackSurface->SetFinish(polished);
  fOpPMTBackSurface->
    SetMaterialPropertiesTable(DSMaterial::Get()->GetPMTBackMPT());

  fOpLumirrorSurface = new G4OpticalSurface("OpLumirrorSurface");
  fOpLumirrorSurface->SetType(dielectric_metal);
  fOpLumirrorSurface->SetModel(unified);
  fOpLumirrorSurface->SetFinish(groundfrontpainted);
  fOpLumirrorSurface->
    SetMaterialPropertiesTable(DSMaterial::Get()->GetLumirrorMPT());
 
  fPhotocathodeSurface = new G4LogicalBorderSurface("PhotocathodeSurface",
						    fMotherVolume,
						    fPhysicPhotocathode,
						    fOpVPhotocathodeSurface);
  fPMTbaseSurface = new G4LogicalBorderSurface("PMTbaseSurface",
					       fPhysicPMTbase,
					       fMotherVolume,
					       fOpLumirrorSurface);
  fPMTbulbSurface = new G4LogicalBorderSurface("PMTbulbSurface",
					       fPhysicPMTbulb,
					       fMotherVolume,
					       fOpPMTBackSurface);
}

/*
 * $Log: DSDetectorPMTNeutronVetoG3.cc,v $
 * Revision 1.1  2014/07/25 14:07:20  perassos
 * Implementation of the DSG3 (TPC + NV). Configuration #8
 *
 * Revision 1.1  2014/05/07 12:21:01  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.13  2013/10/15 14:55:07  dfranco
 * fixed a buf in rdmchain
 *
 * Revision 1.12  2013/08/20 03:25:52  swesterd
 * added G2 TPC geoemtry (not complete) and added monoenergetic energy distribution to generator
 *
 * Revision 1.11  2013/08/05 03:13:52  swesterd
 * some fine tuning of bscint and veto parameters
 *
 * Revision 1.10  2013/06/06 23:00:19  swesterd
 * moved veto PMT numbers from copy number to name and made a separate logical volume for each veto PMT. Assigned each physical volume in the veto a unique copy number 1wxyz, where w=1 for PMT bulbs, w=2 for PMT bases, w=3 for photocathodes, and w=0 for everything else
 *
 * Revision 1.9  2013/06/05 23:03:32  swesterd
 * moved optical boundary MPTs to DSMaterial and gave the trunks optical boundary properties consistent with untreated stainless steel
 *
 * Revision 1.8  2013/06/04 01:02:29  swesterd
 * other than the optical boundary of the trunks, the veto optics appear to be complete and up and running...modulo whatever I may have missed...
 *
 * Revision 1.7  2013/05/27 23:59:02  swesterd
 * added a (currently commented out) Lumirror sheath to the cryostat and introduced DSOpBoundaryProcess to try to figure out why the boundaries are being screwy, with some edits so that it can handle constant and vector properties with freaking out
 *
 * Revision 1.6  2013/05/25 07:58:23  swesterd
 * Got the veto PMT optical boundaries all working along with photocathode optical properties, added PMT quantum efficiency to DSTrackingAction, and added a function to DSTrackingAction that locates and quadratically interpolates points in data, for getting useful QEs
 *
 * Revision 1.5  2013/05/20 13:58:41  swesterd
 * updated the veto PMT geometry
 *
 * Revision 1.4  2013/05/14 04:20:22  swesterd
 * Fixed the veto PMT geometry
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
