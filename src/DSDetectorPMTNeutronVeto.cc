#include "DSDetectorPMTNeutronVeto.hh"
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

DSDetectorPMTNeutronVeto::DSDetectorPMTNeutronVeto(G4VPhysicalVolume *myMotherVolume) {
  
  fMotherVolume = myMotherVolume;
    
  DSLog(routine) << " Constructing PMTNeutronVeto Geometry" << endlog ;

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
  const G4ThreeVector PMT_slide(0,0,(0.5*PMT_cyl_h)*mm);
  const G4ThreeVector PC_slide(0,0,(3*PMT_height-PHOTOCATHODE_h));
  const G4ThreeVector PMTtop_slide(0,0,(2.*PMT_height)*mm);
  G4RotationMatrix* RotationPCcyl = new G4RotationMatrix;
  RotationPCcyl->rotateZ(M_PI/2.);
  G4RotationMatrix* RotationPMTtop = new G4RotationMatrix;
  RotationPMTtop->rotateY(M_PI);
 
  //////////////////////////////////////////////
  //        PMT Volumes definition            //
  //////////////////////////////////////////////

  // Cathode
  fSolidPhotocathode= new G4Paraboloid("Phothode_Solid",PHOTOCATHODE_h,0,PHOTOCATHODE_Rmax);
  fLogicPhotocathode= new G4LogicalVolume(fSolidPhotocathode,DSMaterial::Get()->GetBialkali(),"Photocathode_Logic");
  fLogicPhotocathode->SetVisAttributes(new G4VisAttributes(G4Colour(1,0,0)));

  // IMPORTANT: set here the index of the cathode material
  DSStorage::Get()->SetVetoPMTMaterialIndex(fLogicPhotocathode->GetMaterial()->GetIndex());

  // The Bulb
  fSolidPMTtop      = new G4Paraboloid("PMTtop_Solid",PMT_height,PMT_Rmin,PMT_Rmax);
  fLogicPMTtop      = new G4LogicalVolume(fSolidPMTtop,DSMaterial::Get()->GetBialkali(),"PMTtop_Logic");

  fSolidPMTbottom   = new G4Paraboloid("PMTbottom_Solid",PMT_height,PMT_Rmin,PMT_Rmax);
  fLogicPMTbottom   = new G4LogicalVolume(fSolidPMTbottom,DSMaterial::Get()->GetStainlessSteel(),"PMTbottom_Logic");

  fSolidPMTbulb     = new G4UnionSolid("PMTbulb_Solid",fSolidPMTbottom,fSolidPMTtop,RotationPMTtop,PMTtop_slide);
  //fLogicPMTbulb     = new G4LogicalVolume(fSolidPMTbulb,DSMaterial::Get()->GetStainlessSteel(),"PMTbulb_Logic");

  // Base
  fSolidPMTcyl      = new G4Tubs("PMTcyl_Solid",0,PMT_cyl_d/2.,PMT_cyl_h/2.,0.,twopi*rad);
  fLogicPMTcyl      = new G4LogicalVolume(fSolidPMTcyl,DSMaterial::Get()->GetStainlessSteel(),"PMTcyl_Logic");

  fSolidPMTbase     = new G4SubtractionSolid("PMTbase_Solid",fSolidPMTcyl,fSolidPMTbottom,0,PMT_slide);
  fLogicPMTbase     = new G4LogicalVolume(fSolidPMTbase,DSMaterial::Get()->GetStainlessSteel(),"PMTbase_Logic");
  fLogicPMTbase->SetVisAttributes(new G4VisAttributes(G4Colour(0,0,1)));

  //////////////////////////////////////////////
  //        PMT Volumes placements            //
  //////////////////////////////////////////////

  G4double RadiusPMTtop = 1845*mm;//-PHOTOCATHODE_h; 
  G4double RadiusPMTbulb = (RadiusPMTtop + 2*PMT_height);
  //  G4double RadiusPMTbulb = (RadiusPhotocathode+2*PMT_height+0.5*PHOTOCATHODE_h)*mm;
  //  G4double RadiusPMTbulb = (RadiusPhotocathode + PMT_height - PHOTOCATHODE_h/3. + PC_slide.mag());

  G4double RadiusPMTbase = (RadiusPMTtop+2*PMT_height+0.5*PMT_cyl_h)*mm;

  //  ofstream posfile;
  //  posfile.open("vpmtPos.txt");
  
  G4double theta[8];
  theta[0]=((80*twopi*rad)/360);
  theta[1]=((60*twopi*rad)/360);
  theta[2]=((40*twopi*rad)/360);
  theta[3]=((20*twopi*rad)/360);
  theta[4]=((100*twopi*rad)/360);
  theta[5]=((120*twopi*rad)/360);
  theta[6]=((140*twopi*rad)/360);
  theta[7]=((160*twopi*rad)/360);
  
  G4double start_value[8];
  start_value[0]=9;
  start_value[1]=18;
  start_value[2]=27;
  start_value[3]=22.5;
  start_value[4]=9;
  start_value[5]=6;
  start_value[6]=3;
  start_value[7]=22.5;
  
  G4double coef[8];
  coef[0]=18;
  coef[1]=24;
  coef[2]=30;
  coef[3]=45;
  coef[4]=18;
  coef[5]=24;
  coef[6]=30;
  coef[7]=45;
  
  G4int nmax[8];
  nmax[0]=20;
  nmax[1]=15;
  nmax[2]=12;
  nmax[3]=8;
  nmax[4]=20;
  nmax[5]=15;
  nmax[6]=12;
  nmax[7]=8;
  
  G4RotationMatrix* RotationPC = new G4RotationMatrix;
  RotationPC->rotateX(M_PI);

  //  int nPMT(0);
  vector<G4double> phis = DSParameters::Get()->GetVPMTTheta(); //Inconsistent angle name conventions are being used
  vector<G4double> thetas = DSParameters::Get()->GetVPMTPhi();
  vector<G4double> ids = DSParameters::Get()->GetVPMTNumber();
  for(int nPMT; nPMT < (int)DSParameters::Get()->GetNVPMTs(); nPMT++){
    DSLog(debugging) << "constructing pmt # : " << nPMT << endlog;
    G4double phi = phis[nPMT];
    G4double theta = thetas[nPMT];
    G4double id = ids[nPMT];

    G4RotationMatrix* RotationPMTbulb = new G4RotationMatrix;
    RotationPMTbulb->rotateZ(-phi);
    RotationPMTbulb->rotateY(-(theta+M_PI));
	  
  // G4RotationMatrix* RotationPMTbottom = new G4RotationMatrix;
  RotationPMTtop->rotateZ(-phi);
  RotationPMTtop->rotateY(-theta);
	  
  G4RotationMatrix* RotationPMTbase = new G4RotationMatrix;
  RotationPMTbase->rotateZ(-phi);
  RotationPMTbase->rotateY(-(theta+M_PI));
	  
  G4ThreeVector RelocationPMTbulb(RadiusPMTtop*sin(theta)*cos(phi),RadiusPMTtop*sin(theta)*sin(phi),RadiusPMTtop*cos(theta));
  G4ThreeVector RelocationPMTbottom(RadiusPMTbulb*sin(theta)*cos(phi),RadiusPMTbulb*sin(theta)*sin(phi),RadiusPMTbulb*cos(theta));	 
  G4ThreeVector RelocationPMTbase(RadiusPMTbase*sin(theta)*cos(phi),RadiusPMTbase*sin(theta)*sin(phi),RadiusPMTbase*cos(theta));
	  
  G4String number = G4UIcommand::ConvertToString(id);
  G4String namebulb="VPMTbulb_"+number;
  G4String namephotocathode="VPMT_"+number;
  G4String namebase="VPMTMuMetal_"+number;
	  
  fLogicPMTbulb     = new G4LogicalVolume(fSolidPMTbulb,DSMaterial::Get()->GetStainlessSteel(),namebulb+"_Logic");
  fLogicPMTbulb->SetVisAttributes(new G4VisAttributes(G4Colour(0,1,0)));

  // Back of the cathode
  fPhysicPMTbulb                  = new G4PVPlacement(RotationPMTbulb,
						      RelocationPMTbottom,
						      namebulb,
						      fLogicPMTbulb,
						      fMotherVolume,
						      true,
						      volumeNumber_bulb+nPMT,
						      DSStorage::Get()->GetCheckOverlap());
  
  DSLog(development) << fPhysicPMTbulb->GetName() << " = " << fPhysicPMTbulb->GetCopyNo() << endlog;
	  
  // MuMetal
  fPhysicPMTbase               = new G4PVPlacement(RotationPMTbase,
						   RelocationPMTbase,
						   namebase,
						   fLogicPMTbase,
						   fMotherVolume,
						   true,
						   volumeNumber_base+nPMT,
						   DSStorage::Get()->GetCheckOverlap());
  
  DSLog(development) << fPhysicPMTbase->GetName() << " = " << fPhysicPMTbase->GetCopyNo() << endlog;
  
  // Cathode
  fPhysicPhotocathode             = new G4PVPlacement(RotationPC,
						      PC_slide,
						      namephotocathode,
						      fLogicPhotocathode,
						      fPhysicPMTbulb,
						      true,
						      volumeNumber_photocathode+nPMT,
						      DSStorage::Get()->GetCheckOverlap());
  
  DSLog(development) << fPhysicPhotocathode->GetName() << " = " << fPhysicPhotocathode->GetCopyNo() << endlog;
  
  DSLog(debugging) << "Placing PMT " << id << endl<< endlog ;    
  DSLog(debugging) << "Theta = " << (theta*360)/(twopi*rad) << ". and Phi = " << (phi*360)/(twopi*rad) << "." << endlog; 
  /*
    posfile << nPMT 
    << "\t" << fPhysicPMTbulb->GetObjectTranslation().getX()/m
    << "\t" << fPhysicPMTbulb->GetObjectTranslation().getY()/m
    << "\t" << fPhysicPMTbulb->GetObjectTranslation().getZ()/m << std::endl;
  */
  //  posfile.close();
  /*
    G4RotationMatrix* RotationPC = new G4RotationMatrix;
    RotationPC->rotateX(M_PI);
    
    // Cathode
    fPhysicPhotocathode             = new G4PVPlacement(RotationPC,
    PC_slide,
    "VPMT_0",
    fLogicPhotocathode,
    fPhysicPMTbulb,
    true,
    0,
    DSStorage::Get()->GetCheckOverlap());
  */
  }
  DefineSurfaces();
}

DSDetectorPMTNeutronVeto::~DSDetectorPMTNeutronVeto(){
  ; //delete fMessenger;
}

void DSDetectorPMTNeutronVeto::DefineSurfaces(){

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
  /* 
  fPhotocathodeSurface = new G4LogicalBorderSurface("PhotocathodeSurface",
						    fPhysicPhotocathode,
						    fMotherVolume,
						    fOpVPhotocathodeSurface);*/
 
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
 * $Log: DSDetectorPMTNeutronVeto.cc,v $
 * Revision 1.3  2015/01/07 16:45:53  pagnes
 * changed veto optical properties format from arrays to vectors
 *
 * Revision 1.2  2014/10/13 18:43:52  swesterd
 * fixed veto PMT positions
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
