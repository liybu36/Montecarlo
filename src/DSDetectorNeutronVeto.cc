#include "DSDetectorNeutronVeto.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSDetectorPMTNeutronVeto.hh"
#include "DSMaterial.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "DSParameters.hh"
#include <iostream>
#include <fstream>
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4VisAttributes.hh"
using namespace std;

DSDetectorNeutronVeto::DSDetectorNeutronVeto(G4VPhysicalVolume *myMotherVolume) {

  fMotherVolume = myMotherVolume;

  G4int volumeNumber = 10000;

  DSLog(routine) << " Constructing NeutronVeto Geometry" << endlog ;


  G4bool cellTest = false;
  G4bool Organpipe = false;
  
  /*  
  //Water box
  fSolidNeutronVeto  = new G4Box("NeutronVeto_Solid",3*m,3*m,3*m); 
  fLogicNeutronVeto  = new G4LogicalVolume(fSolidNeutronVeto, DSMaterial::Get()->GetWater(), "NeutronVeto_Logic");
  fPhysicNeutronVeto = new G4PVPlacement(0,
				 G4ThreeVector(0,0,0),
				 "NeutronVeto",
				 fLogicNeutronVeto,
				 fMotherVolume,
				 false,
				 volumeNumber++,
				 DSStorage::Get()->GetCheckOverlap());
  DSLog(routine) << fPhysicNeutronVeto->GetName() << " = " << fPhysicNeutronVeto->GetCopyNo() << endlog;
  */
  //Stainless Steel Vessel

  fSolidSteelVessel  = new G4Orb("SteelVessel_Solid", 2008.*mm);
  fLogicSteelVessel  = new G4LogicalVolume(fSolidSteelVessel, DSMaterial::Get()->GetStainlessSteel(), "SteelVessel_Logic");
  fPhysicSteelVessel = new G4PVPlacement(0,
					 G4ThreeVector(0,0,0),
					 "SteelVessel",
					 fLogicSteelVessel,
					 fMotherVolume,
					 false,
					 ++volumeNumber,
					 DSStorage::Get()->GetCheckOverlap());
      
  DSLog(routine) << fPhysicSteelVessel->GetName() << " = " << fPhysicSteelVessel->GetCopyNo() << endlog;
    

  //Boron Scintillator
  fSolidBScintillator  = new G4Orb("BoronScintillator_Solid", 2000.*mm);
  if(!DSStorage::Get()->GetScintillator()) {
    fLogicBScintillator  = new G4LogicalVolume(fSolidBScintillator, DSMaterial::Get()->GetBoronScintillator(), "BoronScintillator_Logic");
    DSLog(routine) << "Borate Scintillator in the neutron veto" << endlog ;
  } else { 
    fLogicBScintillator  = new G4LogicalVolume(fSolidBScintillator, DSMaterial::Get()->GetGdScintillator(), "BoronScintillator_Logic");
    DSLog(routine) << "Gd Scintillator in the neutron veto" << endlog ;
  }


  if(Organpipe)
    { 
      G4double holelength = 10*cm;
      G4double holeradius = 2*cm;
      G4double hole_x = 120.48*cm;
      G4double hole_y = 0;
      G4double hole_z = 160.64*cm;
      G4double pipelength = 275.5*cm;//345.5*cm
      G4double piperadius = holeradius;
      G4double pipethickness = 2*cm;
      G4double pipe_x = 77.4*2008./2004.*cm;
      G4double pipe_y = 0*cm;
      G4double pipe_z = 184.85*2008./2004.*cm;
      G4double steelthickness = 8*mm;
      G4ThreeVector SurfaceHolePosition(hole_x,hole_y,hole_z);
      G4ThreeVector Stainlesstranslate(-hole_x,0,-hole_z+(holelength-pipelength)/2.);
      DSLog(routine)<<"hole_x= "<<hole_x<<"  "<<"pipe_z= "<<hole_z-holelength/2.+pipelength<<endlog;
           
      fSolidSurfaceHole = new G4Tubs("SurfaceHole_Solid",0, holeradius, holelength/2., 0, twopi*rad);
      fSurfaceHole= new G4IntersectionSolid("SurfaceHole",fSolidSteelVessel,fSolidSurfaceHole,0,SurfaceHolePosition);
      // fSurfaceHole= new G4IntersectionSolid("SurfaceHole",fSolidSurfaceHole,fSolidSteelVessel,0,Stainlesstranslate);
      fLogicSurfaceHole = new G4LogicalVolume(fSurfaceHole, DSMaterial::Get()->GetBoronScintillator(),"SurfaceHole_Logic");
      fPhysicSurfaceHole = new G4PVPlacement(0, G4ThreeVector(0,0,0),"SurfaceHole",fLogicSurfaceHole,fPhysicSteelVessel,false,volumeNumber++,DSStorage::Get()->GetCheckOverlap());
      fLogicSurfaceHole->SetVisAttributes(new G4VisAttributes(G4Color(0,0,1.)));
      DSLog(routine) << fPhysicSurfaceHole->GetName() << " = " << fPhysicSurfaceHole->GetCopyNo() << endlog;
      
      //      G4ThreeVector pipeBodyPosition(pipe_x, pipe_y, pipe_z);
      G4ThreeVector Pipetranslate(-pipe_x,0,-pipe_z+(holelength-pipelength)/2.);
      fSolidpipeBody = new G4Tubs("pipeBody_Solid", 0, piperadius+pipethickness, pipelength/2.,0,twopi*rad);
      fpipeBody = new G4SubtractionSolid("pipeBody",fSolidpipeBody,fSolidSteelVessel,0,Stainlesstranslate);
      fLogicpipe = new G4LogicalVolume(fpipeBody, DSMaterial::Get()->GetStainlessSteel(),"pipe_Logic");
      fPhysicpipe = new G4PVPlacement(0, -Stainlesstranslate, "pipe",fLogicpipe,fMotherVolume, false, volumeNumber++,DSStorage::Get()->GetCheckOverlap());
      fLogicpipe->SetVisAttributes(new G4VisAttributes(G4Color(0,1,0)));
      DSLog(routine) << fPhysicpipe->GetName() << " = " << fPhysicpipe->GetCopyNo() << endlog;
    
      fSolidpipeScintBody = new G4Tubs("pipeScintBody_Solid", 0, piperadius, pipelength/2. ,0,twopi*rad);
      fpipeScintBody = new G4SubtractionSolid("pipeScintBody",fSolidpipeScintBody,fSolidSteelVessel,0,Stainlesstranslate);
      fLogicpipeScint = new G4LogicalVolume(fpipeScintBody,DSMaterial::Get()->GetBoronScintillator(),"pipeScint_logic");
      fPhysicpipeScint = new G4PVPlacement(0 ,G4ThreeVector(0,0,0), "pipeScint",fLogicpipeScint,fPhysicpipe, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());
      fLogicpipeScint->SetVisAttributes(new G4VisAttributes(G4Color(1,0,0)));  
      DSLog(routine) << fPhysicpipeScint->GetName() << " = " << fPhysicpipeScint->GetCopyNo() << endlog;
  }

    fPhysicBScintillator = new G4PVPlacement(0,
					   G4ThreeVector(0,0,0),
					   "BoronScintillatorVolume",
					   fLogicBScintillator,
					   fPhysicSteelVessel,
					   false,
					   911,
					   DSStorage::Get()->GetCheckOverlap());

  DSLog(routine) << fPhysicBScintillator->GetName() << " = " << fPhysicBScintillator->GetCopyNo() << endlog;



  //PMT
  new DSDetectorPMTNeutronVeto(fPhysicBScintillator);
  /*  
  //Lumirror Sheath outside flat part of cryostat
  G4double sheathGap = 1.*nm;
  G4double cryoSheathMaxZ = DSParameters::Get()->GetCryoSheathZ()[0];
  G4double cryoSheathMinZ = DSParameters::Get()->GetCryoSheathZ()[1];
  G4double cryoSheathR = DSParameters::Get()->GetCryoSheathR()[0]+sheathGap;
  G4double sheathThickness = 1.*nm;
  fSolidCryoSheath = new G4Tubs("CryoSheath_Solid", 
				cryoSheathR, 
				cryoSheathR+sheathThickness,
				(cryoSheathMaxZ - cryoSheathMinZ)/2.,
				0, 2*M_PI);
  fLogicCryoSheath = new G4LogicalVolume(fSolidCryoSheath, DSMaterial::Get()->GetBoronScintillator(), "CryoSheath_Logic");
  fPhysicCryoSheath = new G4PVPlacement(0,
					G4ThreeVector(0,0,(cryoSheathMaxZ+cryoSheathMinZ)/2.),
					"CryoSheath",
					fLogicCryoSheath,
					fPhysicBScintillator,
					false,
					volumeNumber++,
					DSStorage::Get()->GetCheckOverlap());
  DSLog(routine) << fPhysicCryoSheath->GetName() << " = " << fPhysicCryoSheath->GetCopyNo() << endlog;

  fLogicCryoSheath->SetVisAttributes(new G4VisAttributes(G4Colour(0,.5,.5)));
  */

  //Cover the top flange with untreated SS
  G4double topFlangeRadius = 45.*cm;
  G4double topFlangeAngle = std::asin(topFlangeRadius/(200.*cm));
  G4double topFlangeThickness = 8.*mm;
  fSolidTopFlange = new G4Sphere("TopFlange_Solid",
				 200.*cm - topFlangeThickness,
				 200.*cm,
				 0,
				 2*M_PI,
				 0,
				 topFlangeAngle);
  fLogicTopFlange = new G4LogicalVolume(fSolidTopFlange,DSMaterial::Get()->GetStainlessSteel(), "TopFlange_Logic");
  fPhysicTopFlange = new G4PVPlacement(0,
				       G4ThreeVector(0,0,0),
				       "TopFlangeVolume",
				       fLogicTopFlange,
				       fPhysicBScintillator,
				       false,
				       volumeNumber++,
				       DSStorage::Get()->GetCheckOverlap());

  DSLog(routine) << fPhysicTopFlange->GetName() << " = " << fPhysicTopFlange->GetCopyNo() << endlog;

//*******************************************************work in progressss******************************************************************				       
//Source holder

  if(DSStorage::Get()->GetHolderSource()){ 
    //Values given from outside
    G4double    myHolderRadius               = (37.5*mm)+(DSStorage::Get()->GetHolderRadius()); //user values are referred to the bottom of the holder
    G4double    myHolderZ                    = DSStorage::Get()->GetHolderZ();
    G4double    myHolderPhi                  = DSStorage::Get()->GetHolderPhi(); 


    G4Tubs           *mySolidSteelHolder      = new G4Tubs("SteelHolder_Solid",0,12.5*mm,75/2.*mm,0,twopi*rad);
    G4LogicalVolume  *myLogicSteelHolder      = new G4LogicalVolume(mySolidSteelHolder, DSMaterial::Get()->GetStainlessSteel(), "SteelHolder_Logic");

    G4RotationMatrix *myHolderRotation        = new G4RotationMatrix;
    myHolderRotation->rotateY( -pi/2*rad );
    myHolderRotation->rotateX(myHolderPhi);

    G4VPhysicalVolume *myPhysicalSteelHolder  = new G4PVPlacement(myHolderRotation, //se vuoi fare una cosa fica definisci queste variabili nel .hh
								  G4ThreeVector(myHolderRadius*std::cos(myHolderPhi),myHolderRadius*std::sin(myHolderPhi),myHolderZ),  //remember to give this coordinate from outside Y=0 sempre
								  "SteelHolder",
								  myLogicSteelHolder,
								  fPhysicBScintillator,
								  false,
								  volumeNumber++,
								  DSStorage::Get()->GetCheckOverlap());

  //Inside of the source holder
    G4Tubs             *mySolidVacuumHolder    = new G4Tubs("VacuumHolder_Solid", 0,10.5*mm,61/2.*mm,0,twopi*rad);
    G4LogicalVolume    *myLogicVacuumHolder    = new G4LogicalVolume(mySolidVacuumHolder, DSMaterial::Get()->GetVacuum(), "VacuumHolder_Logic");
    G4VPhysicalVolume  *myPhysicalVacuumHolder = new G4PVPlacement(0,
				                                   G4ThreeVector(0,0,-5*mm),  
								   "VacuumHolder",
								   myLogicVacuumHolder,
								   myPhysicalSteelHolder,
								   false,
								   volumeNumber++,
								   DSStorage::Get()->GetCheckOverlap());


    G4Tubs            *mySolidSourceDisk      = new G4Tubs("SourceDisk_Solid", 0,10.5*mm,3/2.*mm,0,twopi*rad);
    G4LogicalVolume   *myLogicSourceDisk      = new G4LogicalVolume(mySolidSourceDisk, DSMaterial::Get()->GetTeflon(), "SourceDisk_Logic");
    G4VPhysicalVolume *myPhysicalSourceDisk   = new G4PVPlacement(0,
				                                  G4ThreeVector(0,0,-29*mm),  
								  "SourceDisk",
								  myLogicSourceDisk,
								  myPhysicalVacuumHolder,
								  false,
								  volumeNumber++,
								  DSStorage::Get()->GetCheckOverlap());






    G4VisAttributes* rot       = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
    G4VisAttributes* rosa      = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
    G4VisAttributes* green     = new G4VisAttributes(G4Colour(0.0,1.0,0.0));


    myLogicSteelHolder->SetVisAttributes(green);
    myLogicVacuumHolder->SetVisAttributes(rosa);
    myLogicSourceDisk->SetVisAttributes(rot);
  }
   
   
  DefineSurfaces();
}

DSDetectorNeutronVeto::~DSDetectorNeutronVeto(){
  ; //delete fMessenger;
}
void DSDetectorNeutronVeto::DefineSurfaces(){
  fOpUntreatedStainlessSteelSurface = new G4OpticalSurface("OpUntreatedStainlessSteelSurface");
  fOpUntreatedStainlessSteelSurface->SetType(dielectric_metal);
  fOpUntreatedStainlessSteelSurface->SetModel(unified);
  fOpUntreatedStainlessSteelSurface->SetFinish(groundbackpainted);
  fOpUntreatedStainlessSteelSurface->
    SetMaterialPropertiesTable(DSMaterial::Get()->GetUntreatedStainlessSteelMPT());

  fOpElectropolishedStainlessSteelSurface = new G4OpticalSurface("OpElectropolishedStainlessSteelSurface");
  fOpElectropolishedStainlessSteelSurface->SetType(dielectric_metal);
  fOpElectropolishedStainlessSteelSurface->SetModel(unified);
  fOpElectropolishedStainlessSteelSurface->SetFinish(groundbackpainted);
  fOpElectropolishedStainlessSteelSurface->
    SetMaterialPropertiesTable(DSMaterial::Get()->GetElectropolishedStainlessSteelMPT());

  fOpAluminumFoilSurface = new G4OpticalSurface("OpAluminumFoilSurface");
  fOpAluminumFoilSurface->SetType(dielectric_metal);
  fOpAluminumFoilSurface->SetModel(unified);
  fOpAluminumFoilSurface->SetFinish(groundbackpainted);
  fOpAluminumFoilSurface->
    SetMaterialPropertiesTable(DSMaterial::Get()->GetAluminumFoilMPT());

  fOpLumirrorSurface = new G4OpticalSurface("OpLumirrorSurface");
  fOpLumirrorSurface->SetType(dielectric_metal);
  fOpLumirrorSurface->SetModel(unified);
  fOpLumirrorSurface->SetFinish(groundfrontpainted);
  fOpLumirrorSurface->
    SetMaterialPropertiesTable(DSMaterial::Get()->GetLumirrorMPT());
  
  fSteelInnerSurface = new G4LogicalBorderSurface("SSteelInnerSurface",
						  fPhysicBScintillator,
						  fPhysicSteelVessel,
						  fOpLumirrorSurface);

  fSteelInnerSurfaceFlip = new G4LogicalBorderSurface("SSteelInnerSurface",
						      fPhysicSteelVessel,
						      fPhysicBScintillator,
						      fOpLumirrorSurface);
					   
  fSteelOuterSurface = new G4LogicalBorderSurface("SSteelOuterSurface",
						  fPhysicSteelVessel,
						  fMotherVolume,
						  fOpUntreatedStainlessSteelSurface);
  /*
  fCryoSheathSurface = new G4LogicalBorderSurface("CryoSheathSurface",
						  fPhysicBScintillator,
						  fPhysicCryoSheath,
						  fOpLumirrorSurface);
  */

  fTopFlangeSurface = new G4LogicalBorderSurface("TopFlangeSurface",
						 fPhysicBScintillator,
						 fPhysicTopFlange,
						 fOpLumirrorSurface);
    
  fpipeInnerSurface = new G4LogicalBorderSurface("PipeInnerSurface",fPhysicpipeScint,fPhysicpipe,fOpUntreatedStainlessSteelSurface);
  fpipeInnerSurfaceFlip = new G4LogicalBorderSurface("PipeInnerSurface",fPhysicpipe,fPhysicpipeScint,fOpUntreatedStainlessSteelSurface);
  fpipeOuterSurface = new G4LogicalBorderSurface("PipeOuterSurface",fPhysicpipe,fMotherVolume,fOpUntreatedStainlessSteelSurface);

  //fOpUntreatedStainlessSteelSurface);
}

/*
 * $Log: DSDetectorNeutronVeto.cc,v $
 * Revision 1.5  2015/01/14 16:58:36  dfranco
 * added source vial and correspondent commands by Laura and Bianca; manual updated
 *
 * Revision 1.4  2014/11/06 17:39:45  dfranco
 * added source and fixed a bug: removed gaseous nitrogen from the neutron veto
 *
 * Revision 1.3  2014/10/13 18:43:48  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.2  2014/05/07 14:27:31  dfranco
 * fixed some bugs and added GdScintillator
 *
 * Revision 1.1  2014/05/07 12:21:01  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.13  2014/04/18 16:21:42  swesterd
 * fixed an overlap in the G2 detector
 *
 * Revision 1.12  2013/08/27 04:07:01  swesterd
 * some fine tuning of the boron scintillator kB and scint yield, and some modifications to the DSG2 geometry
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
 * Revision 1.6  2013/05/14 04:20:22  swesterd
 * Fixed the veto PMT geometry
 *
 * Revision 1.5  2013/05/07 23:06:29  swesterd
 * added optical boundaries and Lumirror in the veto
 *
 * Revision 1.4  2013/05/07 16:15:45  swesterd
 * Optical processes now seem to be working in the boron-loaded scintillator
 *
 * Revision 1.3  2013/05/07 09:44:05  dfranco
 * Changed logical names of BoronScintillator volumes
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
