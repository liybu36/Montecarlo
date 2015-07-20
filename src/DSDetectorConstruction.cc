#include <iostream>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <fstream>
#include "G4String.hh"

#include "G4RunManager.hh" 
#include "G4SDManager.hh"


#include "G4CSGSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
//migration : canceled out BREPPolycone and added GenerciPolycone
//#include  "G4GenericPolycone.hh"
#include "G4Polycone.hh"

#include "G4LogicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"


#include "G4RotationMatrix.hh"
#include "G4Point3D.hh"
#include "G4Transform3D.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"


#include "DSDetectorConstruction.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSDetectorConstructionMessenger.hh"

#include "DSDetectorNeutronVeto.hh"
#include "DSDetectorNeutronVetoG3.hh"
#include "DSDetectorWaterTank.hh"
#include "DSDetectorDS10.hh"
#include "DSDetectorDS50.hh"
#include "DSDetectorDSG2.hh"
#include "DSDetectorDSG3.hh"
#include "DSDetectorDS5k.hh"
#include "DSDetectorDS20k.hh"
#include "DSDetectorTester.hh"
#include "DSDetectorCalibrationDevice.hh"
#include "DSScintCell.hh"
#include "DSDetectorTubeTest.hh"
#include "DSDetectorSourceHolder.hh"
#include "DSStorage.hh"
#include "DSDetectorSource.hh"
#include "DSEventHandler.hh"

#include <math.h>
#include <time.h>
#include<sstream>

using namespace std;

extern G4int thePMTconfig;

G4int	NumberOfVetoPMT;

DSDetectorConstruction::DSDetectorConstruction() {
  
  time_t rawtime;
  time ( &rawtime );
  DSLog(routine) << ctime(&rawtime) << endlog ;
  
  fDetectorConfiguration = 0 ;
  
  fIsSource  = false ;
  
  fMessenger = new DSDetectorConstructionMessenger(this);
  
}

DSDetectorConstruction::~DSDetectorConstruction(){
  delete fMessenger;
}

G4VPhysicalVolume* DSDetectorConstruction::Construct() {
  return ConstructDetector();
}

void DSDetectorConstruction::DefineSizes() {

 ;

}


 

G4VPhysicalVolume* DSDetectorConstruction::ConstructDetector() { 

    // World 
  fSolidWorld  = new G4Box("World_Solid",20*m,20*m,20*m);
 
  fLogicWorld  = new G4LogicalVolume(fSolidWorld, DSMaterial::Get()->GetAir(), "World_Logic");
 
  fPhysicWorld = new G4PVPlacement(0,
                                 G4ThreeVector(0,0,0),
                                 "World",
                                 fLogicWorld,
                                 NULL,
                                 false,
                                 0);

  DSEventHandler::Get()->SetDetectorFlag(fDetectorConfiguration);

  if(fDetectorConfiguration == 0) {  // entire detector
    DSLog(routine) << " Detector Configuration - TPC+NV+WT: " <<  fDetectorConfiguration << endlog ;
    //Water Tank	 
    DSDetectorWaterTank* WTVolume   = new DSDetectorWaterTank(fPhysicWorld);
    //Neutron Veto
    DSDetectorNeutronVeto* NVVolume = new DSDetectorNeutronVeto(WTVolume->GetDetectorComponent());  
    //TPC
    new DSDetectorDS50( NVVolume->GetDetectorComponent() );
    
	  //    DSStorage::Get()->SetLArGArBoundaryPosZ( LArGArBoundaryPosZDS50  );    

    if(fIsSource)  new DSDetectorSource( NVVolume->GetDetectorComponent() ) ;


  } else if(fDetectorConfiguration == 2) { // only TPC
    DSLog(routine) << " Detector Configuration - TPC: " <<  fDetectorConfiguration << endlog ;
    //TPC
    new DSDetectorDS50( fPhysicWorld );    


  } else if(fDetectorConfiguration == 1) { // TPC + Veto
    DSLog(routine) << " Detector Configuration - TPC+NV: " <<  fDetectorConfiguration << endlog ;
    //Neutron Veto
    DSDetectorNeutronVeto* NVVolume = new DSDetectorNeutronVeto( fPhysicWorld );
    //TPC
    new DSDetectorDS50( NVVolume->GetDetectorComponent() );


  } else if(fDetectorConfiguration == 3) { // Water Tank + Veto
    DSLog(routine) << " Detector Configuration - NV+WT: " <<  fDetectorConfiguration << endlog ;
    //Water Tank
    DSDetectorWaterTank* WTVolume   = new DSDetectorWaterTank(fPhysicWorld);
    //Neutron Veto
    new DSDetectorNeutronVeto(WTVolume->GetDetectorComponent());


  } else if(fDetectorConfiguration == 4) { // Optical Tester
    DSLog(routine) << " Detector Configuration - Tester: " <<  fDetectorConfiguration << endlog ;
    //Tester
    new DSDetectorTester(fPhysicWorld);


  } else if(fDetectorConfiguration == 5) { // DS10
    DSLog(routine) << " Detector Configuration - TPC: " << fDetectorConfiguration << endlog;
    //TPC
    new DSDetectorDS10( fPhysicWorld );


  } else if(fDetectorConfiguration == 6) { // 3" Scintillator Cell
    DSLog(routine) << " Detector Configuration - 3 inch scintillator cell: " << fDetectorConfiguration << endlog;
    //3" Scintillator Cell
    new DSScintCell(fPhysicWorld);


  } else if(fDetectorConfiguration == 7) {// WT + Veto + DSG2 TPC
    DSLog(routine) << " Detector Configuration - DSG2 TPC + NV + WT: " <<  fDetectorConfiguration << endlog ;
    //Water Tank	 
    DSDetectorWaterTank* WTVolume   = new DSDetectorWaterTank(fPhysicWorld);
    //Neutron Veto
    DSDetectorNeutronVeto* NVVolume = new DSDetectorNeutronVeto(WTVolume->GetDetectorComponent());  
    //DSG2 TPC
    new DSDetectorDSG2( NVVolume->GetDetectorComponent() );
  
  } else if(fDetectorConfiguration == 8) {// G3: Veto + TPC
    DSLog(routine) << " Detector Configuration - DSG3 TPC + NV: " <<  fDetectorConfiguration << endlog ;
    //Neutron Veto
    DSDetectorNeutronVetoG3* NVVolume = new DSDetectorNeutronVetoG3(fPhysicWorld);  
    //DSG2 TPC
    new DSDetectorDSG3( NVVolume->GetDetectorComponent() );


  }else  if(fDetectorConfiguration == 9) {  // DS50TPC + NV + TubeTest +WT
    DSLog(routine) << " Detector Configuration - DS50TPC+NV+TubeTest+WT: " <<  fDetectorConfiguration << endlog ;
    //Water Tank	 
    DSDetectorWaterTank* WTVolume   = new DSDetectorWaterTank(fPhysicWorld);
    //Neutron Veto
    DSDetectorNeutronVeto* NVVolume = new DSDetectorNeutronVeto(WTVolume->GetDetectorComponent());  
    //TPC
    new DSDetectorDS50( NVVolume->GetDetectorComponent() );
    //Tube Test
    new DSDetectorTubeTest( NVVolume->GetDetectorComponent() );

  }else  if(fDetectorConfiguration == 10) {  // DS50TPC + NV + SourceHolder +WT
    DSLog(routine) << " Detector Configuration - DS50TPC+NV+SourceHolder+WT: " <<  fDetectorConfiguration << endlog ;
    //Water Tank	 
    DSDetectorWaterTank* WTVolume   = new DSDetectorWaterTank(fPhysicWorld);
    //Neutron Veto
    DSDetectorNeutronVeto* NVVolume = new DSDetectorNeutronVeto(WTVolume->GetDetectorComponent());  
    //TPC
    new DSDetectorDS50( NVVolume->GetDetectorComponent() );
    //Source Holder
    new DSDetectorSourceHolder( NVVolume->GetDetectorComponent() );


  } else if(fDetectorConfiguration == 11) {// 5k Veto + TPC
    DSLog(routine) << " Detector Configuration - DS5k: " <<  fDetectorConfiguration << endlog ;
    DSStorage::Get()->Set5KGeometry(true) ; 
    new DSDetectorDS5k( fPhysicWorld);


  }  else if(fDetectorConfiguration == 12) {// 20k - with GdWater Veto 
    DSLog(routine) << " Detector Configuration - DS20k: " <<  fDetectorConfiguration << endlog ;
    DSStorage::Get()->Set20KGeometry(true) ; 
    DSDetectorNeutronVetoG3* NVVolume = new DSDetectorNeutronVetoG3(fPhysicWorld);  
    new DSDetectorDS20k( NVVolume->GetDetectorComponent());


  }else if(fDetectorConfiguration == 808) {// WT + NV + DS50 TPC + Collimator - 808 to leave a gap with respect to "real" detector configurations. This is for testing purposes
    DSLog(routine) << " Detector Configuration - TPC+NV+WT+Collimator: " << fDetectorConfiguration << endlog;
    //Water Tank	 
    DSDetectorWaterTank* WTVolume   = new DSDetectorWaterTank(fPhysicWorld);
    //Neutron Veto
    DSDetectorNeutronVeto* NVVolume = new DSDetectorNeutronVeto(WTVolume->GetDetectorComponent());  
    //TPC
    new DSDetectorDS50( NVVolume->GetDetectorComponent() );
    //Collimator
    new DSDetectorCalibrationDevice( NVVolume->GetDetectorComponent() );

  }  

  DSLog(development) << " LAr - GAr boundary z coordinate set to: " << DSStorage::Get()->GetLArGArBoundaryPosZ()/cm << " cm"  << endlog ;

  return fPhysicWorld;
}


void DSDetectorConstruction::PrintDetectorParameters() {
}


void DSDetectorConstruction::UpdateGeometry() {
   G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());
}



/*
 * $Log: DSDetectorConstruction.cc,v $
 * Revision 1.10  2015/04/28 10:16:21  pagnes
 * G3 veto added to DS20k geometry
 *
 * Revision 1.9  2015/04/23 14:04:01  pagnes
 * DS20K geometry added (config 10)
 *
 * Revision 1.8  2015/03/09 15:20:37  pagnes
 * DS 5tons geometry added (conf 9)
 *
 * Revision 1.7  2015/01/14 16:58:35  dfranco
 * added source vial and correspondent commands by Laura and Bianca; manual updated
 *
 * Revision 1.6  2014/11/06 17:39:44  dfranco
 * added source and fixed a bug: removed gaseous nitrogen from the neutron veto
 *
 * Revision 1.5  2014/07/25 14:10:09  perassos
 * Improved handling of the LArGArBoundaryZ variable
 *
 * Revision 1.4  2014/07/25 08:41:40  reinhold
 * minor changes in collimator, better message printed to logfile. From an interface perspective: in order to have the collimator call '/ds/detector/configuration 808' in your macro (instead of '/ds/detector/configuration 8')
 *
 * Revision 1.3  2014/07/25 06:27:58  eedkins
 * Added DSDetectorCalibrationDevice class, featuring a collimator but no additional insertion system features. Parameters are hard coded - macro commands will be added later. added temporarily fDetectorConfiguration 8 for testing purposes.
 *
 * Revision 1.2  2014/07/23 14:55:34  pagnes
 * S2 first tuning
 *
 * Revision 1.1  2014/05/07 12:21:01  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.11  2013/11/19 10:33:20  perassos
 * Added methods to handle the electric field and the liquid/gas interface z coordinate
 *
 * Revision 1.10  2013/08/20 03:25:52  swesterd
 * added G2 TPC geoemtry (not complete) and added monoenergetic energy distribution to generator
 *
 * Revision 1.9  2013/08/06 13:58:20  dfranco
 * Added 3 variables to store the deposited energy in LAr, scintillator, and water. The last two are not yet implemented. g4rooter has been updated with 3 new variables: tpcene, vetoene, and muene
 *
 * Revision 1.8  2013/08/05 03:13:52  swesterd
 * some fine tuning of bscint and veto parameters
 *
 * Revision 1.7  2013/06/19 08:41:40  dfranco
 * removed DSCintScell.hh from DSDectectorConstruction. This file does not exist
 *
 * Revision 1.6  2013/06/17 18:23:16  swesterd
 * added some things for testing scintillator properties, including an incompelte ScintCell geometry
 *
 * Revision 1.5  2013/06/10 11:55:01  perassos
 * DS10 TPC added
 *
 * Revision 1.4  2013/05/31 13:02:15  dfranco
 * Added a detector tester, with simpplified geometry (configuration number = 4) to test optical properties of the materials
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
