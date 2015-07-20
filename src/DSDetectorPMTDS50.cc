#include <iostream>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4PVPlacement.hh"
#include "G4UIcommand.hh"

#include "DSDetectorPMTDS50.hh"
#include "DSStorage.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"




using namespace std;


DSDetectorPMTDS50::DSDetectorPMTDS50( G4VPhysicalVolume* myMotherVolume, G4VPhysicalVolume* myTeflonCapDiskVolume ){

  fMotherVolume        = myMotherVolume;
  fTeflonCapDiskVolume = myTeflonCapDiskVolume;


  const double myTwoPi    = 2*M_PI*rad;
  G4bool   myCheckOverlap = DSStorage::Get()->GetCheckOverlap();



  // Variables to distinguish bewteen top and bottom
  G4int  myOffset = 19;
  if( fMotherVolume->GetRotation() == 0 ) {
    myOffset = 0;
    DSLog(routine) << " Constructing Top TPC PMTs Geometry  " << endlog;
  }
  else DSLog(routine) << " Constructing Bottom TPC PMTs Geometry  " << endlog;

  // Teflon Cap
  G4double myPMTAssembly_h = DSParameters::Get()->GetPMTAssemblyHeight()/mm;
  G4double myTeflonCap_h = DSParameters::Get()->GetTeflonCapHeight()/mm;
  G4double myPMTBody_d = DSParameters::Get()->GetPMTBodyDiameter()/mm;
  G4double myPMTBody_h = DSParameters::Get()->GetPMTBodyHeight()/mm;
  G4double myPMTHead_d = DSParameters::Get()->GetPMTHeadDiameter()/mm;
  G4double myPMTHead_h = DSParameters::Get()->GetPMTHeadHeight()/mm;
  G4double myPMTWall   = DSParameters::Get()->GetPMTWallThickness()/mm; 
  G4double myPMTWindowThickness = DSParameters::Get()->GetPMTWindowThickness()/mm;
  G4double myPMTOffset = DSParameters::Get()->GetPMTOffset()/mm;
 
  DSLog(routine) << "myPMTAssembly_h= " <<myPMTAssembly_h << endlog;
  DSLog(routine) << "myTeflonCap_h= " <<myTeflonCap_h << endlog;
  DSLog(routine) << "myPMTBody_d= " <<myPMTBody_d << endlog;
  DSLog(routine) << "myPMTBody_h= " <<myPMTBody_h << endlog;
  DSLog(routine) << "myPMTHead_d= " <<myPMTHead_h << endlog;
  DSLog(routine) << "myPMTHead_h= " <<myPMTHead_h << endlog;
  DSLog(routine) << "myPMTWall= " <<myPMTWall << endlog;
  DSLog(routine) << "myPMTWindowsThickness= " <<myPMTWindowThickness << endlog;
  DSLog(routine) << "myPMTOffset= " <<myPMTOffset << endlog;  

/*
  cout << myPMTAssembly_h << endl ;
  cout << myTeflonCap_h<< endl ;
  cout << myPMTBody_d/2. << endl ;
  cout << myPMTBody_h << endl ;
  cout << myPMTHead_d/2. << endl ;
  cout << myPMTHead_h << endl ;
  cout << myPMTWall   << endl ;
  cout << myPMTWindowThickness<< endl ;
  cout << myPMTOffset << endl ;
  */
  // PMT Surrounding
  fSolidPMTSurrounding = new G4Tubs( "PMTSurrounding_Solid", 
                                      0,  
                                      myPMTHead_d/2.,	     
                                      myTeflonCap_h/2.,      					
				      0, 
				      myTwoPi );
  // Part of the PMT inside the surrounding (Head, vacuum, window and body) 
  fSolidPMTHead        = new G4Tubs( "PMTHead_Solid",        
                                      0,  
				      myPMTHead_d/2.,              
				      myPMTHead_h/2.,                                                    
				      0, 
				      myTwoPi );
				      
  fSolidPMTHeadVac     = new G4Tubs( "PMTHeadVac_Solid",     
                                      0,  
				      myPMTHead_d/2. - myPMTWall,  myPMTHead_h/2. - myPMTWall/2.,                                     
				      0, 
				      myTwoPi );
				      
  fSolidPMTWindow      = new G4Tubs( "PMTWindow_Solid",      
                                      0,  
				      myPMTHead_d/2. - myPMTWall,  myPMTWindowThickness/2.,                                           
				      0, 
				      myTwoPi );
				      
  fSolidPMTJoinVac     = new G4Tubs( "PMTJoinVac_Solid",     
                                      0,  
				      myPMTBody_d/2. - myPMTWall,  
				      myPMTWall/2.,                                                      
				      0, 
				      myTwoPi );
				      
  fSolidPMTBody        = new G4Tubs( "PMTBody_Solid",        
                                      0, 
				      myPMTBody_d/2.,              
				      myTeflonCap_h/2. - myPMTHead_h/2. - myPMTOffset/2.,                
				      0, 
				      myTwoPi ); 
				      
  fSolidPMTBodyVac     = new G4Tubs( "PMTBodyVac_Solid",     
                                      0,  
				      myPMTBody_d/2. - myPMTWall,  
				      myTeflonCap_h/2. - myPMTHead_h/2. - myPMTOffset/2.,                
				      0, 
				      myTwoPi );
				      
  // Part of the PMT body outside the surrounding
  fSolidPMTTop         = new G4Tubs( "PMTTop_Solid",         
                                      0,  
				      myPMTBody_d/2.,              
				      myPMTBody_h/2. - myTeflonCap_h/2. + myPMTOffset/2.,                
				      0, 
				      myTwoPi );
				      
  fSolidPMTTopVac      = new G4Tubs( "PMTTopVac_Solid",      
                                      0,  
				      myPMTBody_d/2. - myPMTWall,  
				      myPMTBody_h/2. - myTeflonCap_h/2. + myPMTOffset/2. - myPMTWall/2., 
				      0, 
				      myTwoPi );
				      
  fSolidPMTLArDisk     = new G4Tubs( "PMTLArDisk_Solid",     
                                      0,  
				      myPMTHead_d/2.,              
				      myPMTOffset/2.,                                                    
				      0, 
				      myTwoPi );


  fSolidPMTStem        = new G4Tubs( "PMTStem_Solid",
                                      0,
              myPMTBody_d/2. - myPMTWall,
              myPMTWall/2.,
              0,
              myTwoPi );            


  // Positioning
  G4double myPMTSpacing = DSParameters::Get()->GetPMTSpacing()/mm;
  G4ThreeVector myPMTTopZ( 0, 0, -myPMTAssembly_h/2. + myPMTBody_h/2. + myPMTOffset/2. + myTeflonCap_h/2. );
  G4ThreeVector myPMTPos[19];
  G4ThreeVector myPMTRad1( myPMTSpacing, 0, 0 );
  G4ThreeVector myPMTRad2 = myPMTRad1 + G4ThreeVector( myPMTSpacing, 0, 0).rotateZ( myTwoPi/6. );
  G4ThreeVector myPMTRad3 = 2*myPMTRad1;
  
  
  // 3 circles of PMTs at radii myPMTRad1, myPMTRad2, myPMTRad3
  myPMTPos[0].set( 0, 0, 0);
  myPMTPos[1]  = myPMTRad1;
  myPMTPos[2]  = myPMTRad1.rotateZ( myTwoPi/6. );
  myPMTPos[3]  = myPMTRad1.rotateZ( myTwoPi/6. );
  myPMTPos[4]  = myPMTRad1.rotateZ( myTwoPi/6. );
  myPMTPos[5]  = myPMTRad1.rotateZ( myTwoPi/6. );
  myPMTPos[6]  = myPMTRad1.rotateZ( myTwoPi/6. );
  myPMTPos[7]  = myPMTRad2;
  myPMTPos[8]  = myPMTRad2.rotateZ( myTwoPi/6. );
  myPMTPos[9]  = myPMTRad2.rotateZ( myTwoPi/6. );
  myPMTPos[10] = myPMTRad2.rotateZ( myTwoPi/6. );
  myPMTPos[11] = myPMTRad2.rotateZ( myTwoPi/6. );
  myPMTPos[12] = myPMTRad2.rotateZ( myTwoPi/6. );
  myPMTPos[13] = myPMTRad3;
  myPMTPos[14] = myPMTRad3.rotateZ( myTwoPi/6. );
  myPMTPos[15] = myPMTRad3.rotateZ( myTwoPi/6. );
  myPMTPos[16] = myPMTRad3.rotateZ( myTwoPi/6. );
  myPMTPos[17] = myPMTRad3.rotateZ( myTwoPi/6. );
  myPMTPos[18] = myPMTRad3.rotateZ( myTwoPi/6. );

  // shifts for assemblying the pmt
  G4ThreeVector myZeros        ( 0, 0, 0 ); 
  G4ThreeVector myPMTHeadPos   ( 0, 0, -myTeflonCap_h/2. + myPMTHead_h/2. + myPMTOffset );
  G4ThreeVector myPMTHeadVacPos( 0, 0, -myPMTWall/2. );
  G4ThreeVector myPMTStemPos   ( 0, 0, fSolidPMTTop->GetDz() - myPMTWall/2.);
  G4ThreeVector myPMTWindowPos ( 0, 0, -myPMTHead_h/2. + myPMTWall/2. + myPMTWindowThickness/2. );
  G4ThreeVector myPMTJoinVacPos( 0, 0,  myPMTHead_h/2. - myPMTWall/2. );
  G4ThreeVector myPMTBodyPos   ( 0, 0,  myPMTHead_h/2. + myPMTOffset/2. );
  G4ThreeVector myPMTTopPos[19];
  G4ThreeVector myPMTLArDiskPos( 0, 0, -myTeflonCap_h/2. + myPMTOffset/2. );
 
  DSLog(routine)<<"myPMTHeadPos= "<< -myTeflonCap_h/2. + myPMTHead_h/2. + myPMTOffset <<endlog;
  DSLog(routine)<<"myPMTWindowsPos= "<< -myPMTHead_h/2. + myPMTWall/2. + myPMTWindowThickness/2. <<endlog;
 
 
  /*cout << myPMTHeadPos <<  " "
  << myPMTHeadVacPos <<  " "
  << myPMTWindowPos <<  " "
  << myPMTJoinVacPos <<  " "
  << myPMTBodyPos <<  " "
  << myPMTLArDiskPos <<  " " << endl ;
  */
   
  G4String myPMTSurLogicName[19];
  G4String myPMTHeadLogicName[19];
  G4String myPMTHeadVacLogicName[19];
  G4String myPMTWindowLogicName[19];
  G4String myPMTJoinVacLogicName[19];
  G4String myPMTBodyLogicName[19];
  G4String myPMTBodyVacLogicName[19];
  G4String myPMTTopLogicName[19];
  G4String myPMTTopVacLogicName[19];
  G4String myPMTLArDiskLogicName[19];
  G4String myPMTStemLogicName[19];

  G4String myPMTSurName[19];
  G4String myPMTHeadName[19];
  G4String myPMTHeadVacName[19];
  G4String myPMTWindowName[19];
  G4String myPMTJoinVacName[19];
  G4String myPMTBodyName[19];
  G4String myPMTBodyVacName[19];
  G4String myPMTTopName[19];
  G4String myPMTTopVacName[19];
  G4String myPMTLArDiskName[19];
  G4String myPMTStemName[19];


  for( int i = 0; i < 19; i++) {

    //if( i == 1 ) myCheckOverlap = false;

    myPMTSurLogicName[i]     = "PMTSurrounding_" + G4UIcommand::ConvertToString( i + myOffset ) + "_Logic";
    myPMTHeadLogicName[i]    = "PMTHead_"        + G4UIcommand::ConvertToString( i + myOffset ) + "_Logic";
    myPMTHeadVacLogicName[i] = "PMTHeadVac_"     + G4UIcommand::ConvertToString( i + myOffset ) + "_Logic";
    myPMTWindowLogicName[i]  = "PMTWindow_"      + G4UIcommand::ConvertToString( i + myOffset ) + "_Logic";
    myPMTJoinVacLogicName[i] = "PMTJoinVac_"     + G4UIcommand::ConvertToString( i + myOffset ) + "_Logic";
    myPMTBodyLogicName[i]    = "PMTBody_"        + G4UIcommand::ConvertToString( i + myOffset ) + "_Logic";
    myPMTBodyVacLogicName[i] = "PMTBodyVac_"     + G4UIcommand::ConvertToString( i + myOffset ) + "_Logic";
    myPMTTopLogicName[i]     = "PMTTop_"         + G4UIcommand::ConvertToString( i + myOffset ) + "_Logic";
    myPMTTopVacLogicName[i]  = "PMTTopVac_"      + G4UIcommand::ConvertToString( i + myOffset ) + "_Logic";
    myPMTLArDiskLogicName[i] = "PMTLArDisk_"     + G4UIcommand::ConvertToString( i + myOffset ) + "_Logic";
    myPMTStemLogicName[i]    = "PMTStem_"        + G4UIcommand::ConvertToString( i + myOffset ) + "_Logic";

    myPMTSurName[i]     = "PMTSurrounding_" + G4UIcommand::ConvertToString( i + myOffset );
    myPMTHeadName[i]    = "PMTHead_"        + G4UIcommand::ConvertToString( i + myOffset );
    myPMTHeadVacName[i] = "PMTHeadVac_"     + G4UIcommand::ConvertToString( i + myOffset ); 
    myPMTWindowName[i]  = "TPMT_"           + G4UIcommand::ConvertToString( i + myOffset );
    myPMTJoinVacName[i] = "PMTJoinVac_"     + G4UIcommand::ConvertToString( i + myOffset );
    myPMTBodyName[i]    = "PMTBody_"        + G4UIcommand::ConvertToString( i + myOffset );
    myPMTBodyVacName[i] = "PMTBodyVac_"     + G4UIcommand::ConvertToString( i + myOffset );
    myPMTTopName[i]     = "PMTTop_"         + G4UIcommand::ConvertToString( i + myOffset );
    myPMTTopVacName[i]  = "PMTTopVac_"      + G4UIcommand::ConvertToString( i + myOffset );
    myPMTLArDiskName[i] = "PMTLArDisk_"     + G4UIcommand::ConvertToString( i + myOffset );
    myPMTStemName[i]    = "PMTStem_"        + G4UIcommand::ConvertToString( i + myOffset );

    myPMTTopPos[i] = myPMTPos[i] + myPMTTopZ;


    fLogicPMTSurrounding[i] = new G4LogicalVolume( fSolidPMTSurrounding, DSMaterial::Get()->GetTeflon(),      myPMTSurLogicName[i] );
    fLogicPMTHead[i]        = new G4LogicalVolume( fSolidPMTHead,        DSMaterial::Get()->GetKovar(),       myPMTHeadLogicName[i] ); 
    fLogicPMTHeadVac[i]     = new G4LogicalVolume( fSolidPMTHeadVac,     DSMaterial::Get()->GetVacuum(),      myPMTHeadVacLogicName[i] );
    fLogicPMTWindow[i]      = new G4LogicalVolume( fSolidPMTWindow,      DSMaterial::Get()->GetBialkali(),    myPMTWindowLogicName[i] );
    fLogicPMTJoinVac[i]     = new G4LogicalVolume( fSolidPMTJoinVac,     DSMaterial::Get()->GetVacuum(),      myPMTJoinVacLogicName[i] );
    fLogicPMTBody[i]        = new G4LogicalVolume( fSolidPMTBody,        DSMaterial::Get()->GetKovar(),       myPMTBodyLogicName[i] );
    fLogicPMTBodyVac[i]     = new G4LogicalVolume( fSolidPMTBodyVac,     DSMaterial::Get()->GetVacuum(),      myPMTBodyVacLogicName[i] );
    fLogicPMTTop[i]         = new G4LogicalVolume( fSolidPMTTop,         DSMaterial::Get()->GetKovar(),       myPMTTopLogicName[i] );
    fLogicPMTTopVac[i]      = new G4LogicalVolume( fSolidPMTTopVac,      DSMaterial::Get()->GetVacuum(),      myPMTTopVacLogicName[i] ); 
    fLogicPMTLArDisk[i]     = new G4LogicalVolume( fSolidPMTLArDisk,     DSMaterial::Get()->GetLiquidArgon(), myPMTLArDiskLogicName[i] );  
    fLogicPMTStem[i]        = new G4LogicalVolume( fSolidPMTStem,        DSMaterial::Get()->GetFusedSilica(), myPMTStemLogicName[i] );

    // IMPORTANT: set here the index of the cathode material
    DSStorage::Get()->SetPMTMaterialIndex(fLogicPMTWindow[i]->GetMaterial()->GetIndex());

    fPhysicPMTSurrounding[i] = new G4PVPlacement( 0, myPMTPos[i],      myPMTSurName[i],     fLogicPMTSurrounding[i], fTeflonCapDiskVolume,     false, 0, myCheckOverlap ); 
    fPhysicPMTHead[i]        = new G4PVPlacement( 0, myPMTHeadPos,     myPMTHeadName[i],    fLogicPMTHead[i],        fPhysicPMTSurrounding[i], false, 0, myCheckOverlap );
    fPhysicPMTHeadVac[i]     = new G4PVPlacement( 0, myPMTHeadVacPos,  myPMTHeadVacName[i], fLogicPMTHeadVac[i],     fPhysicPMTHead[i],        false, 0, myCheckOverlap );
    fPhysicPMTWindow[i]      = new G4PVPlacement( 0, myPMTWindowPos,   myPMTWindowName[i],  fLogicPMTWindow[i],      fPhysicPMTHeadVac[i],     false, 0, myCheckOverlap );
    fPhysicPMTJoinVac[i]     = new G4PVPlacement( 0, myPMTJoinVacPos,  myPMTJoinVacName[i], fLogicPMTJoinVac[i],     fPhysicPMTHead[i],        false, 0, myCheckOverlap );
    fPhysicPMTBody[i]        = new G4PVPlacement( 0, myPMTBodyPos,     myPMTBodyName[i],    fLogicPMTBody[i],        fPhysicPMTSurrounding[i], false, 0, myCheckOverlap );
    fPhysicPMTBodyVac[i]     = new G4PVPlacement( 0, myZeros,          myPMTBodyVacName[i], fLogicPMTBodyVac[i],     fPhysicPMTBody[i],        false, 0, myCheckOverlap );
    fPhysicPMTTop[i]         = new G4PVPlacement( 0, myPMTTopPos[i],   myPMTTopName[i],     fLogicPMTTop[i],         fMotherVolume,            false, 0, myCheckOverlap );
    fPhysicPMTTopVac[i]      = new G4PVPlacement( 0, myPMTHeadVacPos,  myPMTTopVacName[i],  fLogicPMTTopVac[i],      fPhysicPMTTop[i],         false, 0, myCheckOverlap );
    fPhysicPMTLArDisk[i]     = new G4PVPlacement( 0, myPMTLArDiskPos,  myPMTLArDiskName[i], fLogicPMTLArDisk[i],     fPhysicPMTSurrounding[i], false, 0, myCheckOverlap );
    fPhysicPMTStem[i]        = new G4PVPlacement( 0, myPMTStemPos,     myPMTStemName[i],    fLogicPMTStem[i],        fPhysicPMTTop[i],         false, 0, myCheckOverlap );
    
  }


  DefineSurfaces();

}


DSDetectorPMTDS50::~DSDetectorPMTDS50(){
  ;
}


void  DSDetectorPMTDS50::DefineSurfaces(){

  // Photocathode - Vacuum
  fOpPMTVacuumSurface = new G4OpticalSurface("OpPMTLArSurface");
  for(int i = 0; i < 19; i++)  fPMTVacuumSurface[i] = new G4LogicalBorderSurface("PMTVacuumSurface", fPhysicPMTHeadVac[i], fPhysicPMTWindow[i], fOpPMTVacuumSurface ); 
  fOpPMTVacuumSurface->SetType( dielectric_metal );
  fOpPMTVacuumSurface->SetModel( glisur );
  fOpPMTVacuumSurface->SetFinish( polished );  
  fPMTVacuumSurfProp = new G4MaterialPropertiesTable();
  fPMTVacuumSurfProp->AddConstProperty("REFLECTIVITY", 1.0);
  fPMTVacuumSurfProp->AddConstProperty("EFFICIENCY",   0.0);
  fOpPMTVacuumSurface->SetMaterialPropertiesTable( fPMTVacuumSurfProp );
  

  // Photocathode - LAr
  fOpPMTLArSurface = new G4OpticalSurface("OpPMTLArSurface");
  for(int i = 0; i < 19; i++)  fPMTLArSurface[i] = new G4LogicalBorderSurface("PMTLArSurface", fPhysicPMTLArDisk[i], fPhysicPMTWindow[i], fOpPMTLArSurface );  
  //for(int i = 0; i < 19; i++)  fPMTLArSurface[i] = new G4LogicalBorderSurface("PMTLArSurface", fPhysicPMTWindow[i], fPhysicPMTLArDisk[i],  fOpPMTLArSurface );  
  fOpPMTLArSurface->SetType( dielectric_dielectric );
  fOpPMTLArSurface->SetModel( unified );
  fOpPMTLArSurface->SetFinish( polished );  
  //fOpPMTLArSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetPhotocathodeMPT());  
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
  for( int i = 0; i < 19; i++ )  fTeflonLArSurface[i] = new G4LogicalBorderSurface("TeflonLArSurface", fTeflonCapDiskVolume, fPhysicPMTLArDisk[i], fOpTeflonLArSurface );
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
 * $Log: DSDetectorPMTDS50.cc,v $
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
 * Fixed a bug in the photoelectron absorption in DS50. Back to the previous QE method
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
