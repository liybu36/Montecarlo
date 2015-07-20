#include <iostream>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4PVPlacement.hh"
#include "G4UIcommand.hh"

#include "DSDetectorPMTDSG2.hh"
#include "DSStorage.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"



using namespace std;


DSDetectorPMTDSG2::DSDetectorPMTDSG2( G4VPhysicalVolume* myMotherVolume, G4VPhysicalVolume* myTeflonCapDiskVolume ){

  fMotherVolume        = myMotherVolume;
  fTeflonCapDiskVolume = myTeflonCapDiskVolume;


  const double myTwoPi    = 2*M_PI*rad;
  G4bool   myCheckOverlap = DSStorage::Get()->GetCheckOverlap();



  // Variables to distinguish bewteen top and bottom
  G4int  myOffset = NPMTS;
  if( fMotherVolume->GetRotation() == 0 ) {
    myOffset = 0;
    DSLog(routine) << " Constructing Top TPC PMTs Geometry  " << endlog;
  }
  else DSLog(routine) << " Constructing Bottom TPC PMTs Geometry  " << endlog;

  // Teflon Cap
  G4double myPMTAssembly_h = 130.*mm;//DSParameters::Get()->GetPMTAssemblyHeight()/mm;
  G4double myTeflonCap_h = DSParameters::Get()->GetTeflonCapHeight()/mm;
  G4double myPMTBody_d = DSParameters::Get()->GetPMTBodyDiameter()/mm;
  G4double myPMTBody_h = DSParameters::Get()->GetPMTBodyHeight()/mm;
  G4double myPMTHead_d = DSParameters::Get()->GetPMTHeadDiameter()/mm;
  G4double myPMTHead_h = DSParameters::Get()->GetPMTHeadHeight()/mm;
  G4double myPMTWall   = DSParameters::Get()->GetPMTWallThickness()/mm; 
  G4double myPMTWindowThickness = DSParameters::Get()->GetPMTWindowThickness()/mm;
  G4double myPMTOffset = DSParameters::Get()->GetPMTOffset()/mm;

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
  double mySurrounding_h = myTeflonCap_h/2.;
  double mySurrounding_r = myPMTHead_d/2.;
  fSolidPMTSurrounding = new G4Tubs( "PMTSurrounding_Solid", 
                                      0,  
                                      mySurrounding_r,	     
                                      mySurrounding_h,      				    
				      0, 
				      myTwoPi );
  // Part of the PMT inside the surrounding (Head, vacuum, window and body) 
  fSolidPMTHead        = new G4Tubs( "PMTHead_Solid",        
                                      0,  
				      mySurrounding_r,              
				      myPMTHead_h/2.,
				      0, 
				      myTwoPi );
				      
  fSolidPMTHeadVac     = new G4Tubs( "PMTHeadVac_Solid",     
                                      0,  
				      mySurrounding_r - myPMTWall,  
				      myPMTHead_h/2. - myPMTWall/2.,         
				      0, 
				      myTwoPi );
				      
  fSolidPMTWindow      = new G4Tubs( "PMTWindow_Solid",      
                                      0,  
				      mySurrounding_r - myPMTWall,  
				      myPMTWindowThickness/2., 
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
				      mySurrounding_h - myPMTHead_h/2. - myPMTOffset/2.,   
				      0, 
				      myTwoPi ); 
				      
  fSolidPMTBodyVac     = new G4Tubs( "PMTBodyVac_Solid",     
                                      0,  
				      myPMTBody_d/2. - myPMTWall,  
				      mySurrounding_h - myPMTHead_h/2. - myPMTOffset/2.,                
				      0, 
				      myTwoPi );
				      
  // Part of the PMT body outside the surrounding
  fSolidPMTTop         = new G4Tubs( "PMTTop_Solid",         
                                      0,  
				      myPMTBody_d/2.,              
				      myPMTBody_h/2. - mySurrounding_h + myPMTOffset/2.,                
				      0, 
				      myTwoPi );
				      
  fSolidPMTTopVac      = new G4Tubs( "PMTTopVac_Solid",      
                                      0,  
				      myPMTBody_d/2. - myPMTWall,  
				      myPMTBody_h/2. - mySurrounding_h + myPMTOffset/2. - myPMTWall/2., 
				      0, 
				      myTwoPi );
				      
  fSolidPMTLArDisk     = new G4Tubs( "PMTLArDisk_Solid",     
                                      0,  
				      mySurrounding_r,              
				      myPMTOffset/2.,                                                    
				      0, 
				      myTwoPi );



  // Positioning
  G4double myPMTRGap = 0*6.19*mm;
  G4double myPMTPhiGap = 0*6.5*mm;
  G4double myPMTSpacing = 76.2*mm;// DSParameters::Get()->GetPMTSpacing()/mm;
  myPMTSpacing += myPMTRGap;
  G4ThreeVector myPMTTopZ( 0, 0, -myPMTAssembly_h/2. + myPMTBody_h/2. + myPMTOffset/2. + mySurrounding_h );
  G4ThreeVector myPMTPos[NPMTS];
  G4ThreeVector myPMTRad1( myPMTSpacing, 0, 0 );
  G4ThreeVector myPMTRad2 = myPMTRad1 + G4ThreeVector( myPMTSpacing, 0, 0).rotateZ( myTwoPi/6. );
  G4ThreeVector myPMTRad3 = 2*myPMTRad1;
  G4ThreeVector myPMTRad4 = myPMTRad2 + G4ThreeVector( myPMTSpacing, 0, 0).rotateZ( myTwoPi/6. );
  G4ThreeVector myPMTRad5 = myPMTRad3 + G4ThreeVector( myPMTSpacing, 0, 0).rotateZ( myTwoPi/6. );
  G4ThreeVector myPMTRad6 = myPMTRad4 + G4ThreeVector( myPMTSpacing, 0, 0).rotateZ( myTwoPi/6. );
  G4ThreeVector myPMTRad7 = myPMTRad5 + G4ThreeVector( myPMTSpacing, 0, 0).rotateZ( myTwoPi/6. );
  G4ThreeVector myPMTRad8 = myPMTRad6 + G4ThreeVector( myPMTSpacing, 0, 0).rotateZ( myTwoPi/6. );
  G4ThreeVector myPMTRad( 0, 0, 0);
  // 3 circles of PMTs at radii myPMTRad1, myPMTRad2, myPMTRad3
  myPMTPos[0].set( 0, 0, 0);
  G4int pmtnum = 1;
  G4double radius = 0;
  G4double circumference = 0;
  G4int npmtsInCirc = 0;
  G4double pmtDiam = 76.2*mm;
  G4double teflonIR = 1483.*mm/2.;
  while(pmtnum < NPMTS)
    {
      radius += myPMTSpacing;
      if(radius + mySurrounding_r > teflonIR)
	{
	  DSLog(routine) << "radius too large : " << radius+mySurrounding_r << ", where cap rad = " << teflonIR << endlog; 
	}
      circumference = myTwoPi*radius;
      npmtsInCirc = (int)(circumference/(pmtDiam+myPMTPhiGap));
      if(npmtsInCirc > (NPMTS-pmtnum))
	npmtsInCirc = NPMTS-pmtnum;
      for(int i = 0; i < npmtsInCirc; i++)
	{
	  myPMTPos[pmtnum++] = G4ThreeVector(radius, 0, 0).rotateZ(myTwoPi*i/npmtsInCirc); 
	}
    }

    /*
  while(pmtnum < NPMTS)
    {
      myPMTRad += myPMTRad1;
      std::cout << "PMT : " << pmtnum << "\t(" << myPMTRad.getX() << ", " << myPMTRad.getY() << ", " << myPMTRad.getZ() << ")\t" << myPMTRad.mag() << std::endl;
      myPMTPos[pmtnum++] = myPMTRad;
      for(int i = 0; i < 5; i++)
	{
	  myPMTPos[pmtnum++] = myPMTRad.rotateZ( myTwoPi/6. );
	  std::cout << "PMT : " << pmtnum << "\t(" << myPMTRad.getX() << ", " << myPMTRad.getY() << ", " << myPMTRad.getZ() << ")\t" << myPMTRad.mag() << std::endl;
	}
      //      myPMTRad1 += G4ThreeVector( myPMTSpacing, 0, 0);
      if(pmtnum % == 0)
      	myPMTRad1.rotateZ( myTwoPi/6. );
    }
    */
  /*
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
  myPMTPos[19] = myPMTRad4;
  myPMTPos[20] = myPMTRad4.rotateZ( myTwoPi/6. );
  myPMTPos[21] = myPMTRad4.rotateZ( myTwoPi/6. );
  myPMTPos[22] = myPMTRad4.rotateZ( myTwoPi/6. );
  myPMTPos[23] = myPMTRad4.rotateZ( myTwoPi/6. );
  myPMTPos[24] = myPMTRad4.rotateZ( myTwoPi/6. );
  myPMTPos[25] = myPMTRad5;
  myPMTPos[26] = myPMTRad5.rotateZ( myTwoPi/6. );
  myPMTPos[27] = myPMTRad5.rotateZ( myTwoPi/6. );
  myPMTPos[28] = myPMTRad5.rotateZ( myTwoPi/6. );
  myPMTPos[29] = myPMTRad5.rotateZ( myTwoPi/6. );
  myPMTPos[30] = myPMTRad5.rotateZ( myTwoPi/6. );
  myPMTPos[31] = myPMTRad6;
  myPMTPos[32] = myPMTRad6.rotateZ( myTwoPi/6. );
  myPMTPos[33] = myPMTRad6.rotateZ( myTwoPi/6. );
  myPMTPos[34] = myPMTRad6.rotateZ( myTwoPi/6. );
  myPMTPos[35] = myPMTRad6.rotateZ( myTwoPi/6. );
  myPMTPos[36] = myPMTRad6.rotateZ( myTwoPi/6. );
  myPMTPos[37] = myPMTRad7;
  myPMTPos[38] = myPMTRad7.rotateZ( myTwoPi/6. );
  myPMTPos[39] = myPMTRad7.rotateZ( myTwoPi/6. );
  myPMTPos[40] = myPMTRad7.rotateZ( myTwoPi/6. );
  myPMTPos[41] = myPMTRad7.rotateZ( myTwoPi/6. );
  myPMTPos[42] = myPMTRad7.rotateZ( myTwoPi/6. );
  myPMTPos[43] = myPMTRad8;
  myPMTPos[44] = myPMTRad8.rotateZ( myTwoPi/6. );
  myPMTPos[45] = myPMTRad8.rotateZ( myTwoPi/6. );
  myPMTPos[46] = myPMTRad8.rotateZ( myTwoPi/6. );
  myPMTPos[47] = myPMTRad8.rotateZ( myTwoPi/6. );
  myPMTPos[48] = myPMTRad8.rotateZ( myTwoPi/6. );
  */
  // shifts for assemblying the pmt
  G4ThreeVector myZeros        ( 0, 0, 0 ); 
  G4ThreeVector myPMTHeadPos   ( 0, 0, -mySurrounding_h + myPMTHead_h/2. + myPMTOffset );
  G4ThreeVector myPMTHeadVacPos( 0, 0, -myPMTWall/2. );
  G4ThreeVector myPMTWindowPos ( 0, 0, -myPMTHead_h/2. + myPMTWall/2. + myPMTWindowThickness/2. );
  G4ThreeVector myPMTJoinVacPos( 0, 0,  myPMTHead_h/2. - myPMTWall/2. );
  G4ThreeVector myPMTBodyPos   ( 0, 0,  myPMTHead_h/2. + myPMTOffset/2. );
  G4ThreeVector myPMTTopPos[NPMTS];
  G4ThreeVector myPMTLArDiskPos( 0, 0, -mySurrounding_h + myPMTOffset/2. );
  
  /*cout << myPMTHeadPos <<  " "
  << myPMTHeadVacPos <<  " "
  << myPMTWindowPos <<  " "
  << myPMTJoinVacPos <<  " "
  << myPMTBodyPos <<  " "
  << myPMTLArDiskPos <<  " " << endl ;
  */
   
  G4String myPMTSurLogicName[NPMTS];
  G4String myPMTHeadLogicName[NPMTS];
  G4String myPMTHeadVacLogicName[NPMTS];
  G4String myPMTWindowLogicName[NPMTS];
  G4String myPMTJoinVacLogicName[NPMTS];
  G4String myPMTBodyLogicName[NPMTS];
  G4String myPMTBodyVacLogicName[NPMTS];
  G4String myPMTTopLogicName[NPMTS];
  G4String myPMTTopVacLogicName[NPMTS];
  G4String myPMTLArDiskLogicName[NPMTS];

  G4String myPMTSurName[NPMTS];
  G4String myPMTHeadName[NPMTS];
  G4String myPMTHeadVacName[NPMTS];
  G4String myPMTWindowName[NPMTS];
  G4String myPMTJoinVacName[NPMTS];
  G4String myPMTBodyName[NPMTS];
  G4String myPMTBodyVacName[NPMTS];
  G4String myPMTTopName[NPMTS];
  G4String myPMTTopVacName[NPMTS];
  G4String myPMTLArDiskName[NPMTS];

  for( int i = 0; i < NPMTS; i++) {

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

    myPMTTopPos[i] = myPMTPos[i] + myPMTTopZ;


    fLogicPMTSurrounding[i] = new G4LogicalVolume( fSolidPMTSurrounding, DSMaterial::Get()->GetTeflon(),      myPMTSurLogicName[i] );
    fLogicPMTHead[i]        = new G4LogicalVolume( fSolidPMTHead,        DSMaterial::Get()->GetKovar(),       myPMTHeadLogicName[i] ); 
    fLogicPMTHeadVac[i]     = new G4LogicalVolume( fSolidPMTHeadVac,     DSMaterial::Get()->GetVacuum(),      myPMTHeadVacLogicName[i] );
    fLogicPMTWindow[i]      = new G4LogicalVolume( fSolidPMTWindow,      DSMaterial::Get()->GetBialkali(),    myPMTWindowLogicName[i] );
    fLogicPMTJoinVac[i]     = new G4LogicalVolume( fSolidPMTJoinVac,     DSMaterial::Get()->GetVacuum(),      myPMTJoinVacLogicName[i] );
    fLogicPMTBody[i]        = new G4LogicalVolume( fSolidPMTBody,        DSMaterial::Get()->GetKovar(),       myPMTBodyLogicName[i] );
    fLogicPMTBodyVac[i]     = new G4LogicalVolume( fSolidPMTBodyVac,     DSMaterial::Get()->GetVacuum(),      myPMTBodyVacLogicName[i] );
    fLogicPMTTop[i]         = new G4LogicalVolume( fSolidPMTTop,         DSMaterial::Get()->GetKovar(),       myPMTTopLogicName[i] );
    fLogicPMTTop[i]->SetVisAttributes(new G4VisAttributes(G4Colour(1,0,0)));
    fLogicPMTTopVac[i]      = new G4LogicalVolume( fSolidPMTTopVac,      DSMaterial::Get()->GetVacuum(),      myPMTTopVacLogicName[i] ); 
    fLogicPMTLArDisk[i]     = new G4LogicalVolume( fSolidPMTLArDisk,     DSMaterial::Get()->GetLiquidArgon(), myPMTLArDiskLogicName[i] );  

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

  }


  DefineSurfaces();

}


DSDetectorPMTDSG2::~DSDetectorPMTDSG2(){
  ;
}


void  DSDetectorPMTDSG2::DefineSurfaces(){

  // Photocathode - Vacuum
  fOpPMTVacuumSurface = new G4OpticalSurface("OpPMTLArSurface");
  for(int i = 0; i < NPMTS; i++)  fPMTVacuumSurface[i] = new G4LogicalBorderSurface("PMTVacuumSurface", fPhysicPMTHeadVac[i], fPhysicPMTWindow[i], fOpPMTVacuumSurface ); 
  fOpPMTVacuumSurface->SetType( dielectric_metal );
  fOpPMTVacuumSurface->SetModel( glisur );
  fOpPMTVacuumSurface->SetFinish( polished );  
  fPMTVacuumSurfProp = new G4MaterialPropertiesTable();
  fPMTVacuumSurfProp->AddConstProperty("REFLECTIVITY", 1.0);
  fPMTVacuumSurfProp->AddConstProperty("EFFICIENCY",   0.0);
  fOpPMTVacuumSurface->SetMaterialPropertiesTable( fPMTVacuumSurfProp );
  

  // Photocathode - LAr
  fOpPMTLArSurface = new G4OpticalSurface("OpPMTLArSurface");
  for(int i = 0; i < NPMTS; i++)  fPMTLArSurface[i] = new G4LogicalBorderSurface("PMTLArSurface", fPhysicPMTLArDisk[i], fPhysicPMTWindow[i], fOpPMTLArSurface );  
  //for(int i = 0; i < 19; i++)  fPMTLArSurface[i] = new G4LogicalBorderSurface("PMTLArSurface", fPhysicPMTWindow[i], fPhysicPMTLArDisk[i],  fOpPMTLArSurface );  
  fOpPMTLArSurface->SetType( dielectric_dielectric );
  fOpPMTLArSurface->SetModel( unified );
  fOpPMTLArSurface->SetFinish( polished );  
  //fOpPMTLArSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetPhotocathodeMPT());  
  fPMTLArSurfProp = new G4MaterialPropertiesTable();
  fPMTLArSurfProp->AddConstProperty("REFLECTIVITY", 0.2);
  fPMTLArSurfProp->AddConstProperty("EFFICIENCY",   0.0);
  fOpPMTLArSurface->SetMaterialPropertiesTable( fPMTLArSurfProp );
  


  // Teflon - LAr
  fOpTeflonLArSurface = new G4OpticalSurface("OpTeflonLArSurface");
  for( int i = 0; i < NPMTS; i++ )  fTeflonLArSurface[i] = new G4LogicalBorderSurface("TeflonLArSurface", fTeflonCapDiskVolume, fPhysicPMTLArDisk[i], fOpTeflonLArSurface );
  fOpTeflonLArSurface->SetType(dielectric_dielectric);
  fOpTeflonLArSurface->SetModel(unified);
  fOpTeflonLArSurface->SetSigmaAlpha(0.1);
  //where sigma_alpha is in [rad]
  fOpTeflonLArSurface->SetFinish(polished);
  fTeflonLArSurfProp = new G4MaterialPropertiesTable();
  fTeflonLArSurfProp->AddConstProperty("REFLECTIVITY", 0.2);
  fTeflonLArSurfProp->AddConstProperty("EFFICIENCY",   0.0);
  fOpTeflonLArSurface->SetMaterialPropertiesTable( fTeflonLArSurfProp );

}
