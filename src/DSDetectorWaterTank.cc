#include "DSDetectorWaterTank.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSMaterial.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include <iostream>
#include <fstream>
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4UIcommand.hh"

using namespace std;

DSDetectorWaterTank::DSDetectorWaterTank(G4VPhysicalVolume *myMotherVolume) {
  
  fMotherVolume = myMotherVolume;  
  
  DSLog(routine) << " Constructing WaterTank Geometry" << endlog ;
  
  //Stainless Steel Tank
  
  fSolidSteelTank  = new G4Tubs("SteelTank_Solid", 0,5000*mm,4613*mm,0,twopi*rad);
  fLogicSteelTank  = new G4LogicalVolume(fSolidSteelTank, DSMaterial::Get()->GetStainlessSteel(), "SteelTank_Logic");
  fPhysicSteelTank = new G4PVPlacement(0,
				 G4ThreeVector(0,0,0),
				 "SteelTank",
				 fLogicSteelTank,
				 fMotherVolume,
				 false,
				 0,
				 DSStorage::Get()->GetCheckOverlap());

  //Water volume

  fSolidWaterVolume  = new G4Tubs("WaterVolume_Solid", 0,4995*mm,4608*mm,0,twopi*rad);
  fLogicWaterVolume  = new G4LogicalVolume(fSolidWaterVolume, DSMaterial::Get()->GetWater(), "WaterVolume_Logic");
  fPhysicWaterVolume = new G4PVPlacement(0,
				 G4ThreeVector(0,0,0),
				 "WaterVolume",
				 fLogicWaterVolume,
				 fPhysicSteelTank,
				 false,
				 811,
				 DSStorage::Get()->GetCheckOverlap());

  DSLog(routine) << fPhysicWaterVolume->GetName() << " = " << fPhysicWaterVolume->GetCopyNo() << endlog;

  //Neutron sphere legs and structure

  G4double cube_h=400*mm;
  G4double cube_x=500*mm;
  G4double cube_y=450*mm;
  G4double leg_h=4583*mm; 
  G4double leg_r=340/2.*mm;
  
  G4double Bbar_l=2829*mm;
  G4double Bbar_h=60*mm;
  G4double Bbar_th=30*mm;
  G4double Bbar_z=-4094*mm;

  G4double Hbar_l=2829*mm;
  G4double Hbar_h=120*mm;
  G4double Hbar_th=60*mm;
  G4double Hbar_z=-2360*mm;

  G4double Dbar_l=60*mm;
  //G4double Dbar_h=1558.2*mm;
  G4double Dbar_h=1557.2*mm;
  G4double Dbar_th=30*mm;
  G4double DbarU_z=-1520.9*mm;
  G4double DbarB_z=-3199.1*mm;
  G4double Dbar_theta=41*twopi*rad/360.;

  G4double disc_th=25*mm;
  G4double disc_r=650*mm;

  G4double wing_l=460*mm;
  G4double wing_th=12*mm;


  fSolidNVlegUp    = new G4Box("NVlegUp_solid",cube_x/2.,cube_y/2.,cube_h/2.);
  fLogicNVlegUp    = new G4LogicalVolume(fSolidNVlegUp,DSMaterial::Get()->GetStainlessSteel(),"NVlegUp_Logic");
    
  fSolidNVleg      = new G4Tubs("NVleg_Solid",0,leg_r,(leg_h-cube_h/2.)/2.,0.,twopi*rad);
  fLogicNVleg      = new G4LogicalVolume(fSolidNVleg,DSMaterial::Get()->GetStainlessSteel(),"NVleg_Logic");
  
  //bottom horizontal bar
  fSolidNVbarB    = new G4Box("NVbarB_solid",Bbar_l/2.,Bbar_th/2.,Bbar_h/2.);
  fLogicNVbarB    = new G4LogicalVolume(fSolidNVbarB,DSMaterial::Get()->GetStainlessSteel(),"NVbarB_Logic");
  //upper horizontal bar
  fSolidNVbarH    = new G4Box("NVbarH_solid",Hbar_l/2.,Hbar_th/2.,Hbar_h/2.);
  fLogicNVbarH    = new G4LogicalVolume(fSolidNVbarH,DSMaterial::Get()->GetStainlessSteel(),"NVbarH_Logic");
  //diagonal right bar
  fSolidNVbarDR    = new G4Para("NVbarDR_solid",Dbar_l/2.,Dbar_th/2.,Dbar_h/2.,0,Dbar_theta,0);
  fLogicNVbarDR    = new G4LogicalVolume(fSolidNVbarDR,DSMaterial::Get()->GetStainlessSteel(),"NVbarDR_Logic");
 //diagonal left bar
  fSolidNVbarDL    = new G4Para("NVbarDL_solid",Dbar_l/2.,Dbar_th/2.,Dbar_h/2.,0,Dbar_theta,0);
  fLogicNVbarDL    = new G4LogicalVolume(fSolidNVbarDL,DSMaterial::Get()->GetStainlessSteel(),"NVbarDL_Logic");
  //disc base
  fSolidNVdisc      = new G4Tubs("NVdisc_Solid",0,disc_r,disc_th/2.,0.,twopi*rad);
  fLogicNVdisc      = new G4LogicalVolume(fSolidNVdisc,DSMaterial::Get()->GetStainlessSteel(),"NVdisc_Logic");
  //wing base
  fSolidNVwing      = new G4Trap("NVwing_Solid",wing_l,wing_l,wing_th,wing_th);
  fLogicNVwing      = new G4LogicalVolume(fSolidNVwing,DSMaterial::Get()->GetStainlessSteel(),"NVwing_Logic");


  G4double phi[5];
  phi[0]=67.5*twopi*rad/360.;
  phi[1]=157.5*twopi*rad/360.;
  phi[2]=247.5*twopi*rad/360.;
  phi[3]=337.5*twopi*rad/360.;
  phi[4]=(360+67.5)*twopi*rad/360.;

  G4double rleg=2241*mm;
  
  for(int n=0;n<4;n++)
    {
      G4String number = G4UIcommand::ConvertToString(n+1);
      G4String nameleg="VNleg_"+number;
      G4String namelegUp="VNlegUp_"+number;
      G4String namebarB="VNbarB_"+number;
      G4String namebarH="VNbarH_"+number;
      G4String namebarDRU="VNbarDRU_"+number;
      G4String namebarDRB="VNbarDRB_"+number;
      G4String namebarDLU="VNbarDLU_"+number;
      G4String namebarDLB="VNbarDLB_"+number;
      G4String namedisc="VNdisc_"+number;
      G4String namewing1="VNwing1_"+number;
      G4String namewing2="VNwing2_"+number;
      G4String namewing3="VNwing3_"+number;
      G4String namewing4="VNwing4_"+number;

      G4ThreeVector RelocationLeg(rleg*sin(phi[n]),rleg*cos(phi[n]),-(leg_h+cube_h/2.)/2.);
      G4ThreeVector RelocationLegUp(rleg*sin(phi[n]),rleg*cos(phi[n]),0);
      G4double rbar=rleg*cos((phi[n+1]-phi[n])/2.);
      G4ThreeVector RelocationBarB(rbar*sin((phi[n]+phi[n+1])/2.),rbar*cos((phi[n]+phi[n+1])/2.),Bbar_z);
      G4ThreeVector RelocationBarH(rbar*sin((phi[n]+phi[n+1])/2.),rbar*cos((phi[n]+phi[n+1])/2.),Hbar_z);
      G4double rbarD=sqrt(rbar*rbar + (707.25*mm)*(707.25*mm));
      G4double angleD=acos(rbar/rbarD);      
      G4ThreeVector RelocationBarDRU(rbarD*sin((phi[n]+phi[n+1])/2.+angleD ),rbarD*cos((phi[n]+phi[n+1])/2.+angleD),DbarU_z);
      G4ThreeVector RelocationBarDRB(rbarD*sin((phi[n]+phi[n+1])/2.-angleD ),rbarD*cos((phi[n]+phi[n+1])/2.-angleD),DbarB_z);
      G4ThreeVector RelocationBarDLU(rbarD*sin((phi[n]+phi[n+1])/2.-angleD ),rbarD*cos((phi[n]+phi[n+1])/2.-angleD),DbarU_z);
      G4ThreeVector RelocationBarDLB(rbarD*sin((phi[n]+phi[n+1])/2.+angleD ),rbarD*cos((phi[n]+phi[n+1])/2.+angleD),DbarB_z);
      G4ThreeVector RelocationDisc(rleg*sin(phi[n]),rleg*cos(phi[n]),-leg_h-disc_th/2.);
      G4double rwing=wing_l/2.+leg_r;
      G4ThreeVector RelocationWing1(rleg*sin(phi[n])+rwing*sin(phi[n]),rleg*cos(phi[n])+rwing*cos(phi[n]),(-leg_h+wing_l/2.));
      G4ThreeVector RelocationWing2(rleg*sin(phi[n])+rwing*sin(phi[n]+90*twopi*rad/360.),rleg*cos(phi[n])+rwing*cos(phi[n]+90*twopi*rad/360.),(-leg_h+wing_l/2.));
      G4ThreeVector RelocationWing3(rleg*sin(phi[n])+rwing*sin(phi[n]+180*twopi*rad/360.),rleg*cos(phi[n])+rwing*cos(phi[n]+180*twopi*rad/360.),(-leg_h+wing_l/2.));
      G4ThreeVector RelocationWing4(rleg*sin(phi[n])+rwing*sin(phi[n]+270*twopi*rad/360.),rleg*cos(phi[n])+rwing*cos(phi[n]+270*twopi*rad/360.),(-leg_h+wing_l/2.));


      G4RotationMatrix* RotationLegUp = new G4RotationMatrix;
      RotationLegUp->rotateZ((phi[n]));

      G4RotationMatrix* RotationBar = new G4RotationMatrix;
      RotationBar->rotateZ((phi[n]+phi[n+1])/2.);

      G4RotationMatrix* RotationBarL = new G4RotationMatrix;
      RotationBarL->rotateZ((phi[n]+phi[n+1])/2. + pi*rad);

      G4RotationMatrix* RotationWing1 = new G4RotationMatrix;
      RotationWing1->rotateZ((phi[n]));

      G4RotationMatrix* RotationWing2 = new G4RotationMatrix;
      RotationWing2->rotateZ((phi[n])+90*twopi*rad/360.);

      G4RotationMatrix* RotationWing3 = new G4RotationMatrix;
      RotationWing3->rotateZ((phi[n])+180*twopi*rad/360.);
      
      G4RotationMatrix* RotationWing4 = new G4RotationMatrix;
      RotationWing4->rotateZ((phi[n])+270*twopi*rad/360.);


      fPhysicNVleg   = new G4PVPlacement(0,
					 RelocationLeg,
					 nameleg,
					 fLogicNVleg,
					 fPhysicWaterVolume,
					 true,
					 n,
					 DSStorage::Get()->GetCheckOverlap());
      
      fPhysicNVlegUp  = new G4PVPlacement(RotationLegUp,
					  RelocationLegUp,
					  namelegUp,
					  fLogicNVlegUp,
					  fPhysicWaterVolume,
					  true,
					  n,
					  DSStorage::Get()->GetCheckOverlap());

      fPhysicNVbarB   = new G4PVPlacement(RotationBar,
					  RelocationBarB,
					  namebarB,
					  fLogicNVbarB,
					  fPhysicWaterVolume,
					  true,
					  n,
					  DSStorage::Get()->GetCheckOverlap());

      fPhysicNVbarH   = new G4PVPlacement(RotationBar,
					  RelocationBarH,
					  namebarH,
					  fLogicNVbarH,
					  fPhysicWaterVolume,
					  true,
					  n,
					  DSStorage::Get()->GetCheckOverlap());

      fPhysicNVbarDR   = new G4PVPlacement(RotationBar,
					  RelocationBarDRU,
					  namebarDRU,
					  fLogicNVbarDR,
					  fPhysicWaterVolume,
					  true,
					  n,
					  DSStorage::Get()->GetCheckOverlap());

      fPhysicNVbarDR   = new G4PVPlacement(RotationBar,
					   RelocationBarDRB,
					   namebarDRB,
					   fLogicNVbarDR,
					   fPhysicWaterVolume,
					   true,
					   n,
					   DSStorage::Get()->GetCheckOverlap());

      fPhysicNVbarDL   = new G4PVPlacement(RotationBarL,
					   RelocationBarDLU,
					   namebarDLU,
					   fLogicNVbarDL,
					   fPhysicWaterVolume,
					   true,
					   n,
					   DSStorage::Get()->GetCheckOverlap());

      fPhysicNVbarDL   = new G4PVPlacement(RotationBarL,
					   RelocationBarDLB,
					   namebarDLB,
					   fLogicNVbarDL,
					   fPhysicWaterVolume,
					   true,
					   n,
					   DSStorage::Get()->GetCheckOverlap());

      fPhysicNVdisc  = new G4PVPlacement(0,
					 RelocationDisc,
					 namedisc,
					 fLogicNVdisc,
					 fPhysicWaterVolume,
					 true,
					 n,
					 DSStorage::Get()->GetCheckOverlap());

      fPhysicNVwing  = new G4PVPlacement(RotationWing1,
					 RelocationWing1,
					 namewing1,
					 fLogicNVwing,
					 fPhysicWaterVolume,
					 true,
					 n,
					 DSStorage::Get()->GetCheckOverlap());
      
      fPhysicNVwing  = new G4PVPlacement(RotationWing2,
					 RelocationWing2,
					 namewing2,
					 fLogicNVwing,
					 fPhysicWaterVolume,
					 true,
					 n,
					 DSStorage::Get()->GetCheckOverlap()); 

      fPhysicNVwing  = new G4PVPlacement(RotationWing3,
					 RelocationWing3,
					 namewing3,
					 fLogicNVwing,
					 fPhysicWaterVolume,
					 true,
					 n,
					 DSStorage::Get()->GetCheckOverlap()); 

      fPhysicNVwing  = new G4PVPlacement(RotationWing4,
					 RelocationWing4,
					 namewing4,
					 fLogicNVwing,
					 fPhysicWaterVolume,
					 true,
					 n,
					 DSStorage::Get()->GetCheckOverlap());    

}
  
}
DSDetectorWaterTank::~DSDetectorWaterTank(){
  ; //delete fMessenger;
}


/*
 * $Log: DSDetectorWaterTank.cc,v $
 * Revision 1.1  2014/05/07 12:21:02  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.7  2013/03/26 13:57:38  meregaglia
 * fixed small overlap between diagonal bars of the NV support
 *
 * Revision 1.6  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
