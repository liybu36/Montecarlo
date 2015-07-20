#include "DSMaterial.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "DSParameters.hh"
#include "DSEventHandler.hh"
#include <fstream>

using namespace std;

DSMaterial* DSMaterial::me =0;

DSMaterial::DSMaterial() {
  DefineMaterials();
  DefineBoundaries();
  DefineProperties();
  SetMaterialIndexes();
}

DSMaterial::~DSMaterial(){}

DSMaterial* DSMaterial::Get() {
  if (!me) 
    me = new DSMaterial();
       
  return me;
}

void DSMaterial::DefineMaterials() {

  G4String  name,
            symbol;
  G4int     z,   
            ncomponents,     
            natoms; 
  G4double  a,
            density,
            fractionmass,
            temperature,
            pressure;



  G4Element* H  = new G4Element(name="Hydrogen",  symbol="H" , z=1,  a=1.00794*g/mole );
  G4Element* C  = new G4Element(name="Carbon",    symbol="C",  z=6,  a=12.0107*g/mole ); 
  //G4Element* B =  new G4Element(name="Boron",     symbol="B",  z=5,  a=10.811*g/mole  );  
  G4Element* N  = new G4Element(name="Nitrogen",  symbol="N",  z=7,  a=14.00067*g/mole);
  G4Element* O  = new G4Element(name="Oxygen",    symbol="O",  z=8,  a=15.9994*g/mole ); 
  G4Element* Fl = new G4Element(name="Fluorine",  symbol="F",  z=9,  a=18.9984*g/mole );
  G4Element* Na = new G4Element(name="Sodium",    symbol="Na", z=11, a=22.9898*g/mole );
  G4Element* Mg = new G4Element(name="Magnesium", symbol="Mg", z=12, a=24.305*g/mole  );
  G4Element* Al = new G4Element(name="Aluminum",  symbol="Al", z=13, a=26.9815*g/mole );
  G4Element* Si = new G4Element(name="Silicium",  symbol="Si", z=14, a=28.0855*g/mole ); 
  G4Element* S  = new G4Element(name="Sulfur",    symbol="S",  z=16, a=32.065*g/mole  );
  G4Element* Ar = new G4Element(name="Argon",     symbol="Ar", z=18, a=39.948*g/mole  );  
  G4Element* K  = new G4Element(name="Potassium", symbol="K",  z=19, a=39.0989*g/mole );
  G4Element* Ca = new G4Element(name="Calcium",   symbol="Ca", z=20, a=40.078*g/mole  );
  G4Element* Ti = new G4Element(name="Titanium",  symbol="Ti", z=22, a=47.867*g/mole  );
  G4Element* Cr = new G4Element(name="Chromium",  symbol="Cr", z=24, a=51.9961*g/mole );  
  G4Element* Mn = new G4Element(name="Manganese", symbol="Mn", z=25, a=54.93805*g/mole);  
  G4Element* Fe = new G4Element(name="Iron",      symbol="Fe", z=26, a=55.845*g/mole  );
  G4Element* Co = new G4Element(name="Cobalt",    symbol="Co", z=27, a=58.9331*g/mole );
  G4Element* Ni = new G4Element(name="Nickel",    symbol="Ni", z=28, a=58.70*g/mole   );
  G4Element* Cu = new G4Element(name="Copper",    symbol="Cu", z=29, a=63.546*g/mole  );
  G4Element* P  = new G4Element(name="Phosphorus",symbol="P",  z=30, a=30.97398*g/mole);
  G4Element* Rb = new G4Element(name="Rubidium",  symbol="Rb", z=37, a=85.4678*g/mole );
  G4Element* In = new G4Element(name="Indium",    symbol="In", z=49, a=114.818*g/mole );
  G4Element* Sn = new G4Element(name="Selenium",  symbol="Sn", z=50, a=118.710*g/mole ); // Tin????
  G4Element* Sb = new G4Element(name="Antimony",  symbol="Sb", z=51, a=121.760*g/mole );
  G4Element* Pb = new G4Element(name="Lead",      symbol="Pb", z=82, a=207.2*g/mole   );
  G4Element* I  = new G4Element(name="Iodine",    symbol="I",  z=53, a=126.9045*g/mole);
  G4Element* Xe = new G4Element(name="Xenon",     symbol="Xe", z=54, a=131.293*g/mole );
  G4Element* Cs = new G4Element(name="Cesium",    symbol="Cs", z=55, a=132.91*g/mole  );
  
  G4Isotope* B10= new G4Isotope(name="B10", 5, 10, a=10.013*g/mole);
  G4Isotope* B11= new G4Isotope(name="B11", 5, 11, a=11.009*g/mole);
  G4Element* NewB=new G4Element(name="NewB", symbol="NewB", 2);
  NewB->AddIsotope(B10, 19.9*perCent);
  NewB->AddIsotope(B11, 80.1*perCent);


  G4Isotope* Gd152 = new G4Isotope("Gd152", 64, 152, 152.0*g/mole);
  G4Isotope* Gd154 = new G4Isotope("Gd154", 64, 154, 154.0*g/mole);
  G4Isotope* Gd155 = new G4Isotope("Gd155", 64, 155, 155.0*g/mole);
  G4Isotope* Gd156 = new G4Isotope("Gd156", 64, 156, 156.0*g/mole);
  G4Isotope* Gd157 = new G4Isotope("Gd157", 64, 157, 157.0*g/mole);
  G4Isotope* Gd158 = new G4Isotope("Gd158", 64, 158, 158.0*g/mole);
  G4Isotope* Gd160 = new G4Isotope("Gd160", 64, 160, 160.0*g/mole);

  G4Element* Gd = new G4Element("Gadolinium","Gd",7);
  Gd->AddIsotope(Gd152,  0.2*perCent);
  Gd->AddIsotope(Gd154,  2.2*perCent);
  Gd->AddIsotope(Gd155, 14.9*perCent);   //beware: it is abundance,
  Gd->AddIsotope(Gd156, 20.6*perCent);   //        not fractionMass
  Gd->AddIsotope(Gd157, 15.7*perCent);
  Gd->AddIsotope(Gd158, 24.7*perCent);
  Gd->AddIsotope(Gd160, 21.7*perCent);

  //---------------------//
  // Material Definition //
  //---------------------//
  
  // PC
  density = 0.882*g/cm3; //should be corrected
  fPC = new G4Material(name="PC", density, ncomponents=2);
  fPC->AddElement(C, natoms=9);
  fPC->AddElement(H, natoms=12);
  
  //TMB
  density = 0.932*g/cm3; 
  fTMB= new G4Material(name="TMB", density, ncomponents=4);
  fTMB->AddElement(H, natoms=9);
  fTMB->AddElement(O, natoms=3);
  fTMB->AddElement(C, natoms=3);
  fTMB->AddElement(NewB, natoms=1);
  //fTMB->AddElement(B, natoms=3);
  
  double fraction = DSStorage::Get()->GetTMBfraction();
  DSLog(trace) << "TMB/PC ratio: "<<  fraction << endlog;
  density = 0.882*(1.-fraction)*g/cm3+fraction*0.932*g/cm3 ;
  //density=0.91*g/cm3;
   fBoronScintillator = new G4Material(name="BoronScintillator", density, ncomponents=2);
   fBoronScintillator->AddMaterial(fTMB,fractionmass= fraction);
   fBoronScintillator->AddMaterial(fPC,fractionmass=(1.-fraction) );
   //fBoronScintillator->AddMaterial(fTMB,fractionmass=50*perCent);
   //fBoronScintillator->AddMaterial(fPC,fractionmass=50*perCent );
   
  
  // boron loaded scintillator
  //  fBoronScintillator = new G4Material(name="BoronScintillator",density=0.876*g/cm3, ncomponents=4);
  // fBoronScintillator = new G4Material(name="BoronScintillator",density=0.904*g/cm3, ncomponents=4);
  //fBoronScintillator->AddElement(C,fractionmass=0.623);
  //fBoronScintillator->AddElement(O,fractionmass=0.231);
  //fBoronScintillator->AddElement(H,fractionmass=0.094);
  //fBoronScintillator->AddElement(B,fractionmass=0.052);


  // Air
  density = 1.290*mg/cm3;
  fAir = new G4Material(name="Air", density, ncomponents=2);
  fAir->AddElement(N, fractionmass=0.7);
  fAir->AddElement(O, fractionmass=0.3);
  
  // Vacuum
  density = 1.29e-20*g/cm3;
  fVacuum = new G4Material(name="Vacuum", density, ncomponents=2);
  fVacuum->AddElement(N, fractionmass=0.7);
  fVacuum->AddElement(O, fractionmass=0.3);

  // PPO
  density = 1.000*g/cm3; //should be corrected
  fPPO = new G4Material(name="PPO", density, ncomponents=4);
  fPPO->AddElement(C, natoms=15);
  fPPO->AddElement(H, natoms=11);
  fPPO->AddElement(N, natoms=1);
  fPPO->AddElement(O, natoms=1);
  
  // DMP
  density = 1.000*g/cm3; //should be corrected
  fDMP = new G4Material(name="DMP", density, ncomponents=4);
  fDMP->AddElement(C, natoms=10);
  fDMP->AddElement(H, natoms=10);
  fDMP->AddElement(N, natoms=1);
  fDMP->AddElement(O, natoms=1);


  //Acrylic
  density = 1.19*g/cm3;
  fAcrylic = new G4Material(name="Acrylic",density, ncomponents=3, kStateSolid);
  fAcrylic->AddElement(C, 5);
  fAcrylic->AddElement(O, 2);
  fAcrylic->AddElement(H, 8);
 
  // Liquid Argon
  density      = 1.40*g/cm3;
  temperature  = 87*kelvin;
  fLiquidArgon = new G4Material(name="LiquidArgon",density,ncomponents=1,kStateLiquid, temperature);
  fLiquidArgon->AddElement(Ar,1.0);

  // Pseudo Argon
  density      = 1.40*g/cm3;
  temperature  = 87*kelvin;
  fPseudoArgon = new G4Material(name="PseudoArgon",density,ncomponents=1,kStateLiquid, temperature);
  fPseudoArgon->AddElement(Ar,1.0);

  // Non Scintillating Liquid Argon
  density      = 1.40*g/cm3;
  temperature  = 87*kelvin;
  fNSLiquidArgon = new G4Material(name="NSLiquidArgon",density,ncomponents=1,kStateLiquid, temperature);
  fNSLiquidArgon->AddElement(Ar,1.0);

  // Gaseous Argon
  density       = 5.4*mg/cm3;
  temperature   = 87.00*kelvin;
  pressure      = 0.93*atmosphere;
  fGaseousArgon = new G4Material(name="GaseousArgon",density,ncomponents=1,kStateGas,temperature,pressure);
  fGaseousArgon->AddElement(Ar,1.0);

  // Liquid Xenon
  density = 3.057*g/cm3;
  temperature = 165*kelvin; 
  fLiquidXenon = new G4Material(name="LiquidXenon",density,ncomponents=1,kStateLiquid,temperature);
  fLiquidXenon->AddElement(Xe,1.0);

  // Gaseous Xenon
  density = 5.89*mg/cm3;
  temperature=165*kelvin;
  pressure=0.93*atmosphere;
  fGaseousXenon = new G4Material(name="GaseousXenon",density,ncomponents=1,kStateGas,temperature,pressure);
  fGaseousXenon->AddElement(Xe,1.0);

   // Stainless steel (Medical Physics, Vol 25, No 10, Oct 1998)
  density = 8.290*g/cm3; 
  fSteel = new G4Material(name="Steel", density, ncomponents=5);
  fSteel->AddElement(Mn, 0.02);
  fSteel->AddElement(Si, 0.01);
  fSteel->AddElement(Cr, 0.19);
  fSteel->AddElement(Ni, 0.10);
  fSteel->AddElement(Fe, 0.68);


  // Stainless Steel
  density=7.7*g/cm3;
  fStainlessSteel = new G4Material(name="StainlessSteel",density,ncomponents=3);
  fStainlessSteel->AddElement(Cr,fractionmass= 0.20);
  fStainlessSteel->AddElement(Fe,fractionmass= 0.68);
  fStainlessSteel->AddElement(Ni,fractionmass= 0.12);

  // Grid Steel
  density=7.7*g/cm3;
  fGridSteel = new G4Material(name="GridSteel",density,ncomponents=3);
  fGridSteel->AddElement(Cr,fractionmass= 0.20);
  fGridSteel->AddElement(Fe,fractionmass= 0.68);
  fGridSteel->AddElement(Ni,fractionmass= 0.12);

  // Water
  density = 1.00*g/cm3;
  fWater = new G4Material(name="Water", density, ncomponents=2);
  fWater->AddElement(H,natoms=2);
  fWater->AddElement(O,natoms=1);

  // Metal Lead
  density = 11.34*g/cm3;
  fMetalLead =  new G4Material(name="MetalLead",density,ncomponents=1);
  fMetalLead->AddElement(Pb,1.0);

  //  Teflon
  density=2.165*g/cm3;
  fTeflon = new G4Material(name="Teflon", density=2.165*g/cm3, ncomponents=2);
  fTeflon->AddElement(Fl,fractionmass= 0.76);
  fTeflon->AddElement(C, fractionmass= 0.24);


  //  Rock: hep-ex/0312050v2 Wulandari et al.
  density = 2.71*g/cm3; //
  fRock = new G4Material(name="Rock", density, ncomponents=7);
  fRock->AddElement(C,  fractionmass=0.1188);
  fRock->AddElement(O,  fractionmass=0.4892); // <-- here I added 1.01% to obtain a fractionmass = 1
  fRock->AddElement(Mg, fractionmass=0.0558);
  fRock->AddElement(Al, fractionmass=0.0103);
  fRock->AddElement(Si, fractionmass=0.0127);
  fRock->AddElement(K,  fractionmass=0.0103);
  fRock->AddElement(Ca, fractionmass=0.3029);


  //  Norite (SNOlab rock) Density from Thomas Cement AB (Thomas Concrete Group)
  density=3.15*g/cm3;
  fNorite = new G4Material(name="Norite",density,ncomponents=12,kStateSolid);
  fNorite->AddElement(H, fractionmass=0.15*perCent);
  fNorite->AddElement(C, fractionmass=0.04*perCent);
  fNorite->AddElement(O, fractionmass=46.0*perCent);
  fNorite->AddElement(Na,fractionmass=2.2*perCent );
  fNorite->AddElement(Mn,fractionmass=3.3*perCent );
  fNorite->AddElement(Al,fractionmass=9.0*perCent );
  fNorite->AddElement(Si,fractionmass=26.2*perCent);
  fNorite->AddElement(K, fractionmass=1.2*perCent );
  fNorite->AddElement(Ca,fractionmass=5.2*perCent );
  fNorite->AddElement(Mn,fractionmass=0.1*perCent );
  fNorite->AddElement(Fe,fractionmass=6.2*perCent );
  fNorite->AddElement(Ti,fractionmass=0.5*perCent );


   // High Density PolyEthylene
  density=0.942*g/cm3;
  fHDPE = new G4Material(name="HDPE", density, ncomponents=2);
  fHDPE->AddElement(H,0.135);
  fHDPE->AddElement(C,0.865);



   // Fake Styrofoam
  density=0.0996996029*g/cm3;
  fStyrofoam = new G4Material(name="Styrofoam",density, ncomponents=2);
  fStyrofoam->AddElement(H,0.135);
  fStyrofoam->AddElement(C,0.865);


  // Metal Copper
  density=8.96*g/cm3;
  fMetalCopper = new G4Material(name="MetalCopper",density,ncomponents=1);
  fMetalCopper->AddElement(Cu,1.0);

   // Metal Titanium
  density=4.506*g/cm3;
  fMetalTitanium = new G4Material(name="MetalTitanium",density, ncomponents=1);
  fMetalTitanium->AddElement(Ti,1.0);
  //Ship Steel
  density=7.85*g/cm3;
  fShipSteel = new G4Material(name="ShipSteel",density,ncomponents=2);
  fShipSteel->AddElement(C, fractionmass=0.002);
  fShipSteel->AddElement(Fe,fractionmass=0.998);


  //TPB (Tetra-Phenyl Butadiene)
  //Density copied from PseudoCumene
  density=0.9*g/cm3;
  fTPB = new G4Material (name="TPB",density,ncomponents=2,kStateSolid);
  fTPB->AddElement(C,28);
  fTPB->AddElement(H,22);


  //ThreeMFoil (Polyethylene terephthalate)
  //http://multimedia.3m.com/mws/mediawebserver?mwsId=66666UuZjcFSLXTtM8T6oXTVEVuQEcuZgVs6EVs6E666666--
  //Density and composition from wiki
  density=1.4*g/cm3;
  fThreeMFoil = new G4Material(name="ThreeMFoil",density,ncomponents=3);
  fThreeMFoil->AddElement(C,10);
  fThreeMFoil->AddElement(H,8);
  fThreeMFoil->AddElement(O,4);

  //Fused Silica (same composition as quartz)
  //Density from Wikipedia
  density=2.203*g/cm3;
  fFusedSilica = new G4Material(name="FusedSilica",density,ncomponents=2,kStateSolid);
  fFusedSilica->AddElement(Si,1);
  fFusedSilica->AddElement(O,2);

  // Sodium Iodide
  density=3.67*g/cm3;
  fSodiumIodide = new G4Material(name="SodiumIodide",density,ncomponents=2);
  fSodiumIodide->AddElement(Na,natoms=1);
  fSodiumIodide->AddElement(I,natoms=1);


  //ITO
  //Indium corporation product data sheet
  density=7.2*g/cm3;
  fITO = new G4Material(name="ITO",density,ncomponents=3);
  fITO->AddElement(In,fractionmass=0.7444);
  fITO->AddElement(Sn,fractionmass=0.0788);
  fITO->AddElement(O,fractionmass=0.1768);

  // Liquid Nitrogen
  density = 0.807*g/cm3;
  fLiquidNitrogen = new G4Material(name="LiquidNitrogen",density,ncomponents=1,kStateLiquid,temperature = 77*kelvin);
  fLiquidNitrogen->AddElement(N,1.0);

  // Gaseous Nitrogen
  density = 0.000125*g/cm3;
  fGaseousNitrogen = new G4Material(name="GaseousNitrogen",density,ncomponents=1,kStateGas,temperature = 298*kelvin);
  fGaseousNitrogen->AddElement(N,1.0);

  // Kovar
  density=8.0*g/cm3;
  fKovar = new G4Material(name="Kovar",density,ncomponents=3);
  fKovar->AddElement(Fe,fractionmass=0.54);
  fKovar->AddElement(Ni,fractionmass=0.29);
  fKovar->AddElement(Co,fractionmass=0.17);


  //Bialkali (Wiki: Photocathode)
  density=3.39*g/cm3;
  fBialkali =  new G4Material(name="Bialkali",density,ncomponents=3);
  fBialkali->AddElement(Sb,fractionmass=0.34);
  fBialkali->AddElement(Rb,fractionmass=0.33);
  fBialkali->AddElement(Cs,fractionmass=0.33);

  //Kapton (Wiki: Kapton)
  density=1.0*g/cm3;
  fKapton = new G4Material(name="Kapton",density,ncomponents=4);
  fKapton->AddElement(C,22);
  fKapton->AddElement(N,2);
  fKapton->AddElement(O,5);
  fKapton->AddElement(H,10);


  //Concrete: hep-ex/0312050v2 Wulandari et al.
  density = 2.4*g/cm3; //
  fConcrete = new G4Material(name="Concrete", density, ncomponents=13);
  fConcrete->AddElement(H,  fractionmass=0.0089);
  fConcrete->AddElement(C,  fractionmass=0.0799);
  fConcrete->AddElement(O,  fractionmass=0.4971);// <-- here I add 1.28% to obtain a fractionmass = 1
  fConcrete->AddElement(Na, fractionmass=0.0006);
  fConcrete->AddElement(Mg, fractionmass=0.0085);
  fConcrete->AddElement(Al, fractionmass=0.0009);
  fConcrete->AddElement(Si, fractionmass=0.0386);
  fConcrete->AddElement(P,  fractionmass=0.0004);
  fConcrete->AddElement(S,  fractionmass=0.0016);
  fConcrete->AddElement(K,  fractionmass=0.0054);
  fConcrete->AddElement(Ca, fractionmass=0.3534);// <-- here I add 1.28% to obtain a fractionmass = 1
  fConcrete->AddElement(Ti, fractionmass=0.0004);
  fConcrete->AddElement(Fe, fractionmass=0.0043);

  // Borex materials

  // Aluminum
  density = 2.700*g/cm3;
  fAluminum = new G4Material(name="Aluminum", density, ncomponents=1);
  fAluminum->AddElement(Al, fractionmass=1.);  
  
  // Glass
  density = 2.2*g/cm3; //should be corrected
  fGlass = new G4Material(name="Glass", density, ncomponents=2);
  fGlass->AddElement(Si, natoms=1);
  fGlass->AddElement(O, natoms=2);
 
  // Nylon
  density = 1.1*g/cm3; //should be corrected
  fNylon = new G4Material(name="Nylon", density, ncomponents=3);
  fNylon->AddElement(H, fractionmass=0.08);
  fNylon->AddElement(C, fractionmass=0.60);
  fNylon->AddElement(O, fractionmass=0.32);

  // NylonL
  fNylonL = new G4Material(name="NylonL", density, ncomponents=3);
  fNylonL->AddElement(H, fractionmass=0.08);
  fNylonL->AddElement(C, fractionmass=0.60);
  fNylonL->AddElement(O, fractionmass=0.32);

  // Bialkali
  density = 4.000*g/cm3;  //should be corrected, but not important
  fOldBialkali = new G4Material(name="Bialkali", density, ncomponents=2);
  fOldBialkali->AddElement(K, fractionmass=0.5);
  fOldBialkali->AddElement(Cs, fractionmass=0.5);

  // Paint
  density = 4.000*g/cm3;  //should be corrected, but not important
  fPaint = new G4Material(name="Paint", density, ncomponents=2);
  fPaint->AddElement(K, fractionmass=0.5);
  fPaint->AddElement(Cs, fractionmass=0.5);

  // fBorexScintillator
  density = 0.8774*g/cm3; 
  fBorexScintillator = new G4Material(name="BorexScintillator", density, ncomponents=2);
  fBorexScintillator->AddMaterial(fPC, fractionmass=99.83*perCent);
  fBorexScintillator->AddMaterial(fPPO, fractionmass=0.17*perCent);
  
  // DMPBuffer
  //Density of the buffer reported by Frank Calaprice is by 0.2 kg/m^3 higher than the scintillator 
  //2.5 g/l DMP - October 2009     
  density = 0.8776*g/cm3;
  fDMPbuffer = new G4Material(name="DMPbuffer", density, ncomponents=2);
  fDMPbuffer->AddMaterial(fPC, fractionmass=99.66*perCent);
  fDMPbuffer->AddMaterial(fDMP, fractionmass=0.34*perCent);

  
  // Quartz
  density=2.200*g/cm3;
  fQuartz = new G4Material (name="Quartz", density, ncomponents=2);
  fQuartz->AddElement(Si, 1);
  fQuartz->AddElement(O , 2);
  
  //  Beryllium
  density = 1.848*g/cm3;
  fBe = new G4Material(name="Beryllium", density, ncomponents=2);
  fBe->AddElement(H, fractionmass=0.067);
  fBe->AddElement(C, fractionmass=0.933); 
  

 
  // Derlin
  density = 1.47*g/cm3;
  fDerlin = new G4Material(name="Derlin", density, ncomponents=2);
  fDerlin->AddElement(H, fractionmass=0.067);
  fDerlin->AddElement(C, fractionmass=0.933);
 
  // BlackHole
  density = 100*g/cm3;
  fBlackHole = new G4Material(name="BlackHole", density, ncomponents=1);
  fBlackHole->AddElement(C, 2);
   
  // Gadolinium Oxde
  density = 7.41*g/cm3; 
  fGdOxide = new G4Material(name="fGdOxide", density, ncomponents=2);
  fGdOxide->AddElement(Gd, natoms=2);
  fGdOxide->AddElement(O, natoms=3);

  // GdScintillator
  density = 0.8776*g/cm3;
  fGdScintillator = new G4Material(name="GdScintillator", density, ncomponents=2);
  fGdScintillator->AddMaterial(fBorexScintillator,  fractionmass=(877.6-1.)/877.6);
  fGdScintillator->AddMaterial(fGdOxide, fractionmass=1./877.6 );

  // GdWater
  density = 1.0*g/cm3;
  fGdWater = new G4Material(name="GdWater", density, ncomponents=2);
  fGdWater->AddMaterial(fWater,  fractionmass=(1000.-1.)/1000);
  fGdWater->AddMaterial(fGdOxide, fractionmass=1./1000. );

  // Metal Silicon
  density = 2.33*g/cm3;
  fMetalSilicon =  new G4Material(name="MetalSilicon",density,ncomponents=1);
  fMetalSilicon->AddElement(Si,1.0);

 

}


void DSMaterial::DefineProperties() {
  //---------------------------------------------------------------------------------
  // Scintillator
  //---------------------------------------------------------------------------------
  const G4double hc = 1239.84172; //eV*nm

  // local variables
  G4float   myvalue, myene;
  G4int   dim ;
  const G4int BScint_NUMENTRIES = 401;
  const G4int PPOScint_NUMENTRIES = 401;
  const G4int PPOABS_NUMENTRIES = 241;
  //const G4int PCABS_NUMENTRIES  = 401;
  //const G4int TMBABS_NUMENTRIES = 401;
  G4double BScint_Energy[BScint_NUMENTRIES];
  G4double BScint_SCINT[BScint_NUMENTRIES];
  G4double BScint_ABSL_Energy[BScint_NUMENTRIES];
  G4double BScint_ABSL[BScint_NUMENTRIES];
  G4double BScint_RINDEX[BScint_NUMENTRIES];
  G4double BScint_const_rind = 1.4; // Choosing 1.4, since it is between 1.3 (TMB) and 1.5 (PC)
  G4double ppo_Energy[PPOScint_NUMENTRIES];
  G4double ppo_SCINT[BScint_NUMENTRIES];
  G4double ppo_ABSL_Energy[BScint_NUMENTRIES];
  G4double ppo_ABSL[BScint_NUMENTRIES];

  // Fill in values for SS_RIND now, since it needs the same wavelengths as the scint
  G4double StainlessSteel_RIND[BScint_NUMENTRIES];
  G4double StainlessSteel_const_rind =  1.4;//1.56;
  
  // Load BScint properties, as measured by Aldo
  char bscintFileName[128] = "../data/detector/bscintSpectrum.dat";
  std::ifstream bscint_data(bscintFileName);
  G4int il = 0;
  if(bscint_data.is_open())
    {
      while(bscint_data >> BScint_Energy[BScint_NUMENTRIES-1-il] >> BScint_SCINT[BScint_NUMENTRIES-1-il])
	{
	  // Convert BScint_Energy values fromnm to eV
	  BScint_Energy[BScint_NUMENTRIES-1-il] = (hc/BScint_Energy[BScint_NUMENTRIES-1-il])*eV;
	  il++;
	}
      for(int j = 0; j < BScint_NUMENTRIES; j++)
	{
	  DSLog(debugging) << "BSCINT ENERGY = " << BScint_Energy[j]/eV << " eV\tSCINT = " << BScint_SCINT[j] << endlog;
	}
    }
  else
    DSLog(warning) << "ERROR: COULD NOT OPEN BSCINT SCINT FILE " << bscintFileName << endlog;
  if(!bscint_data.eof())
    DSLog(warning) << "ERROR: CLOSED BSCINT SCINT FILEBEFORE REACHING END" << endlog;
  bscint_data.close();

  char bscintAbslFileName[128] = "../data/detector/bscintAbsl.dat";
  std::ifstream bscintabsl_data(bscintAbslFileName);
  il = 0;
  if(bscintabsl_data.is_open())
    {
      while(bscintabsl_data >> BScint_ABSL_Energy[BScint_NUMENTRIES-1-il] >> BScint_ABSL[BScint_NUMENTRIES-1-il])
	{
	  // Convert BScint_Energy values fromnm to eV
	  BScint_ABSL_Energy[BScint_NUMENTRIES-1-il] = (hc/BScint_ABSL_Energy[BScint_NUMENTRIES-1-il])*eV;
	  BScint_ABSL[BScint_NUMENTRIES-1-il] *= m;
	  BScint_RINDEX[il]      = BScint_const_rind;
	  il++;
	}
      for(int j = 0; j < BScint_NUMENTRIES; j++)
	{
	  DSLog(debugging) << "BSCINT ENERGY = " << BScint_ABSL_Energy[j]/eV << " eV\tABSL = " << BScint_ABSL[j]/m << " m \tRINDEX = " << BScint_RINDEX[j] << endlog;
	}
    }
  else
    DSLog(warning) << "ERROR: COULD NOT OPEN BSCINT ABSL FILE " << bscintAbslFileName << endlog;
  if(!bscintabsl_data.eof())
    DSLog(warning) << "ERROR: CLOSED BSCINT ABSL FILE BEFORE REACHING END" << endlog;
  bscintabsl_data.close();

  // Load PPO Properties
  char ppoFileName[128] = "../data/detector/ppoEmissionSpectrum.dat";
  std::ifstream ppo_data(ppoFileName);
  il = 0;
  if(ppo_data.is_open())
    {
      while(ppo_data >> ppo_Energy[BScint_NUMENTRIES-1-il] >> ppo_SCINT[BScint_NUMENTRIES-1-il])
	{
	  // Convert ppo_Energy values from nm to eV
	  ppo_Energy[BScint_NUMENTRIES-1-il] = (hc/ppo_Energy[BScint_NUMENTRIES-1-il])*eV;
	  StainlessSteel_RIND[BScint_NUMENTRIES-1-il] = StainlessSteel_const_rind;
	  il++;
	}
    }
  else DSLog(warning) << "ERROR: COULD NOT OPEN FILE : " << ppoFileName << endlog;
  if(!ppo_data.eof())
    DSLog(warning) << "ERROR: CLOSED PPO DATA FILE BEFORE REACHING END" << endlog;
  ppo_data.close();
  
  char ppoAbsFileName[128] = "../data/detector/ppoAbsorptionLength.dat";
  std::ifstream ppoabs_data(ppoAbsFileName);
  il = 0;
  G4double wavelength, absl;
  if(ppoabs_data.is_open())
    {
      while((ppoabs_data >> wavelength >> absl) && (il < PPOABS_NUMENTRIES))
	{
	  // Factor of 1/2 because measurement was taken with 1.5 g/L rather than 3 g/L PPO
	  ppo_ABSL[BScint_NUMENTRIES-1-il] = absl/2.*cm;
	  ppo_ABSL_Energy[BScint_NUMENTRIES-1-il] = (hc/wavelength)*eV;
	  //BScint_RIND[BScint_NUMENTRIES-1-i] = BScint_const_rind;
	  il++;
	}
      while(il < BScint_NUMENTRIES)
	{
	  wavelength = 800.-250.*(il-PPOABS_NUMENTRIES)/(BScint_NUMENTRIES-PPOABS_NUMENTRIES);
	  ppo_ABSL[il-PPOABS_NUMENTRIES] = 28509.5*cm;
	  ppo_ABSL_Energy[il-PPOABS_NUMENTRIES] = hc/wavelength*eV;
	  //BScint_RIND[i-PPOABS_NUMENTRIES] = BScint_const_rind;
	  il++;
	}

/*      DSLog(debugging) << "***BSCINT RINDEX***" << endlog;
      for(int j=0; j < BScint_NUMENTRIES; j++)
	{
	  DSLog(debugging) << "WL = " << hc/ppo_ABSL_Energy[j] << " nm, RINDEX = " << BScint_RIND[j] << endlog;
	}
*/	
      ppo_ABSL_Energy[il-1] = 4.9*eV;
    }
  else DSLog(warning) << "ERROR: COULD NOT OPEN FILE : " << ppoAbsFileName << endlog;
  if(!ppoabs_data.eof())
    DSLog(warning) << "ERROR : CLOSED PPO ABS DATA FILE BEFORE REACHING END" << endlog;
  ppoabs_data.close();
  
  /*
  // Define scintillator absorption lengths
  char pcAbslFileName[128] = "../data/detector/pcAbsorptionLength.dat";
  std::ifstream pcabs_data(pcAbslFileName);
  char tmbAbslFileName[128] = "../data/detector/tmbAbsorptionLength.dat";
  std::ifstream tmbabs_data(tmbAbslFileName);
  i = 0;
  G4double pcWavelength, tmbWavelength;
  G4double pcAbsl, tmbAbsl;
  if(pcabs_data.is_open())
    {
      if(tmbabs_data.is_open())
	{
	  while((pcabs_data >> pcWavelength >> pcAbsl) &&
		(tmbabs_data >> tmbWavelength >> tmbAbsl) &&
		(i < PCABS_NUMENTRIES) &&
		(i < TMBABS_NUMENTRIES))
	    {
	      if(pcWavelength == tmbWavelength)
		{
		  BScint_ABSL_Energy[i] = (hc/pcWavelength)*eV;
		  BScint_ABSL[i]        = (2*pcAbsl*tmbAbsl/(pcAbsl+tmbAbsl))*m;
		  BScint_RINDEX[i]      = BScint_const_rind;
		}
	      i++;
	    }
	}
      else
	DSLog(warning) << "ERROR : COULD NOT OPEN TMB ABSL FILE : " << tmbAbslFileName << endlog;
    }
  else
    DSLog(warning) << "ERROR : COULD NOT OPEN PC ABSL FILE : " << pcAbslFileName << endlog;
  */


  //  G4double scintYield = 8868/MeV; // With POPOP
  //  G4double scintYield = 0.71*13249.8/MeV; // with bis-MSB
  G4double scintYield = 13249.8/MeV; // Our original guess
  scintYield *= DSParameters::Get()->GetVPmtMaxQe_adjusted();
  // Enter Material Properties into table
  G4MaterialPropertiesTable *BScint_MPT = new G4MaterialPropertiesTable();
  BScint_MPT->AddProperty("FASTCOMPONENT", BScint_Energy, BScint_SCINT, BScint_NUMENTRIES)->SetSpline(true);
  BScint_MPT->AddProperty("SLOWCOMPONENT", BScint_Energy, BScint_SCINT, BScint_NUMENTRIES)->SetSpline(true);
  BScint_MPT->AddProperty("RINDEX",        BScint_Energy, BScint_RINDEX,BScint_NUMENTRIES )->SetSpline(true);
  BScint_MPT->AddProperty("ABSLENGTH",     BScint_ABSL_Energy, BScint_ABSL, BScint_NUMENTRIES)->SetSpline(true);
  BScint_MPT->AddConstProperty("SCINTILLATIONYIELD", scintYield);
  BScint_MPT->AddConstProperty("YIELDRATIO", 0.93);
  BScint_MPT->AddConstProperty("RESOLUTIONSCALE",1.0);
  //  BScint_MPT->AddConstProperty("FASTTIMECONSTANT", 3.57*ns);
  BScint_MPT->AddConstProperty("FASTTIMECONSTANT", 2.0*ns);
  BScint_MPT->AddConstProperty("SLOWTIMECONSTANT", 17.61*ns);
  BScint_MPT->AddProperty("WLSCOMPONENT", ppo_Energy, ppo_SCINT, BScint_NUMENTRIES);
  BScint_MPT->AddProperty("WLSABSLENGTH", ppo_ABSL_Energy, ppo_ABSL, BScint_NUMENTRIES);
  BScint_MPT->AddConstProperty("WLSTIMECONSTANT",1.8*ns);
  //  BScint_MPT->AddConstProperty("WLSEFFICIENCY",0.82);
  BScint_MPT->AddConstProperty("WLSEFFICIENCY",1.);
  fBoronScintillator->SetMaterialPropertiesTable(BScint_MPT);
  fGdScintillator->SetMaterialPropertiesTable(BScint_MPT);
  // Value of Birks's constant from Richard's thesis
  //  G4double pcBirk = 0.01154*cm/MeV;
  // Value from Borexino
  G4double pcBirk = 0.0109*cm/MeV;
  //  G4double pcBirk = 0.0059*cm/MeV;
  fBoronScintillator->GetIonisation()->SetBirksConstant(pcBirk);
  fGdScintillator->GetIonisation()->SetBirksConstant(pcBirk);
  
  //---------------------------------------------------------------------------------
  // Stainless Steel
  //---------------------------------------------------------------------------------
  G4double mySSEnergy[2] = {0.1*eV, 10*eV};
  G4double mySSABSL[2]   = {0.1*nm, 0.1*nm};
  //G4double mySSRI[2]     = {1.47, 1.47};
  G4MaterialPropertiesTable *ss_MPT = new G4MaterialPropertiesTable();
  //ss_MPT->AddProperty("RINDEX",mySSEnergy,mySSRI,2)->SetSpline(true); 
  ss_MPT->AddProperty("ABSLENGTH",mySSEnergy,mySSABSL,2);
  fStainlessSteel->SetMaterialPropertiesTable(ss_MPT);
  
  //---------------------------------------------------------------------------------
  // Gaseous Nitrogen
  //---------------------------------------------------------------------------------
  G4double myGN2Energy[2] = {0.1*eV, 10.*eV};
  G4double myGN2ABSL[2]   = {100.*m, 100.*nm};
  G4double myGN2RI[2]     = {1., 1.};
  G4MaterialPropertiesTable *gN2_MPT = new G4MaterialPropertiesTable();
  gN2_MPT->AddProperty("RINDEX",myGN2Energy,myGN2RI,2); 
  gN2_MPT->AddProperty("ABSLENGTH",myGN2Energy,myGN2ABSL,2);
  fGaseousNitrogen->SetMaterialPropertiesTable(gN2_MPT);

  //---------------------------------------------------------------------------------
  // Fused Silica
  //---------------------------------------------------------------------------------
  G4MaterialPropertiesTable *myFusedSilica = new G4MaterialPropertiesTable();
  G4double myFSEnergy[4],  myFSRI[4], myFSAbs[4], myRLength[4];  
  myFSEnergy[0] = 1.0*eV  ;   myFSEnergy[1] = 8.0*eV ; 
  myFSRI[0] = DSParameters::Get()->GetFusedSilicaVisRind();  
  myFSAbs[0] = DSParameters::Get()->GetFusedSilicaVisAbs()*m; 
  myRLength[0] = DSParameters::Get()->GetFSilicaRaylVisLength()*m;
  myFSRI[1] = myFSRI[0];  myFSAbs[1] =myFSAbs[0] ; myRLength[1] =myRLength[0];
  myFSEnergy[2] = 8.3*eV ;  myFSEnergy[3] = 20.0*eV ; 
  myFSRI[2] = DSParameters::Get()->GetFusedSilicaUVRind();  
  myFSAbs[2] = DSParameters::Get()->GetFusedSilicaUVAbs()*m;;
  myRLength[2] = DSParameters::Get()->GetFSilicaRaylUVLength()*m;
  myFSRI[3] = myFSRI[2];  myFSAbs[3] = myFSAbs[2];  myRLength[3] =myRLength[2];
  myFusedSilica->AddProperty("RINDEX",    myFSEnergy, myFSRI,   4);
  myFusedSilica->AddProperty("ABSLENGTH", myFSEnergy, myFSAbs , 4);
  myFusedSilica->AddProperty("RAYLEIGH", myFSEnergy, myRLength , 4);
  fFusedSilica->SetMaterialPropertiesTable(myFusedSilica);


  //---------------------------------------------------------------------------------
  // Teflon 
  //---------------------------------------------------------------------------------

  G4MaterialPropertiesTable *myTeflon = new G4MaterialPropertiesTable;    
  G4double myTeflonEnergy[3],  myTeflonRI[3], myTeflonAbs[3];  
  myTeflonEnergy[0] = 1.0*eV  ; myTeflonRI[0] = 1.4;  myTeflonAbs[0] = .1*mm;
  myTeflonEnergy[1] = 5.0*eV ;  myTeflonRI[1] = 1.4;  myTeflonAbs[1] = .1*mm;
  myTeflonEnergy[2] = 10.0*eV ; myTeflonRI[2] = 1.4;  myTeflonAbs[2] = .1*mm;
  //myTeflon->AddProperty("RINDEX", myTeflonEnergy,    myTeflonRI, 3 );
  myTeflon->AddProperty("ABSLENGTH", myTeflonEnergy, myTeflonAbs , 3);
  fTeflon->SetMaterialPropertiesTable(myTeflon);

  //---------------------------------------------------------------------------------
  // Bialkali
  //---------------------------------------------------------------------------------
  // this is used instead of Bialkali for the cathode as fake material to stop photons
  // FusedSilica can not be used, because of the window disk 

  G4MaterialPropertiesTable *myBialkali = new G4MaterialPropertiesTable();
  G4double myBialkaliEnergy[4],   myBialkaliRI[4], myBialkaliAbs[4];  
  myBialkaliEnergy[0] = 1.0*eV  ; 
  myBialkaliRI[0] = DSParameters::Get()->GetPhotocathodeVisRind();    myBialkaliAbs[0] = 0.1*nm;
  myBialkaliEnergy[1] = 8.0*eV ;  
  myBialkaliRI[1] = DSParameters::Get()->GetPhotocathodeVisRind();  ;  myBialkaliAbs[1] = 0.1*nm;
  myBialkaliEnergy[2] = 8.3*eV ;  
  myBialkaliRI[2] = DSParameters::Get()->GetPhotocathodeUVRind();  myBialkaliAbs[2] = 0.1*nm;
  myBialkaliEnergy[3] = 20.0*eV ; 
  myBialkaliRI[3] = DSParameters::Get()->GetPhotocathodeUVRind();  myBialkaliAbs[3] = 0.1*nm;
  myBialkali->AddProperty("RINDEX",    myBialkaliEnergy, myBialkaliRI,  4);
  myBialkali->AddProperty("ABSLENGTH", myBialkaliEnergy, myBialkaliAbs, 4);
  fBialkali->SetMaterialPropertiesTable(myBialkali);

  //---------------------------------------------------------------------------------
  // Metal Silicon
  //---------------------------------------------------------------------------------
  // same goal as bialkali
  G4MaterialPropertiesTable *myMetalSilicon= new G4MaterialPropertiesTable();
  //http://refractiveindex.info/?shelf=main&book=Si&page=Pierce 
  //myMetalSilicon->AddProperty("IMAGINARYRINDEX",  mySiEnergy  ,mySiImRind ,   6);
  myMetalSilicon->AddProperty("RINDEX",    myBialkaliEnergy, myBialkaliRI,  4);
  myMetalSilicon->AddProperty("ABSLENGTH", myBialkaliEnergy, myBialkaliAbs, 4);
  fMetalSilicon->SetMaterialPropertiesTable(myMetalSilicon);
  
  //---------------------------------------------------------------------------------
  // Liquid Argon
  //---------------------------------------------------------------------------------
  G4MaterialPropertiesTable *myLiquidArgon = new G4MaterialPropertiesTable();
  dim = 0;
  G4double LArRefIndex[1000], LArRefIndexEne[1000],LArRayLength[1000],LArAbsLength[1000] ;  
  for(int ij=0;ij<55; ij++) {
    LArRefIndexEne[ij] = (1.0 + ij*0.2)*eV;
    double lambda     = h_Planck*c_light/LArRefIndexEne[ij];
    LArRefIndex[ij]    = GetLArRefIndex(lambda/nm);
    LArRayLength[ij]   = GetLArRayLength(lambda/nm)*DSParameters::Get()->GetLArRayleighScale();
    // L700 why no to modify LAr absorption length for VUV light? 
    if (lambda < 150*nm )  LArAbsLength[ij]   = DSParameters::Get()->GetLiquidArgonUVAbs()*m;
    if (lambda > 150*nm )  LArAbsLength[ij]   = DSParameters::Get()->GetLiquidArgonVisAbs()*m;
    cout << lambda/nm <<" " << LArRayLength[ij]/m  << endl ; 
    dim = ij;
  }    
  myLiquidArgon->AddProperty("RINDEX",   LArRefIndexEne , LArRefIndex,  dim+1);
  myLiquidArgon->AddProperty("RAYLEIGH", LArRefIndexEne , LArRayLength, dim+1);
  myLiquidArgon->AddProperty("ABSLENGTH",LArRefIndexEne , LArAbsLength, dim+1);  
  
  //field for NEST
  myLiquidArgon->AddConstProperty( "ELECTRICFIELD", 1000*volt/cm ); //for missed nooks and crannies
  myLiquidArgon->AddConstProperty( "TOTALNUM_INT_SITES", -1 );

  fLiquidArgon->SetMaterialPropertiesTable(myLiquidArgon);
  fNSLiquidArgon->SetMaterialPropertiesTable(myLiquidArgon);

  //---------------------------------------------------------------------------------
  // Pseudo Argon
  //---------------------------------------------------------------------------------
  G4MaterialPropertiesTable *myPseudoArgon = new G4MaterialPropertiesTable();
  dim = 0;
  G4double PArRefIndex[1000], PArRefIndexEne[1000],PArRayLength[1000],PArAbsLength[1000] ;
  for(int ij=0;ij<55; ij++) {
    PArRefIndexEne[ij] = (1.0 + ij*0.2)*eV;
    double lambda     = h_Planck*c_light/PArRefIndexEne[ij];
    PArRefIndex[ij]    = DSParameters::Get()->GetPArRind();
    PArRayLength[ij]   = GetLArRayLength(lambda/nm)*DSParameters::Get()->GetLArRayleighScale();
    // L700 why no to modify LAr absorption length for VUV light? 
    if (lambda < 150*nm )  PArAbsLength[ij]   = DSParameters::Get()->GetLiquidArgonUVAbs()*m;
    if (lambda > 150*nm )  PArAbsLength[ij]   = DSParameters::Get()->GetLiquidArgonVisAbs()*m;
    dim = ij;
  }
  myPseudoArgon->AddProperty("RINDEX",   PArRefIndexEne , PArRefIndex,  dim+1);
  myPseudoArgon->AddProperty("RAYLEIGH", PArRefIndexEne , PArRayLength, dim+1);
  myPseudoArgon->AddProperty("ABSLENGTH",PArRefIndexEne , PArAbsLength, dim+1);

  //field for NEST
  myPseudoArgon->AddConstProperty( "ELECTRICFIELD", 1000*volt/cm ); //for missed nooks and crannies
  myPseudoArgon->AddConstProperty( "TOTALNUM_INT_SITES", -1 );

  fPseudoArgon->SetMaterialPropertiesTable(myPseudoArgon);

  //---------------------------------------------------------------------------------
  // Gaseous Argon
  //---------------------------------------------------------------------------------
   
  G4MaterialPropertiesTable *myGaseousArgon = new G4MaterialPropertiesTable();
  G4double GArRefIndex[1000], GArRefIndexEne[1000], GArRayLength[1000], GArAbsLength[1000];
  dim = 0;   
  for(int ij=0;ij<55; ij++) {
    GArRefIndexEne[ij] = (1.4 + ij*0.2)*eV;
    double lambda     = h_Planck*c_light/GArRefIndexEne[ij];
//    GArRefIndex[ij]    = GetGArRefIndex(lambda/nm)*DSParameters::Get()->GetGArRindexScale();
    GArRefIndex[ij]    = GetGArRefIndex(lambda/nm);
    GArRayLength[ij]   = GetGArRayLength(lambda/nm);
    if (lambda < 150*nm )  GArAbsLength[ij]   = DSParameters::Get()->GetGaseousArgonUVAbs()*m;
    if (lambda > 150*nm )  GArAbsLength[ij]   = DSParameters::Get()->GetGaseousArgonVisAbs()*m;
    dim = ij;
  }
  
  myGaseousArgon->AddProperty("RINDEX",    GArRefIndexEne, GArRefIndex,  dim+1); 
  myGaseousArgon->AddProperty("RAYLEIGH",  GArRefIndexEne, GArRayLength, dim+1); 
  myGaseousArgon->AddProperty("ABSLENGTH", GArRefIndexEne, GArAbsLength, dim+1); 

 //field for NEST

  myGaseousArgon->AddConstProperty( "ELECTRICFIELD", 3000*volt/cm ); //for missed nooks and crannies
  myGaseousArgon->AddConstProperty( "TOTALNUM_INT_SITES", -1 );

  fGaseousArgon->SetMaterialPropertiesTable(myGaseousArgon);

  //---------------------------------------------------------------------------------
  // TPB
  //---------------------------------------------------------------------------------
  
  G4MaterialPropertiesTable *myTPB = new G4MaterialPropertiesTable();
  
  G4double TPB_ENE[1000], TPB_EMISSION_VAL[1000], TPB_ABSORPTION_VAL[1000];
  dim = 0 ;
  G4double TPBNorma = 0 ;
  ifstream ftpb_emission("../data/detector/tpb_emission_spectrum.dat");
  while(!ftpb_emission.eof()) {
    ftpb_emission >> myene >> myvalue ;
    if(ftpb_emission.eof()) break;
    TPB_ENE[dim] = myene*eV; 
    if(hc/ (TPB_ENE[dim]/nm) > 600) {
      TPB_EMISSION_VAL[dim] = 0 ;
    } else  {
      TPB_EMISSION_VAL[dim] = myvalue;
      TPBNorma += myvalue ;
    }
    dim++;
  }  

  G4double  TPB_WLS_ABSORPTION_VAL[1000], TPB_RINDEX[1000]; 
  dim = 0 ;
  ifstream ftpb_absorption("../data/detector/tpb_absorption_length.dat");
  while(!ftpb_absorption.eof()) {
    ftpb_absorption >> myene >> myvalue ;
    if(ftpb_absorption.eof()) break;
    TPB_WLS_ABSORPTION_VAL[dim] = myvalue*m*DSParameters::Get()->GetWLSAbsorptionFactor();
    if (myene > 8.3)  {
      TPB_ABSORPTION_VAL[dim]  = DSParameters::Get()->GetTPBUVAbs()*m;
      TPB_RINDEX[dim]          = DSParameters::Get()->GetTPBUVRind();
     }
    else {
      TPB_ABSORPTION_VAL[dim]  = DSParameters::Get()->GetTPBVisAbs()*m;
      TPB_RINDEX[dim]          = DSParameters::Get()->GetTPBVisRind();
    }
    TPB_EMISSION_VAL[dim]       /= TPBNorma;
    
    dim++;
  }

  myTPB->AddProperty("WLSCOMPONENT",TPB_ENE, TPB_EMISSION_VAL,       dim);
  myTPB->AddProperty("WLSABSLENGTH",TPB_ENE, TPB_WLS_ABSORPTION_VAL, dim);
  myTPB->AddProperty("ABSLENGTH",   TPB_ENE, TPB_ABSORPTION_VAL,     dim);
  myTPB->AddProperty("RINDEX",      TPB_ENE, TPB_RINDEX,             dim);
  myTPB->AddConstProperty("WLSMEANNUMBERPHOTONS",DSParameters::Get()->GetWLSMeanNumberPhotons());
  myTPB->AddConstProperty("WLSTIMECONSTANT",DSParameters::Get()->GetWLSTimeConstant_ns()*ns);
  myTPB->AddConstProperty("WLSEFFICIENCY",DSParameters::Get()->GetWLSEfficiency());
  fTPB->SetMaterialPropertiesTable(myTPB);
  

  //---------------------------------------------------------------------------------
  // GridSteel
  //---------------------------------------------------------------------------------
  G4MaterialPropertiesTable *myGridSteel = new G4MaterialPropertiesTable();
  //G4double myGridSteelEnergy[4],  myGridSteelAbs[4], myGridRefIndex[4];  
  dim = 0;
  G4double GridRefIndex[1000], GridRefIndexEne[1000],GridAbsLength[1000] ;  
  for(int ij=0;ij<55; ij++) {
    GridRefIndexEne[ij] = (1.0 + ij*0.2)*eV;
    double lambda     = h_Planck*c_light/GridRefIndexEne[ij];
    GridRefIndex[ij]   = GetLArRefIndex(lambda/nm)*DSParameters::Get()->GetGridSteelRindScale();
    if (lambda < 150*nm ) GridAbsLength[ij]  = DSParameters::Get()->GetGridSteelUVAbs()*m;
      else             GridAbsLength[ij]  = DSParameters::Get()->GetGridSteelVisAbs()*m;
    dim = ij;
  }    

  myGridSteel->AddProperty("RINDEX",   GridRefIndexEne , GridRefIndex,  dim+1);
  myGridSteel->AddProperty("ABSLENGTH",GridRefIndexEne , GridAbsLength, dim+1);  
  fGridSteel->SetMaterialPropertiesTable(myGridSteel);


  //---------------------------------------------------------------------------------
  // Kovar
  //---------------------------------------------------------------------------------
  G4MaterialPropertiesTable *myKovar = new G4MaterialPropertiesTable();
  G4double myKovarEnergy[3],  myKovarAbs[3];  
  myKovarEnergy[0] = 0.1*eV  ; myKovarAbs[0] = 0.01*nm;
  myKovarEnergy[1] = 5.0*eV ;  myKovarAbs[1] = 0.01*nm;
  myKovarEnergy[2] = 20.0*eV ; myKovarAbs[2] = 0.01*nm;
  myKovar->AddProperty("ABSLENGTH", myKovarEnergy, myKovarAbs , 3);
  fKovar->SetMaterialPropertiesTable(myKovar);


  //---------------------------------------------------------------------------------
  // Black Hole
  //---------------------------------------------------------------------------------
  G4MaterialPropertiesTable *myBlackHole = new G4MaterialPropertiesTable();
  G4double myBHEnergy[3],  myBHRI[3], myBHAbs[3];  
  myBHEnergy[0] = 0.1*eV  ; myBHRI[0] = 1.0;  myBHAbs[0] = 0.001*mm;
  myBHEnergy[1] = 5.0*eV ;  myBHRI[1] = 1.0;  myBHAbs[1] = 0.001*mm;
  myBHEnergy[2] = 20.0*eV ; myBHRI[2] = 1.0;  myBHAbs[2] = 0.001*mm;
  myBlackHole->AddProperty("RINDEX",    myBHEnergy, myBHRI,   3);
  myBlackHole->AddProperty("ABSLENGTH", myBHEnergy, myBHAbs , 3);
  fBlackHole->SetMaterialPropertiesTable(myBlackHole);


}
//---------------------------------------------------------------------------------------
void DSMaterial::DefineBoundaries() {
      // Load properties into material properties table
    G4int nEntries = DSParameters::Get()->GetLumirrorNumEntries();
    fLumirrorMPT = new G4MaterialPropertiesTable();
    fLumirrorMPT->AddProperty("SPECULARLOBECONSTANT",
			      &DSParameters::Get()->GetLumirrorEnergy()[0],
			      &DSParameters::Get()->GetLumirrorSpecularLobe()[0],
			      nEntries);
    fLumirrorMPT->AddProperty("SPECULARSPIKECONSTANT",
			      &DSParameters::Get()->GetLumirrorEnergy()[0],
			      &DSParameters::Get()->GetLumirrorSpecularSpike()[0],
			      nEntries);
    fLumirrorMPT->AddProperty("BACKSCATTERCONSTANT",
			      &DSParameters::Get()->GetLumirrorEnergy()[0],
			      &DSParameters::Get()->GetLumirrorBackscatter()[0],
			      nEntries);
    fLumirrorMPT->AddProperty("REFLECTIVITY",
			      &DSParameters::Get()->GetLumirrorEnergy()[0],
			      &DSParameters::Get()->GetLumirrorReflectance()[0],
			      nEntries)->SetSpline(true);
    fLumirrorMPT->AddProperty("EFFICIENCY",
			      &DSParameters::Get()->GetLumirrorEnergy()[0],
			      &DSParameters::Get()->GetLumirrorEfficiency()[0],
			      nEntries);
    fLumirrorMPT->AddProperty("RINDEX",
			      &DSParameters::Get()->GetLumirrorEnergy()[0],
			      &DSParameters::Get()->GetLumirrorRindex()[0],
			      nEntries);
    
    fElectropolishedStainlessSteelMPT = new G4MaterialPropertiesTable();
    fElectropolishedStainlessSteelMPT->AddProperty("SPECULARLOBECONSTANT",
						   &DSParameters::Get()->GetElectropolishedSSEnergy()[0],
						   &DSParameters::Get()->GetElectropolishedSSSpecularLobe()[0],
						   nEntries);
    fElectropolishedStainlessSteelMPT->AddProperty("SPECULARSPIKECONSTANT",
						   &DSParameters::Get()->GetElectropolishedSSEnergy()[0],
						   &DSParameters::Get()->GetElectropolishedSSSpecularSpike()[0],
						   nEntries);
    fElectropolishedStainlessSteelMPT->AddProperty("BACKSCATTERCONSTANT",
						   &DSParameters::Get()->GetElectropolishedSSEnergy()[0],
						   &DSParameters::Get()->GetElectropolishedSSBackscatter()[0],
						   nEntries);
    fElectropolishedStainlessSteelMPT->AddProperty("REFLECTIVITY",
						   &DSParameters::Get()->GetElectropolishedSSEnergy()[0],
						   &DSParameters::Get()->GetElectropolishedSSReflectance()[0],
						   nEntries)->SetSpline(true);
    fElectropolishedStainlessSteelMPT->AddProperty("EFFICIENCY",
						   &DSParameters::Get()->GetElectropolishedSSEnergy()[0],
						   &DSParameters::Get()->GetElectropolishedSSEfficiency()[0],
						   nEntries);
    fElectropolishedStainlessSteelMPT->AddProperty("RINDEX",
						   &DSParameters::Get()->GetElectropolishedSSEnergy()[0],
						   &DSParameters::Get()->GetElectropolishedSSRindex()[0],
						   nEntries);

    fUntreatedStainlessSteelMPT = new G4MaterialPropertiesTable();
    fUntreatedStainlessSteelMPT->AddProperty("SPECULARLOBECONSTANT",
					     &DSParameters::Get()->GetUntreatedSSEnergy()[0],
					     &DSParameters::Get()->GetUntreatedSSSpecularLobe()[0],
					     nEntries);
    fUntreatedStainlessSteelMPT->AddProperty("SPECULARSPIKECONSTANT",
					     &DSParameters::Get()->GetUntreatedSSEnergy()[0],
					     &DSParameters::Get()->GetUntreatedSSSpecularSpike()[0],
					     nEntries);
    fUntreatedStainlessSteelMPT->AddProperty("BACKSCATTERCONSTANT",
					     &DSParameters::Get()->GetUntreatedSSEnergy()[0],
					     &DSParameters::Get()->GetUntreatedSSBackscatter()[0],
					     nEntries);
    fUntreatedStainlessSteelMPT->AddProperty("REFLECTIVITY",
					     &DSParameters::Get()->GetUntreatedSSEnergy()[0],
					     &DSParameters::Get()->GetUntreatedSSReflectance()[0],
					     nEntries)->SetSpline(true);
    fUntreatedStainlessSteelMPT->AddProperty("EFFICIENCY",
					     &DSParameters::Get()->GetUntreatedSSEnergy()[0],
					     &DSParameters::Get()->GetUntreatedSSEfficiency()[0],
					     nEntries);
    fUntreatedStainlessSteelMPT->AddProperty("RINDEX",
					     &DSParameters::Get()->GetUntreatedSSEnergy()[0],
					     &DSParameters::Get()->GetUntreatedSSRindex()[0],
					     nEntries);

    fAluminumFoilMPT = new G4MaterialPropertiesTable();
    fAluminumFoilMPT->AddProperty("SPECULARLOBECONSTANT",
				  &DSParameters::Get()->GetAluminumFoilEnergy()[0],
				  &DSParameters::Get()->GetAluminumFoilSpecularLobe()[0],
				  nEntries);
    fAluminumFoilMPT->AddProperty("SPECULARSPIKECONSTANT",
				  &DSParameters::Get()->GetAluminumFoilEnergy()[0],
				  &DSParameters::Get()->GetAluminumFoilSpecularSpike()[0],
				  nEntries);
    fAluminumFoilMPT->AddProperty("BACKSCATTERCONSTANT",
				  &DSParameters::Get()->GetAluminumFoilEnergy()[0],
				  &DSParameters::Get()->GetAluminumFoilBackscatter()[0],
				  nEntries);
    fAluminumFoilMPT->AddProperty("REFLECTIVITY",
				  &DSParameters::Get()->GetAluminumFoilEnergy()[0],
				  &DSParameters::Get()->GetAluminumFoilReflectance()[0],
				  nEntries);
    fAluminumFoilMPT->AddProperty("EFFICIENCY",
				  &DSParameters::Get()->GetAluminumFoilEnergy()[0],
				  &DSParameters::Get()->GetAluminumFoilEfficiency()[0],
				  nEntries);
    fAluminumFoilMPT->AddProperty("RINDEX",
				  &DSParameters::Get()->GetAluminumFoilEnergy()[0],
				  &DSParameters::Get()->GetAluminumFoilRindex()[0],
				  nEntries);

    fVPhotocathodeMPT = new G4MaterialPropertiesTable();
    fVPhotocathodeMPT->AddConstProperty("DOTRANSMISSION",1);
    fVPhotocathodeMPT->AddProperty("SPECULARLOBECONSTANT",
				   &DSParameters::Get()->GetVPhotocathodeEnergy()[0],
				   &DSParameters::Get()->GetVPhotocathodeSpecularLobe()[0],
				   DSParameters::Get()->GetVQENumEntries());
    fVPhotocathodeMPT->AddProperty("SPECULARSPIKECONSTANT",
				   &DSParameters::Get()->GetVPhotocathodeEnergy()[0],
				   &DSParameters::Get()->GetVPhotocathodeSpecularSpike()[0],
				   DSParameters::Get()->GetVQENumEntries());
    fVPhotocathodeMPT->AddProperty("BACKSCATTERCONSTANT",
				   &DSParameters::Get()->GetVPhotocathodeEnergy()[0],
				   &DSParameters::Get()->GetVPhotocathodeBackscatter()[0],
				   DSParameters::Get()->GetVQENumEntries());
    fVPhotocathodeMPT->AddProperty("REFLECTIVITY",
				   &DSParameters::Get()->GetVPhotocathodeEnergy()[0],
				   &DSParameters::Get()->GetVPhotocathodeReflectance()[0],
				   DSParameters::Get()->GetVQENumEntries());
    fVPhotocathodeMPT->AddProperty("EFFICIENCY",
				   &DSParameters::Get()->GetVPhotocathodeEnergy()[0],
				   &DSParameters::Get()->GetVPhotocathodeEfficiency()[0],
				   DSParameters::Get()->GetVQENumEntries());
    fVPhotocathodeMPT->AddProperty("RINDEX",
				   &DSParameters::Get()->GetVPhotocathodeEnergy()[0],
				   &DSParameters::Get()->GetVPhotocathodeRindex()[0],
				   DSParameters::Get()->GetVQENumEntries());

    fPMTBackMPT = new G4MaterialPropertiesTable();
    fPMTBackMPT->AddProperty("SPECULARLOBECONSTANT",
			     &DSParameters::Get()->GetPMTBackEnergy()[0],
			     &DSParameters::Get()->GetPMTBackSpecularLobe()[0],
			     nEntries);
    fPMTBackMPT->AddProperty("SPECULARSPIKECONSTANT",
			     &DSParameters::Get()->GetPMTBackEnergy()[0],
			     &DSParameters::Get()->GetPMTBackSpecularSpike()[0],
			     nEntries);
    fPMTBackMPT->AddProperty("BACKSCATTERCONSTANT",
			     &DSParameters::Get()->GetPMTBackEnergy()[0],
			     &DSParameters::Get()->GetPMTBackBackscatter()[0],
			     nEntries);
    fPMTBackMPT->AddProperty("REFLECTIVITY",
			     &DSParameters::Get()->GetPMTBackEnergy()[0],
			     &DSParameters::Get()->GetPMTBackReflectance()[0],
			     nEntries);
    fPMTBackMPT->AddProperty("EFFICIENCY",
			     &DSParameters::Get()->GetPMTBackEnergy()[0],
			     &DSParameters::Get()->GetPMTBackEfficiency()[0],
			     nEntries);
    fPMTBackMPT->AddProperty("RINDEX",
			     &DSParameters::Get()->GetPMTBackEnergy()[0],
			     &DSParameters::Get()->GetPMTBackRindex()[0],
			     nEntries);


    fPhotocathodeMPT = new G4MaterialPropertiesTable();
    fPhotocathodeMPT->AddProperty("SPECULARLOBECONSTANT",
				  &DSParameters::Get()->GetPhotocathodeEnergy()[0],
				  &DSParameters::Get()->GetPhotocathodeSpecularLobe()[0],
				  DSParameters::Get()->GetQENumEntries());
    fPhotocathodeMPT->AddProperty("SPECULARSPIKECONSTANT",
				  &DSParameters::Get()->GetPhotocathodeEnergy()[0],
				  &DSParameters::Get()->GetPhotocathodeSpecularSpike()[0],
				  DSParameters::Get()->GetQENumEntries());
    fPhotocathodeMPT->AddProperty("BACKSCATTERCONSTANT",
				  &DSParameters::Get()->GetPhotocathodeEnergy()[0],
				  &DSParameters::Get()->GetPhotocathodeBackscatter()[0],
				  DSParameters::Get()->GetQENumEntries());
    fPhotocathodeMPT->AddProperty("REFLECTIVITY",
				  &DSParameters::Get()->GetPhotocathodeEnergy()[0],
				  &DSParameters::Get()->GetPhotocathodeReflectance()[0],
				  DSParameters::Get()->GetQENumEntries());
    fPhotocathodeMPT->AddProperty("EFFICIENCY",
				  &DSParameters::Get()->GetPhotocathodeEnergy()[0],
				  &DSParameters::Get()->GetPhotocathodeEfficiency()[0],
				  DSParameters::Get()->GetQENumEntries());
    fPhotocathodeMPT->AddProperty("RINDEX",
				  &DSParameters::Get()->GetPhotocathodeEnergy()[0],
				  &DSParameters::Get()->GetPhotocathodeRindex()[0],
				  DSParameters::Get()->GetQENumEntries());

}

//---------------------------------------------------------------------------------------

void DSMaterial::SetMaterialIndexes() { 
  DSLog(development) <<  "----------------------------------------"  <<   endlog ;
  DSLog(development) <<  "              Material indexes          "  <<   endlog ;
  DSLog(development) <<  "----------------------------------------"  <<   endlog ;
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  for(int i=0;i< (int) G4Material::GetNumberOfMaterials(); ++i ) {
    G4Material* aMaterial = (*theMaterialTable)[i];
    DSLog(routine) <<  aMaterial->GetName()  << ": " <<  (int)aMaterial->GetIndex() << endlog ;
    
    if(aMaterial->GetName() == "BoronScintillator") DSStorage::Get()->SetBoronScintillatorIndex((int)aMaterial->GetIndex()) ;
    if(aMaterial->GetName() == "LiquidArgon") DSStorage::Get()->SetLiquidArgonIndex((int)aMaterial->GetIndex()) ;
    if(aMaterial->GetName() == "GaseousArgon") DSStorage::Get()->SetGaseousArgonIndex((int)aMaterial->GetIndex()) ;   
  }
  DSLog(development) <<  "----------------------------------------" <<  endlog ;
  
  DSEventHandler::Get()->SetLArIndex(GetLiquidArgon()->GetIndex()); 
  if     (DSStorage::Get()->GetScintillator()==0)  
    DSEventHandler::Get()->SetScintillatorIndex(GetBoronScintillator()->GetIndex()); 
  else if(DSStorage::Get()->GetScintillator()==1)  
    DSEventHandler::Get()->SetScintillatorIndex(GetGdScintillator()->GetIndex()); 
}

//---------------------------------------------------------------------------------------

G4double DSMaterial::GetLArRefIndex(G4double lambda) {

  G4double epsilon;
  G4double LArRho = 1.390;
  G4double ArRho  = 0.001784;
  if (lambda <= 107.05) return 1.0e4; // lambda MUST be > 107.05 nm
  epsilon = lambda / 1000.0; // switch to micrometers
  epsilon = 1.0 / (epsilon * epsilon); // 1 / (lambda)^2
  epsilon = 1.2055e-2 * ( 0.2075 / (91.012 - epsilon) +
                          0.0415 / (87.892 - epsilon) +
                          4.3330 / (214.02 - epsilon) );
  epsilon *= (8./12.); // Bideau-Sellmeier -> Clausius-Mossotti
  epsilon *= (LArRho / ArRho); // density correction (Ar gas -> LAr liquid)
  if ((epsilon < 0.0) || (epsilon > 0.999999)) return 4.0e6;
  epsilon = (1.0 + 2.0 * epsilon) / (1.0 - epsilon); // solve Clausius-Mossotti

  return sqrt(epsilon);
}
//---------------------------------------------------------------------------------------

G4double DSMaterial::GetGArRefIndex(G4double lambda) {

  G4double epsilon;

  if (lambda <= 107.05) return 1.0e4; // lambda MUST be > 107.05 nm
  epsilon = lambda / 1000.0; // switch to micrometers
  epsilon = 1.0 / (epsilon * epsilon); // 1 / (lambda)^2
  epsilon = 1.2055e-2 * ( 0.2075 / (91.012 - epsilon) +
                          0.0415 / (87.892 - epsilon) +
                          4.3330 / (214.02 - epsilon) );
  epsilon *= (8./12.); // Bideau-Sellmeier -> Clausius-Mossotti
  if ((epsilon < 0.0) || (epsilon > 0.999999)) return 4.0e6;
  epsilon = (1.0 + 2.0 * epsilon) / (1.0 - epsilon); // solve Clausius-Mossotti

  return sqrt(epsilon);
}

//---------------------------------------------------------------------------------------

G4double DSMaterial::GetLArRayLength(G4double lambda) {
  G4double LArT   = 87.0;                 // the actual temperature of LAr in detector
  G4double LArKT  = 2.18e-10 * cm*cm/(1.e-5*newton); // LAr isothermal compressibility
  G4double k      = 1.380658e-23;            // the Boltzmann constant

  G4double h = GetLArEpsilon(lambda);
  if (h < 1.00000001) h = 1.00000001; // just a precaution
  h = (h - 1.0) * (h + 2.0); // the "dielectric constant" dependance
  h *= h; // take the square
  h *= LArKT * LArT *  k  ; // compressibility * temperature * Boltzmann constant
  h /= lambda * lambda * lambda * lambda * 1.0e-36; // (lambda)^4
  h *= 9.18704494231105429; // (2 * Pi / 3)^3
//   if ( h < (1.0 / (10.0 * km)) ) h = 1.0 / (10.0 * km); // just a precaution
//   if ( h > (1.0 / (0.1 * nanometer)) ) h = 1.0 / (0.1 * nanometer); // just a precaution
  //  return ( 100.0 / h );
  return ( 100./h )*um;
}

//---------------------------------------------------------------------------------------

G4double DSMaterial::GetGArRayLength(G4double lambda) {
  G4double ArT   = 87.0;                 // the actual temperature of LAr in detector
  G4double ArKT  = 2.18e-5*cm*cm/(1.e-5*newton); // LAr isothermal compressibility
  G4double k     = 1.380658e-23;            // the Boltzmann constant
  
  G4double h = GetGArEpsilon(lambda);
  if (h < 1.00000001) h = 1.00000001; // just a precaution
  h = (h - 1.0) * (h + 2.0); // the "dielectric constant" dependance
  h *= h; // take the square
  h *= ArKT * ArT * k; // compressibility * temperature * Boltzmann constant
  h /= lambda * lambda * lambda * lambda * 1.0e-36; // (lambda)^4
  h *= 9.18704494231105429; // (2 * Pi / 3)^3
//   if ( h < (1.0 / (10.0 * km)) ) h = 1.0 / (10.0 * km); // just a precaution
//   if ( h > (1.0 / (0.1 * nanometer)) ) h = 1.0 / (0.1 * nanometer); // just a precaution
  return ( 100.0 / h )*um;
}


//---------------------------------------------------------------------------------------
G4double DSMaterial::GetLArEpsilon(G4double lambda) {
  G4double LArRho = 1.390; 
  G4double ArRho  = 0.001784;
  G4double epsilon;

  if (lambda <= 107.05) return 1.0e4; // lambda MUST be > 107.05 nm
  epsilon = lambda / 1000.0; // switch to micrometers
  epsilon = 1.0 / (epsilon * epsilon); // 1 / (lambda)^2
  epsilon = 1.2055e-2 * ( 0.2075 / (91.012 - epsilon) +
                          0.0415 / (87.892 - epsilon) +
                          4.3330 / (214.02 - epsilon) );
  epsilon *= (8./12.); // Bideau-Sellmeier -> Clausius-Mossotti
  epsilon *= (LArRho / ArRho); // density correction (Ar gas -> LAr liquid)
  if ((epsilon < 0.0) || (epsilon > 0.999999)) return 4.0e6;
  epsilon = (1.0 + 2.0 * epsilon) / (1.0 - epsilon); // solve Clausius-Mossotti

  return epsilon;
}
//---------------------------------------------------------------------------------------
G4double DSMaterial::GetGArEpsilon(G4double lambda) {

  G4double epsilon;

  if (lambda <= 107.05) return 1.0e4; // lambda MUST be > 107.05 nm
  epsilon = lambda / 1000.0; // switch to micrometers
  epsilon = 1.0 / (epsilon * epsilon); // 1 / (lambda)^2
  epsilon = 1.2055e-2 * ( 0.2075 / (91.012 - epsilon) +
                          0.0415 / (87.892 - epsilon) +
                          4.3330 / (214.02 - epsilon) );
  epsilon *= (8./12.); // Bideau-Sellmeier -> Clausius-Mossotti
  if ((epsilon < 0.0) || (epsilon > 0.999999)) return 4.0e6;
  epsilon = (1.0 + 2.0 * epsilon) / (1.0 - epsilon); // solve Clausius-Mossotti

  return epsilon;
}


/*
 * $Log: DSMaterial.cc,v $
 * Revision 1.10  2015/04/28 10:14:29  pagnes
 * Implementation of the DS50 optics completed
 *
 * Revision 1.9  2015/04/17 14:51:25  dfranco
 * added water loaded with gadolinium
 *
 * Revision 1.8  2015/03/09 15:20:37  pagnes
 * DS 5tons geometry added (conf 9)
 *
 * Revision 1.7  2015/01/17 11:31:53  pagnes
 * PAr model added form optical tuning
 *
 * Revision 1.6  2015/01/07 16:45:54  pagnes
 * changed veto optical properties format from arrays to vectors
 *
 * Revision 1.5  2014/11/20 15:32:05  dfranco
 * added a command to remove scintillation process from liquid argon between TPC and cryostat
 *
 * Revision 1.4  2014/11/05 15:47:14  pagnes
 * temporary optics tuning
 *
 * Revision 1.3  2014/05/08 10:59:19  pagnes
 * Scintillator Index added in binary header
 *
 * Revision 1.2  2014/05/07 14:27:31  dfranco
 * fixed some bugs and added GdScintillator
 *
 * Revision 1.1  2014/05/07 12:21:04  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.54  2014/04/18 16:21:42  swesterd
 * fixed an overlap in the G2 detector
 *
 * Revision 1.53  2014/04/11 12:33:29  pagnes
 * command to set TMB/PC ratio inside veto added
 *
 * Revision 1.52  2014/04/11 08:38:37  meregaglia
 * corrected B11 mass and TMB composition
 *
 * Revision 1.51  2014/04/10 16:42:18  pagnes
 * TMB added to BoronScintillator material
 *
 * Revision 1.50  2014/03/31 17:38:17  dfranco
 * latest optics config + DSLigh2 scintillation scaling + DSMaterial fixes
 *
 * Revision 1.49  2014/03/19 16:37:26  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.48  2014/03/11 16:49:56  dfranco
 * read optical parameters from external file (DSOptics.dat)
 *
 * Revision 1.47  2014/03/11 09:51:23  dfranco
 * preliminary optical optimization
 *
 * Revision 1.46  2013/10/20 16:30:14  swesterd
 * updated the veto scintillator optical properties based on Aldos measurements
 *
 * Revision 1.45  2013/08/27 04:07:01  swesterd
 * some fine tuning of the boron scintillator kB and scint yield, and some modifications to the DSG2 geometry
 *
 * Revision 1.44  2013/08/20 03:25:52  swesterd
 * added G2 TPC geoemtry (not complete) and added monoenergetic energy distribution to generator
 *
 * Revision 1.43  2013/08/05 12:25:13  swesterd
 * Added G4OpBoundaryProcess.hh so that light can pass through some optical boundaries
 *
 * Revision 1.42  2013/08/05 10:56:14  perassos
 * Grid added to DS50; GridSteel Rindex defined; GridSteel set as Grid Material
 *
 * Revision 1.41  2013/08/05 03:13:52  swesterd
 * some fine tuning of bscint and veto parameters
 *
 * Revision 1.40  2013/06/22 10:27:07  dfranco
 * Changed grid material in DS10. Now the LY at 10 keV is about 9 pe/keV
 *
 * Revision 1.39  2013/06/22 07:21:21  dfranco
 * Fixed a bug in the photoelectron absorption in DS50. Back to the previous QE method
 *
 * Revision 1.38  2013/06/21 13:09:18  dfranco
 * Added DS10 optical surfaces
 *
 * Revision 1.37  2013/06/19 18:35:28  swesterd
 * added DSScintCelll and made tpc PMTs' QE and reflections work like veto PMTs
 *
 * Revision 1.36  2013/06/06 23:19:35  swesterd
 * added quenching to BScint, currently commented out in DSPhysicsList. Have not compared with physical measurements...but might be giving reasonable results? Hard to tell....
 *
 * Revision 1.35  2013/06/05 23:03:32  swesterd
 * moved optical boundary MPTs to DSMaterial and gave the trunks optical boundary properties consistent with untreated stainless steel
 *
 * Revision 1.34  2013/06/05 16:28:30  dfranco
 * Improved optics in the test detector
 *
 * Revision 1.33  2013/06/04 16:56:41  dfranco
 * Added fake BlackHole material to abosrb photons
 *
 * Revision 1.32  2013/06/04 14:39:35  dfranco
 * fixed the teflon rindex
 *
 * Revision 1.31  2013/06/04 12:17:59  dfranco
 * Added a cut to the TPB emission spectrum after 600 nm
 *
 * Revision 1.30  2013/06/04 01:02:29  swesterd
 * other than the optical boundary of the trunks, the veto optics appear to be complete and up and running...modulo whatever I may have missed...
 *
 * Revision 1.29  2013/05/31 17:03:37  dfranco
 * Removed SetSpline(true) (properties inteporlator) since sometimes it fails, as reported in http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/467.html?inline=-1. Thinking about remove it completely from each material property. TPB refrective index temporary set to 1.2, instead of 1.47, in order to avoid total internal relfections, which strongly slow down the simulation
 *
 * Revision 1.28  2013/05/31 16:19:46  dfranco
 * removed ABSLENGTH from TPB properties since it is in contrast with WLSABSLENGTH
 *
 * Revision 1.27  2013/05/31 15:40:06  dfranco
 * fixed some optical properties of the TPC
 *
 * Revision 1.26  2013/05/31 13:28:26  dfranco
 * Added WLSEFFICIENCY to TPB
 *
 * Revision 1.25  2013/05/31 13:02:15  dfranco
 * Added a detector tester, with simpplified geometry (configuration number = 4) to test optical properties of the materials
 *
 * Revision 1.24  2013/05/30 12:34:53  dfranco
 * Fixing the optical properties of the TPC. Not yet concluded
 *
 * Revision 1.23  2013/05/29 16:40:01  dfranco
 * changed optical properties of the materials in order to have consistent array lengths within each material
 *
 * Revision 1.22  2013/05/27 23:59:02  swesterd
 * added a (currently commented out) Lumirror sheath to the cryostat and introduced DSOpBoundaryProcess to try to figure out why the boundaries are being screwy, with some edits so that it can handle constant and vector properties with freaking out
 *
 * Revision 1.21  2013/05/27 11:28:20  dfranco
 * Segmentation fault bug fixed, for the TPC case only. Useless files removed. Fixed the RINDEX, which was not defined in the whole range, for the scintillator.
 *
 * Revision 1.20  2013/05/25 07:58:23  swesterd
 * Got the veto PMT optical boundaries all working along with photocathode optical properties, added PMT quantum efficiency to DSTrackingAction, and added a function to DSTrackingAction that locates and quadratically interpolates points in data, for getting useful QEs
 *
 * Revision 1.19  2013/05/14 14:35:31  dfranco
 * Changed log message levels from trace to warnign or debugging
 *
 * Revision 1.18  2013/05/07 16:15:45  swesterd
 * Optical processes now seem to be working in the boron-loaded scintillator
 *
 * Revision 1.17  2013/05/06 14:59:52  perassos
 * Updates on the TPC surface properties and geometry
 *
 * Revision 1.16  2013/05/02 16:10:01  dfranco
 * TPB WLS fixed
 *
 * Revision 1.15  2013/05/01 08:20:23  swesterd
 * added boron-loaded scintillator optical properties
 *
 * Revision 1.14  2013/04/30 14:33:15  dfranco
 * fixed bug with rayleigh scattering in GAr
 *
 * Revision 1.13  2013/04/30 14:31:07  dfranco
 * fixed bug with Rayleigh scattering in LAr
 *
 * Revision 1.12  2013/04/19 16:24:14  dfranco
 * Added Rayleigh scattering to liquid and gaseous argon
 *
 * Revision 1.11  2013/04/19 15:49:36  dfranco
 * WLS absorption length updated
 *
 * Revision 1.10  2013/04/19 13:36:49  dfranco
 * Added absorption length
 *
 * Revision 1.9  2013/04/19 10:46:22  dfranco
 * added TPB properties
 *
 * Revision 1.8  2013/04/19 10:22:12  meregaglia
 * DSLight first major cleaning
 *
 * Revision 1.7  2013/04/19 09:26:38  dfranco
 * added TPB emission spectrum
 *
 * Revision 1.6  2013/04/11 09:10:17  meregaglia
 * work in progress to clean S1
 *
 * Revision 1.5  2013/04/10 21:06:05  meregaglia
 * added S1S2 from NEST and possibility to switch on or off this physics. Work in progress on the physics itsel
 *
 * Revision 1.4  2013/03/22 17:48:02  dfranco
 * added refrective indexes for gasous and liquid argon materials
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
