//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The NEST program is intended for use with the Geant4 software,   *
// * which is copyright of the Copyright Holders of the Geant4        *
// * Collaboration. This additional software is copyright of the NEST *
// * development team. As such, it is subject to the terms and        *
// * conditions of both the Geant4 License, included with your copy   *
// * of Geant4 and available at http://cern.ch/geant4/license, as     *
// * well as the NEST License included with the download of NEST and  *
// * available at http://nest.physics.ucdavis.edu/                    *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutions, nor the agencies providing financial support for   *
// * this work make any representation or warranty, express or        *
// * implied, regarding this software system, or assume any liability *
// * for its use. Please read the pdf license or view it online       *
// * before download for the full disclaimer and lack of liability.   *
// *                                                                  *
// * This code implementation is based on work by Peter Gumplinger    *
// * and his fellow collaborators on Geant4 and is distributed with   *
// * the express written consent of the Geant4 collaboration. By      *
// * using, copying, modifying, or sharing the software (or any work  *
// * based on the software) you agree to acknowledge use of both NEST *
// * and Geant4 in resulting scientific publications, and you         *
// * indicate your acceptance of all the terms and conditions of the  *
// * licenses, which must always be included with this code.          *
// ********************************************************************
//
//
// Based on NEST and optimized for LAr and DS (A.Meregaglia)
//
////////////////////////////////////////////////////////////////////////
// S1 S2  Scintillation Light Class Implementation
////////////////////////////////////////////////////////////////////////

#include "G4ParticleTypes.hh" //lets you refer to G4OpticalPhoton, etc.
#include "G4EmProcessSubType.hh" //lets you call this process Scintillation
#include "G4Version.hh" //tells you what Geant4 version you are running
#include "DSMaterial.hh"
#include "DSStorage.hh"
#include "DSLight.hh"
#include "DSLogger.hh"
#include "DSEventHandler.hh"
#include <iostream>

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#define MIN_ENE -1*eV //lets you turn NEST off BELOW a certain energy
#define MAX_ENE 1.*TeV //lets you turn NEST off ABOVE a certain energy
#define HIENLIM 5*MeV //energy at which Doke model used exclusively

#define R_TOL 0.2*mm //tolerance (for edge events)
#define R_MAX 1*km //for corraling diffusing electrons


using namespace std ;
G4double thrLAr; //offset in linear light yield formula for S2 ph/e-
G4double E_eVLAr; //energy of single photon in gas, eV
G4double tau1LAr, tau3LAr, MolarMassLAr, ConvertEffLAr;


G4bool diffusion = true;
G4bool SinglePhase=false, ThomasImelTail=true, OutElectrons=true;

G4double biExc = 0.77; //for alpha particles (bi-excitonic collisions)

#define Density_LAr 1.393
#define TaueLAr 15.8*ms//electron lifetime due to impurities //http://arxiv.org/pdf/0910.5087.pdf



///////////////////////////////////////////////////
DSLight::DSLight(const G4String& processName,G4ProcessType type)
  : G4VRestDiscreteProcess(processName, type)
{
  thrLAr = 0.190; E_eVLAr = 9.7; MolarMassLAr = 39.948;
  tau1LAr = 6*ns; tau3LAr = 1600*ns; ConvertEffLAr=0.7857;
  
  //SetProcessSubType(fScintillation);
  SetProcessType(fUserDefined);
  
  fTrackSecondariesFirst = false;
  //particles die first, then scintillation is generated
  
  DSLog(debugging) << GetProcessName() << " is created "<< endlog ;    

}
///////////////////////////////////////////////////
DSLight::~DSLight(){} //destructor needed to avoid linker error

///////////////////////////////////////////////////
G4VParticleChange* DSLight::AtRestDoIt(const G4Track& aTrack, const G4Step& aStep)
// This routine simply calls the equivalent PostStepDoIt since all the
// necessary information resides in aStep.GetTotalEnergyDeposit()
{
  return DSLight::PostStepDoIt(aTrack, aStep);
}

///////////////////////////////////////////////////
G4VParticleChange* DSLight::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
// this is the most important function, where all light & charge yields happen!
{
 
 
  aParticleChange.Initialize(aTrack);
	
  if ( !YieldFactor ) //set YF=0 when you want DSLight off in your sim
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  

  // if parent is the primary and this is the first step, initialize variables
  if( aTrack.GetParentID() == 0 && aTrack.GetCurrentStepNumber() == 1 ) {
    fExcitedNucleus = false; //an initialization or reset
    fVeryHighEnergy = false; //initializes or (later) resets this
    fAlpha = false; //ditto
    fMultipleScattering = false;
  }

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4ParticleDefinition *pDef = aParticle->GetDefinition();
  G4String particleName = pDef->GetParticleName();
  const G4Material* aMaterial = aStep.GetPreStepPoint()->GetMaterial();
  const G4Material* bMaterial = aStep.GetPostStepPoint()->GetMaterial();

  // skip neutrons   	
  if((particleName == "neutron" || particleName == "antineutron") && aStep.GetTotalEnergyDeposit() <= 0)
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);


	
  // code for determining whether the present/next material is noble
  // element, or, in other words, for checking if either is a valid NEST
  // scintillating material, and save Z for later L calculation, or
  // return if no valid scintillators are found on this step, which is
  // protection against G4Exception or seg. fault/violation
  G4Element *ElementA = NULL, *ElementB = NULL;
  if (aMaterial) {
    const G4ElementVector* theElementVector1 = aMaterial->GetElementVector();
    ElementA = (*theElementVector1)[0];
  }
  if (bMaterial) {
    const G4ElementVector* theElementVector2 = bMaterial->GetElementVector();
    ElementB = (*theElementVector2)[0];
  }
  
  G4int z1,z2,TOTALNUM_INT_SITES=1; 
  G4bool NobleNow=false,
    NobleLater=false;
  if (ElementA) 
    z1 = (G4int)(ElementA->GetZ()); 
  else 
    z1 = -1;

  if (ElementB) 
    z2 = (G4int)(ElementB->GetZ()); 
  else 
    z2 = -1;
  if (z1==18) {  //only Argon
    NobleNow = true;
    TOTALNUM_INT_SITES = (G4int)aMaterial->GetMaterialPropertiesTable()->GetConstProperty("TOTALNUM_INT_SITES"); //get current number
    
    if ( TOTALNUM_INT_SITES < 0 ) {
      InitMatPropValues(aMaterial->GetMaterialPropertiesTable());
      TOTALNUM_INT_SITES = 0; //no sites yet
    } //material properties initialized
  } //end of atomic number check
  if (z2==18){ //only Argon
    NobleLater = true;
    TOTALNUM_INT_SITES = (G4int)bMaterial->GetMaterialPropertiesTable()->GetConstProperty("TOTALNUM_INT_SITES");
    if ( TOTALNUM_INT_SITES < 0 ) {
      InitMatPropValues(bMaterial->GetMaterialPropertiesTable());
      TOTALNUM_INT_SITES = 0; //no sites yet
    } //material properties initialized
  } //end of atomic number check
	
  if ( !NobleNow && !NobleLater )
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);


	
  // retrieval of the particle's position, time, attributes at both the 
  // beginning and the end of the current step along its track
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
  G4ThreeVector x1 = pPostStepPoint->GetPosition();
  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  G4double evtStrt = pPreStepPoint->GetGlobalTime();
  G4double      t0 = pPreStepPoint->GetLocalTime();
  G4double      t1 = pPostStepPoint->GetLocalTime();


	
  // now check if we're entering a scintillating material (inside) or
  // leaving one (outside), in order to determine (later on in the code,
  // based on the booleans inside & outside) whether to add/subtract
  // energy that can potentially be deposited from the system
  G4bool outside = false, inside = false, InsAndOuts = false;
  G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();
  if ( NobleNow && !NobleLater ) outside = true;
  else if ( !NobleNow && NobleLater ) {
    aMaterial = bMaterial;
    inside = true; 
    z1 = z2;
    aMaterialPropertiesTable = bMaterial->GetMaterialPropertiesTable();
  }
  else if ( NobleNow && NobleLater && aMaterial->GetDensity() != bMaterial->GetDensity() )
    InsAndOuts = true;
  
    
	
  // retrieve scintillation-related material properties
  G4double Density = aMaterial->GetDensity()/(g/cm3);
  G4double nDensity = Density*AVO; //molar mass factor applied below
  G4int Phase = aMaterial->GetState(); //solid, liquid, or gas?
  G4double ElectricField, FieldSign; //for field quenching of S1
  ElectricField = aMaterialPropertiesTable-> GetConstProperty("ELECTRICFIELD"); 
 
  
  if ( ElectricField >= 0 ) 
    FieldSign = 1; 
  else 
    FieldSign = -1;

  ElectricField = fabs((1e3*ElectricField)/(kilovolt/cm));//ANS test units
  G4double Temperature = aMaterial->GetTemperature();
  G4double ScintillationYield, ResolutionScale, R0 = 1.0*um, DokeBirks[3], ThomasImel = 0.00, delta = 1*mm;
  DokeBirks[0] = 0.00; DokeBirks[2] = 1.00;
  G4double PhotMean = 7*eV, PhotWidth = 1.0*eV; //photon properties
  G4double SingTripRatioR, SingTripRatioX, tau1, tau3, tauR = 0*ns;

  //Argon  //ANS parameters to be tuned
  //ScintillationYield = 1.34077 * 1 / (19.5*eV); // scale factor for increasing the yield
  ScintillationYield = 0.642762 * 1 / (19.5*eV); // scale factor for increasing the yield
  ExcitationRatio = 0.21; //Aprile et. al book
  //ResolutionScale = 0.107; //Doke 1976
  ResolutionScale = 0.20; //Doke 1984
  R0 = 1.568*um; //Mozumder 1995
  if(ElectricField) {
    ThomasImel = 0.156977*pow(ElectricField,-0.1);
    DokeBirks[0] = 0.07*pow((ElectricField/1.0e3),-0.85);
//    DokeBirks[2] = 0.00;
    DokeBirks[2] = 0.80;
  }
  else {
    ThomasImel = 0.099;
    DokeBirks[0] = 0.0003;
    DokeBirks[2] = 0.75;
  } 
  nDensity /= 39.948; //molar mass in grams per mole
  PhotMean = 9.69*eV; 
  PhotWidth = 0.22*eV;
  tau1 = 7*ns;//Hitachi, A. et al., Effect of ionization density on the time dependence of luminescence from liquid argon and xenon, Phys. Rev. B 27 (1983), no. 9 5279.
  tau3 = 1600*ns;//same as above 
  tauR = 0.8*ns; //Kubota 1979
  biExc = 0.6;


  // log present and running tally of energy deposition in this section
  G4double anExcitationEnergy = ((const G4Ions*)(pDef))->GetExcitationEnergy(); //grab nuclear energy level
  G4double TotalEnergyDeposit = aMaterialPropertiesTable->GetConstProperty( "ENERGY_DEPOSIT_TOT" ); //total energy deposited so far
  G4bool convert = false, annihil = false;
  
  //set up special cases for pair production and positron annihilation
  if(pPreStepPoint->GetKineticEnergy()>=(2*electron_mass_c2) && 
     !pPostStepPoint->GetKineticEnergy() && 
     !aStep.GetTotalEnergyDeposit() && aParticle->GetPDGcode()==22) {
    convert = true; TotalEnergyDeposit = electron_mass_c2;
  }
  if(pPreStepPoint->GetKineticEnergy() && 
     !pPostStepPoint->GetKineticEnergy() && 
     aParticle->GetPDGcode()==-11) {
    annihil = true; TotalEnergyDeposit += aStep.GetTotalEnergyDeposit();
  }

  G4bool either = false;
  if(inside || outside || convert || annihil || InsAndOuts) 
    either=true;
    
  //conditions for returning when energy deposits too low
  if( anExcitationEnergy<100*eV && aStep.GetTotalEnergyDeposit()<1*eV && !either && !fExcitedNucleus )
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

  //add current deposit to total energy budget
  if ( !annihil ) TotalEnergyDeposit += aStep.GetTotalEnergyDeposit();
  if ( !convert ) aMaterialPropertiesTable->AddConstProperty( "ENERGY_DEPOSIT_TOT", TotalEnergyDeposit );

  //save current deposit for determining number of quanta produced now
  TotalEnergyDeposit = aStep.GetTotalEnergyDeposit();
	

  // check what the current "goal" E is for dumping scintillation,
  // often the initial kinetic energy of the parent particle, and deal
  // with all other energy-related matters in this block of code
  G4double InitialKinetEnergy = aMaterialPropertiesTable->GetConstProperty( "ENERGY_DEPOSIT_GOL" );
  //if zero, add up initial potential and kinetic energies now
  if ( InitialKinetEnergy == 0 ) {
    G4double tE = pPreStepPoint->GetKineticEnergy()+anExcitationEnergy;
    aMaterialPropertiesTable->AddConstProperty ( "ENERGY_DEPOSIT_GOL", tE );
    
    //excited nucleus is special case where accuracy reduced for total
    //energy deposition because of G4 inaccuracies and scintillation is
    //forced-dumped when that nucleus is fully de-excited
    if ( anExcitationEnergy ) 
      fExcitedNucleus = true;
  }

  //if a particle is leaving, remove its kinetic energy from the goal
  //energy, as this will never get deposited (if depositable)
  if(outside){ 
    aMaterialPropertiesTable->AddConstProperty("ENERGY_DEPOSIT_GOL",InitialKinetEnergy-pPostStepPoint->GetKineticEnergy());
    if(aMaterialPropertiesTable->GetConstProperty("ENERGY_DEPOSIT_GOL")<0)
      aMaterialPropertiesTable->AddConstProperty("ENERGY_DEPOSIT_GOL",0);
  }

  //if a particle is coming back into your scintillator, then add its
  //energy to the goal energy
  if(inside) { 
    aMaterialPropertiesTable->AddConstProperty("ENERGY_DEPOSIT_GOL",InitialKinetEnergy+pPreStepPoint->GetKineticEnergy());
    if ( TotalEnergyDeposit > 0 && InitialKinetEnergy == 0 ) {
      aMaterialPropertiesTable->AddConstProperty("ENERGY_DEPOSIT_GOL",0);
      TotalEnergyDeposit = .000000;
    }
  }

  if ( InsAndOuts ) {
    aMaterialPropertiesTable->AddConstProperty("ENERGY_DEPOSIT_GOL",(-0.1*keV)+InitialKinetEnergy-pPostStepPoint->GetKineticEnergy());
    InitialKinetEnergy = bMaterial->GetMaterialPropertiesTable()->GetConstProperty("ENERGY_DEPOSIT_GOL");
    bMaterial->GetMaterialPropertiesTable()->AddConstProperty("ENERGY_DEPOSIT_GOL",(-0.1*keV)+InitialKinetEnergy+pPreStepPoint->GetKineticEnergy());
    if(aMaterialPropertiesTable->GetConstProperty("ENERGY_DEPOSIT_GOL")<0)
      aMaterialPropertiesTable->AddConstProperty("ENERGY_DEPOSIT_GOL",0);
    if ( bMaterial->GetMaterialPropertiesTable()->GetConstProperty("ENERGY_DEPOSIT_GOL") < 0 )
      bMaterial->GetMaterialPropertiesTable()->AddConstProperty ( "ENERGY_DEPOSIT_GOL", 0 );
  }

  InitialKinetEnergy = aMaterialPropertiesTable->GetConstProperty("ENERGY_DEPOSIT_GOL"); //grab current goal E
  if ( annihil ) //if an annihilation occurred, add energy of two gammas
    InitialKinetEnergy += 2*electron_mass_c2;
  //if pair production occurs, then subtract energy to cancel with the
  //energy that will be added in the line above when the e+ dies
  if ( convert )
    InitialKinetEnergy -= 2*electron_mass_c2;
  //update the relevant material property (goal energy)
  aMaterialPropertiesTable->AddConstProperty("ENERGY_DEPOSIT_GOL",InitialKinetEnergy);
  
  if (anExcitationEnergy < 1e-100 && aStep.GetTotalEnergyDeposit()==0 &&
      aMaterialPropertiesTable->GetConstProperty("ENERGY_DEPOSIT_GOL")==0 &&
      aMaterialPropertiesTable->GetConstProperty("ENERGY_DEPOSIT_TOT")==0)
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);


  G4String procName;
  if ( aTrack.GetCreatorProcess() )
    procName = aTrack.GetCreatorProcess()->GetProcessName();
  else
    procName = "NULL";
  if ( procName == "eBrem" && outside && !OutElectrons ) //ANS check OutElectrons set true by hand???
    fMultipleScattering = true;

  // next 2 codeblocks deal with position-related things
  //ANS code not really clear...
  if ( fAlpha ) 
    delta = 1000.*km;

  G4int i, k, counter = 0; G4double pos[3];
  if ( outside ) { //leaving
    if ( aParticle->GetPDGcode() == 11 && !OutElectrons ) //ANS check OutElectrons set true by hand???
      fMultipleScattering = true;
    x1 = x0; //prevents generation of quanta outside active volume
  } //no scint. for e-'s that leave

  char xCoord[80]; char yCoord[80]; char zCoord[80];
  G4bool exists = false; //for querying whether set-up of new site needed
  for(i=0;i<TOTALNUM_INT_SITES;i++) { //loop over all saved interaction sites
    counter = i; //save site# for later use in storing properties
    sprintf(xCoord,"POS_X_%d",i); 
    sprintf(yCoord,"POS_Y_%d",i);
    sprintf(zCoord,"POS_Z_%d",i);
    pos[0] = aMaterialPropertiesTable->GetConstProperty(xCoord);
    pos[1] = aMaterialPropertiesTable->GetConstProperty(yCoord);
    pos[2] = aMaterialPropertiesTable->GetConstProperty(zCoord);
    if ( sqrt(pow(x1[0]-pos[0],2.)+pow(x1[1]-pos[1],2.)+pow(x1[2]-pos[2],2.)) < delta ) {
      exists = true; 
      break; //we find interaction is close to an old one
    }
  }
  
  if(!exists && TotalEnergyDeposit) { //current interaction too far away
    counter = TOTALNUM_INT_SITES;
    sprintf(xCoord,"POS_X_%i",TOTALNUM_INT_SITES); 
    sprintf(yCoord,"POS_Y_%i",TOTALNUM_INT_SITES); 
    sprintf(zCoord,"POS_Z_%i",TOTALNUM_INT_SITES);
    //save 3-space coordinates of the new interaction site
    aMaterialPropertiesTable->AddConstProperty( xCoord, x1[0] );
    aMaterialPropertiesTable->AddConstProperty( yCoord, x1[1] );
    aMaterialPropertiesTable->AddConstProperty( zCoord, x1[2] );
    TOTALNUM_INT_SITES++; //increment number of sites
    aMaterialPropertiesTable->AddConstProperty( "TOTALNUM_INT_SITES", TOTALNUM_INT_SITES );//save
  }
  
///////// up to now he saved all the steppoints of interaction.  
///////// He will check if the points are too close one to each other
////////  From this points, he will generate "thermalelectrons" for 
////////  the S2 process (or recombination)







  // this is where nuclear recoil "L" factor is handled: total yield is
  // reduced for nuclear recoil as per Lindhard theory
	
  //we assume you have a mono-elemental scintillator only
  //now, grab A's and Z's of current particle and of material (avg)
  G4double a1 = ElementA->GetA();
  z2 = pDef->GetAtomicNumber(); 
  G4double a2 = (G4double)(pDef->GetAtomicMass());
  if ( particleName == "alpha" || (z2 == 2 && a2 == 4) )
    fAlpha = true; //used later to get S1 pulse shape correct for alpha
  if ( fAlpha || abs(aParticle->GetPDGcode()) == 2112 )
    a2 = a1; //get average A for element at hand
  G4double epsilon = 11.5*(TotalEnergyDeposit/keV)*pow(z1,(-7./3.));
  G4double gamma = 3.*pow(epsilon,0.15)+0.7*pow(epsilon,0.6)+epsilon;
  G4double kappa = 0.133*pow(z1,(2./3.))*pow(a2,(-1./2.))*(2./3.);
  
  //check if we are dealing with nuclear recoil (Z same as material)
  if ( (z1 == z2 && pDef->GetParticleType() == "nucleus") || particleName == "neutron" || particleName == "antineutron" ) {
    YieldFactor=(kappa*gamma)/(1+kappa*gamma); //Lindhard factor
    if ( z1 == 18 && Phase == kStateLiquid )
      YieldFactor=0.23*(1+exp(-5*epsilon)); //liquid argon L_eff
    //just a few safety checks, like for recombProb below
    if ( YieldFactor > 1 ) YieldFactor = 1;
    if ( YieldFactor < 0 ) YieldFactor = 0;
    if ( ElectricField == 0 && Phase == kStateLiquid ) {
      if ( z1 == 18 ) ThomasImel = 0.25;//LAr
    } //special TIB parameters for nuclear recoil only, in LAr
    ExcitationRatio = 0.69337 + 0.3065*exp(-0.008806*pow(ElectricField,0.76313));
  }
  else YieldFactor = 1.000; //default
  

	
  // determine ultimate #quanta from current E-deposition (ph+e-)
  G4double MeanNumberOfQuanta = ScintillationYield*TotalEnergyDeposit; //total mean number of exc/ions
  
  //the total number of either quanta produced is equal to product of the
  //work function, the energy deposited, and yield reduction, for NR
  G4double sigma = sqrt(ResolutionScale*MeanNumberOfQuanta); //Fano
  G4int NumQuanta = G4int(floor(G4RandGauss::shoot(MeanNumberOfQuanta,sigma)+0.5));//stochastic variation in NumQuanta
  G4double LeffVar = G4RandGauss::shoot(YieldFactor,0.25*YieldFactor);
  
  if (LeffVar > 1) LeffVar = 1.00000; 
  if (LeffVar < 0) LeffVar = 0;
  if ( YieldFactor < 1 ) NumQuanta = BinomFluct(NumQuanta,LeffVar);
	
  //if E below work function, can't make any quanta, and if NumQuanta
  //less than zero because Gaussian fluctuated low, update to zero
  if(TotalEnergyDeposit < 1/ScintillationYield || NumQuanta < 0)
    NumQuanta = 0;



///////  fino a qui: ok
////////
/////////////


 
  // next section binomially assigns quanta to excitons and ions
  // ExcitationRatio/(1+ExcitationRatio) = probablity to produce excitons
  G4int NumExcitons = BinomFluct(NumQuanta,ExcitationRatio/(1+ExcitationRatio));
  G4int NumIons     = NumQuanta - NumExcitons;

  // this section calculates recombination following the modified Birks'
  // Law of Doke, deposition by deposition, and may be overridden later
  // in code if a low enough energy necessitates switching to the 
  // Thomas-Imel box model for recombination instead (determined by site)
  G4double dE, dx=0, LET=0, recombProb;
  dE = TotalEnergyDeposit/MeV;

  //in other words, if it's a gamma,ion,proton,alpha,pion,et al. do not
  //use the step length provided by Geant4 because it's not relevant,
  //instead calculate an estimated LET and range of the electrons that
  //would have been produced if Geant4 could track them
  if ( particleName != "e-" && particleName != "e+" && z1 != z2 && particleName != "mu-" && particleName != "mu+" ) {
    LET = CalculateElectronLET( 1000*dE);//LAr
    if(LET) dx = dE/(Density*LET); //find the range based on the LET
    if(abs(aParticle->GetPDGcode())==2112) dx=0;
  }
  //normal case of an e-/+ energy deposition recorded by Geant
  else { 
    dx = aStep.GetStepLength()/cm;
    if(dx) 
      LET = (dE/dx)*(1/Density); //lin. energy xfer (prop. to dE/dx)
    if ( LET > 0 && dE > 0 && dx > 0 ) {
      G4double ratio = CalculateElectronLET(dE*1e3)/LET;//LAr
      if ( TOTALNUM_INT_SITES == 1 && ratio < 0.7 && !ThomasImelTail && particleName == "e-" ) {
	dx /= ratio; 
	LET *= ratio; }
    }
  }
  DokeBirks[1] = DokeBirks[0]/(1-DokeBirks[2]); //B=A/(1-C) (see paper)
  //Doke/Birks' Law as spelled out in the NEST paper
  recombProb = (DokeBirks[0]*LET)/(1+DokeBirks[1]*LET)+DokeBirks[2];
  if ( Phase == kStateLiquid ) {
    if ( z1 == 18 ) recombProb *= (Density/Density_LAr);
  }
  
  //check against unphysicality resulting from rounding errors
  if(recombProb<0) recombProb=0;
  if(recombProb>1) recombProb=1;
  
  //use binomial distribution to assign photons, electrons, where photons
  //are excitons plus recombined ionization electrons, while final
  //collected electrons are the "escape" (non-recombined) electrons
  G4int NumPhotons = NumExcitons + BinomFluct(NumIons,recombProb);
  G4int NumElectrons = NumQuanta - NumPhotons;
	
  // next section increments the numbers of excitons, ions, photons, and
  // electrons for the appropriate interaction site; it only appears to
  // be redundant by saving seemingly no longer needed exciton and ion
  // counts, these having been already used to calculate the number of ph
  // and e- above, whereas it does need this later for Thomas-Imel model
  char numExc[80]; char numIon[80]; char numPho[80]; char numEle[80];
  sprintf(numExc,"N_EXC_%i",counter);
  sprintf(numIon,"N_ION_%i",counter);
  NumExcitons += (G4int)aMaterialPropertiesTable->GetConstProperty( numExc );
  NumIons+= (G4int)aMaterialPropertiesTable->GetConstProperty( numIon );
  aMaterialPropertiesTable->AddConstProperty( numExc, NumExcitons );
  aMaterialPropertiesTable->AddConstProperty( numIon, NumIons     );
  sprintf(numPho,"N_PHO_%i",counter); 
  sprintf(numEle,"N_ELE_%i",counter);

  NumPhotons   +=(G4int)aMaterialPropertiesTable->GetConstProperty( numPho );

  NumElectrons +=(G4int)aMaterialPropertiesTable->GetConstProperty( numEle );
  aMaterialPropertiesTable->AddConstProperty( numPho, NumPhotons   );
  aMaterialPropertiesTable->AddConstProperty( numEle, NumElectrons );
	
  // increment and save the total track length, and save interaction
  // times for later, when generating the scintillation quanta
  char trackL[80]; char time00[80]; char time01[80]; char energy[80];
  sprintf(trackL,"TRACK_%i",counter); 
  sprintf(energy,"ENRGY_%i",counter);
  sprintf(time00,"TIME0_%i",counter); 
  sprintf(time01,"TIME1_%i",counter);
  delta = aMaterialPropertiesTable->GetConstProperty( trackL );
  G4double energ = aMaterialPropertiesTable->GetConstProperty( energy );
  delta += dx*cm; 
  energ += dE*MeV;
  aMaterialPropertiesTable->AddConstProperty( trackL, delta );
  aMaterialPropertiesTable->AddConstProperty( energy, energ );
  
  if ( TotalEnergyDeposit > 0 ) {
    G4double deltaTime = aMaterialPropertiesTable->GetConstProperty( time00 );
    //for charged particles, which continuously lose energy, use initial
    //interaction time as the minimum time, otherwise use only the final
    if (aParticle->GetCharge() != 0) {
      if (t0 < deltaTime)
	aMaterialPropertiesTable->AddConstProperty( time00, t0 );
    }
    else {
      if (t1 < deltaTime)
	aMaterialPropertiesTable->AddConstProperty( time00, t1 );
    }
    deltaTime = aMaterialPropertiesTable->GetConstProperty( time01 );
    //find the maximum possible scintillation "birth" time
    if (t1 > deltaTime)
      aMaterialPropertiesTable->AddConstProperty( time01, t1 );
  }
	
  // begin the process of setting up creation of scint./ionization
  TotalEnergyDeposit=aMaterialPropertiesTable->GetConstProperty("ENERGY_DEPOSIT_TOT"); //get the total E deposited
  InitialKinetEnergy=aMaterialPropertiesTable->GetConstProperty("ENERGY_DEPOSIT_GOL"); //E that should have been
  if(InitialKinetEnergy > HIENLIM && abs(aParticle->GetPDGcode()) != 2112)
    fVeryHighEnergy=true;
  G4double safety; //margin of error for TotalE.. - InitialKinetEnergy
  if (fVeryHighEnergy && !fExcitedNucleus) safety = 0.2*keV;
  else safety = 2.*eV;
	
  //force a scintillation dump for NR and for full nuclear de-excitation
  if( !anExcitationEnergy && pDef->GetParticleType() == "nucleus" && aTrack.GetTrackStatus() != fAlive && !fAlpha )
    InitialKinetEnergy = TotalEnergyDeposit;
  if ( particleName == "neutron" || particleName == "antineutron" )
    InitialKinetEnergy = TotalEnergyDeposit;
	
  //force a dump of all saved scintillation under the following
  //conditions: energy goal reached, and current particle dead, or an 
  //error has occurred and total has exceeded goal (shouldn't happen)
  if( fabs(TotalEnergyDeposit-InitialKinetEnergy)<safety || TotalEnergyDeposit>=InitialKinetEnergy ){
    dx = 0; dE = 0;
    //calculate the total number of quanta from all sites and all
    //interactions so that the number of secondaries gets set correctly
    NumPhotons = 0; NumElectrons = 0;
    for(i=0;i<TOTALNUM_INT_SITES;i++) {
      sprintf(numPho,"N_PHO_%d",i); 
      sprintf(numEle,"N_ELE_%d",i);
      NumPhotons  +=(G4int)aMaterialPropertiesTable->GetConstProperty( numPho );
      NumElectrons+=(G4int)aMaterialPropertiesTable->GetConstProperty( numEle );
      sprintf(trackL,"TRACK_%d",i); 
      sprintf(energy,"ENRGY_%d",i);
      //add up track lengths of all sites, for a total LET calc (later)
      dx += aMaterialPropertiesTable->GetConstProperty(trackL);
      dE += aMaterialPropertiesTable->GetConstProperty(energy);
    }
    G4int buffer = 100; 
    if ( fVeryHighEnergy ) buffer = 1;
    aParticleChange.SetNumberOfSecondaries(buffer*(NumPhotons+NumElectrons));
    if (fTrackSecondariesFirst) {
      if (aTrack.GetTrackStatus() == fAlive )
	aParticleChange.ProposeTrackStatus(fSuspend);
    }

    // begin the loop over all sites which generates all the quanta
    for(i=0;i<TOTALNUM_INT_SITES;i++) {
      // get the position X,Y,Z, exciton and ion numbers, total track 
      // length of the site, and interaction times
      
      sprintf(xCoord,"POS_X_%d",i); 
      sprintf(yCoord,"POS_Y_%d",i);
      sprintf(zCoord,"POS_Z_%d",i);
      sprintf(numExc,"N_EXC_%d",i);
      sprintf(numIon,"N_ION_%d",i);
      sprintf(numPho,"N_PHO_%d",i);
      sprintf(numEle,"N_ELE_%d",i);
      NumExcitons = (G4int)aMaterialPropertiesTable->GetConstProperty( numExc );
      NumIons     = (G4int)aMaterialPropertiesTable->GetConstProperty( numIon );
      sprintf(trackL,"TRACK_%d",i);
      sprintf(energy,"ENRGY_%d",i);
      sprintf(time00,"TIME0_%d",i);
      sprintf(time01,"TIME1_%d",i);
      delta = aMaterialPropertiesTable->GetConstProperty( trackL );
      energ = aMaterialPropertiesTable->GetConstProperty( energy );
      t0 = aMaterialPropertiesTable->GetConstProperty( time00 );
      t1 = aMaterialPropertiesTable->GetConstProperty( time01 );
      //if site is small enough, override the Doke/Birks' model with
      //Thomas-Imel, but not if we're dealing with super-high energy 
      //particles, and if it's NR force Thomas-Imel (though NR should be
      //already short enough in track even up to O(100) keV)
      
// This block starts here ---->
      if ( (delta < R0 && !fVeryHighEnergy) || z2 == z1 || fAlpha ) {
	if ( Phase == kStateLiquid ) {
	  if ( z1 == 18 ) ThomasImel *= pow((Density/Density_LAr),0.3);
	}
	 
	//calculate the Thomas-Imel recombination probability, which
	//depends on energy via NumIons, but not on dE/dx, and protect
	//against seg fault by ensuring a positive number of ions
	if (NumIons > 0) {
	  G4double xi;
	  xi = (G4double(NumIons)/4.)*ThomasImel;
	  if ( InitialKinetEnergy == 9.4*keV ) {//ANS to take out? Xe? Kr?
	    G4double NumIonsEff = 1.066e7*pow(t0/ns,-1.303)*(0.17163+162.32/(ElectricField+191.39));
	    if ( NumIonsEff > 1e6 ) NumIonsEff = 1e6;
	    xi = (G4double(NumIonsEff)/4.)*ThomasImel;
	  }
	  recombProb = 1-log(1+xi)/xi;
	  if(recombProb<0) recombProb=0;
	  if(recombProb>1) recombProb=1;
	}
	//just like Doke: simple binomial distribution
	NumPhotons = NumExcitons + BinomFluct(NumIons,recombProb);
	NumElectrons = (NumExcitons + NumIons) - NumPhotons;
	//override Doke NumPhotons and NumElectrons
	aMaterialPropertiesTable->AddConstProperty( numPho, NumPhotons   );
	aMaterialPropertiesTable->AddConstProperty( numEle, NumElectrons );
      }
// ....and ends here <------------------

      // grab NumPhotons/NumElectrons, which come from Birks if
      // the Thomas-Imel block of code above was not executed
      NumPhotons  = (G4int)aMaterialPropertiesTable->GetConstProperty( numPho );
      NumElectrons =(G4int)aMaterialPropertiesTable->GetConstProperty( numEle );
	    
     
     
      /* DF: in this part they evaluate the extra FanoFactor
             which is never used!!! I comment it
      // extra Fano factor caused by recomb. fluct.
      G4double FanoFactor =0; //ionization channel
      if(Phase == kStateLiquid && YieldFactor == 1) {
	FanoFactor = 2575.9*pow((ElectricField+15.154),-0.64064)-1.4707;
	if ( (dE/keV) <= 100 && ElectricField >= 0 ) {
	  G4double keVee = (TotalEnergyDeposit/(100.*keV));
	  if ( keVee <= 0.06 )
	    FanoFactor *= -0.00075+0.4625*keVee+34.375*pow(keVee,2.);
	  else
	    FanoFactor *= 0.069554+1.7322*keVee-.80215*pow(keVee,2.);
	}
      }
      
      if ( Phase == kStateGas && Density>0.5 ) 
	FanoFactor =0.42857-4.7857*Density+7.8571*pow(Density,2.);
      if( FanoFactor <= 0 || fVeryHighEnergy ) 
	FanoFactor = 0;
      */	
	
	
	
	
      NumQuanta = NumPhotons + NumElectrons; 
      // Next line commented by DF
      //NumPhotons = NumQuanta - NumElectrons;
      if ( NumElectrons <= 0 ) 
	NumElectrons = 0;
      if (   NumPhotons <= 0 )   
	NumPhotons = 0;
      else { //other effects
	if ( fAlpha ) //bi-excitonic quenching due to high dE/dx
	  NumPhotons = BinomFluct(NumPhotons,biExc);
        
	// Next line commented by DF
	//NumPhotons = BinomFluct(NumPhotons,QE_EFF);
      } 
      NumElectrons = G4int(floor(NumElectrons+0.5));
	    
      //if ( SinglePhase ) //for a 1-phase det. don't propagate e-'s
      //NumElectrons = 0; //saves simulation time
	    
      // reset material properties numExc, numIon, numPho, numEle, as
      // their values have been used or stored elsewhere already
      aMaterialPropertiesTable->AddConstProperty( numExc, 0 );
      aMaterialPropertiesTable->AddConstProperty( numIon, 0 );
      aMaterialPropertiesTable->AddConstProperty( numPho, 0 );
      aMaterialPropertiesTable->AddConstProperty( numEle, 0 );
	    
      // start particle creation loop
      if( InitialKinetEnergy < MAX_ENE && InitialKinetEnergy > MIN_ENE && !fMultipleScattering )
	NumQuanta = NumPhotons + NumElectrons;
      else NumQuanta = 0;

      
      // DF: save the correspondent S1 and S2 energies
      DSEventHandler::Get()->SetS1Energy(
        DSEventHandler::Get()->GetS1Energy()+ (float(NumPhotons)/ScintillationYield)/keV 
      );
	
      DSEventHandler::Get()->SetS2Energy(
        DSEventHandler::Get()->GetS2Energy()+ (float(NumElectrons)/ScintillationYield)/keV 
      );

      
      // DF: kill S1 and S2 if required
      if(DSStorage::Get()->GetKillS1S2()) NumQuanta = 0 ;
      
      
      
      DSLog(debugging)<<" quanta "<<NumQuanta<<" photons "<<NumPhotons<<" electrons "<<NumElectrons<<endlog;      
      for(k = 0; k < NumQuanta; k++) {
	G4double sampledEnergy;
	G4DynamicParticle* aQuantum;
	      
	// Generate random direction
	G4double cost = 1. - 2.*G4UniformRand();
	G4double sint = std::sqrt((1.-cost)*(1.+cost));
	G4double phi = twopi*G4UniformRand();
	G4double sinp = std::sin(phi); G4double cosp = std::cos(phi);
	G4double px = sint*cosp; G4double py = sint*sinp;
	G4double pz = cost;
	      
	// Create momentum direction vector
	G4ParticleMomentum photonMomentum(px, py, pz);
	  
	  
	
	      
	// case of photon-specific stuff
	if (k < NumPhotons) {
	  // Determine polarization of new photon 
	  G4double sx = cost*cosp;
	  G4double sy = cost*sinp; 
	  G4double sz = -sint;
	  G4ThreeVector photonPolarization(sx, sy, sz);
	  G4ThreeVector perp = photonMomentum.cross(photonPolarization);
	  phi = twopi*G4UniformRand();
	  sinp = std::sin(phi);
	  cosp = std::cos(phi);
	  photonPolarization = cosp * photonPolarization + sinp * perp;
	  photonPolarization = photonPolarization.unit();
		
	  // Generate a new photon or electron:
	  sampledEnergy = G4RandGauss::shoot(PhotMean,PhotWidth);
	  aQuantum = new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(),photonMomentum);
	  aQuantum->SetPolarization(photonPolarization.x(),
				    photonPolarization.y(),
				    photonPolarization.z());
	}
	else // electrons drifted. Nout really used but defined to get the correct time
	  {
	    if(ElectricField) {
	      // point all electrons straight up, for drifting
	      G4ParticleMomentum electronMomentum(0, 0, -FieldSign);
	      aQuantum = new G4DynamicParticle(G4Electron::Electron(),electronMomentum);//useless
	      sampledEnergy = 0.5*EMASS*pow(GetLArDriftVelocity(Temperature,ElectricField)*cm/s,2.);
	    }
	    else {
	      // use "photonMomentum" for the electrons in the case of zero
	      // electric field, which is just randomized vector we made
	      aQuantum = new G4DynamicParticle(G4Electron::Electron(),photonMomentum);//useless
	      sampledEnergy = 1.38e-23*(joule/kelvin)*Temperature;
	    }
	  }
		
	//assign energy to make particle real
	aQuantum->SetKineticEnergy(sampledEnergy);
	DSLog(debugging) <<"sampledEnergy = " << sampledEnergy << endlog;
	
	
	// Generate new G4Track object:
	// emission time distribution
	
	// first an initial birth time is provided that is typically
	// <<1 ns after the initial interaction in the simulation, then
	// singlet, triplet lifetimes, and recombination time, are
	// handled here, to create a realistic S1 pulse shape/timing
	G4double aSecondaryTime = t0+G4UniformRand()*(t1-t0)+evtStrt;
	if ( aQuantum->GetDefinition()->GetParticleName()=="opticalphoton" ) {
	  if ( abs(z2-z1) && !fAlpha && abs(aParticle->GetPDGcode()) != 2112 ) {//electron recoil
	    LET = (energ/MeV)/(delta/cm)*(1/Density); //avg LET over all
	    //in future, this will be done interaction by interaction
	    // Next, find the recombination time, which is LET-dependent
	    // via ionization density (Kubota et al. Phys. Rev. B 20
	    // (1979) 3486). We find the LET-dependence by fitting to the
	    // E-dependence (Akimov et al. Phys. Lett. B 524 (2002) 245).
	   
	    //field dependence based on fitting Fig. 9 of Dawson et al.
	    //NIM A 545 (2005) 690
	    //singlet-triplet ratios adapted from Kubota 1979, converted
	    //into correct units to work here, and separately done for
	    //excitation and recombination processes for electron recoils
	    //and assumed same for all LET (may vary)
	    SingTripRatioX = G4RandGauss::shoot(0.17,0.05);
	    SingTripRatioR = G4RandGauss::shoot(0.8,0.2);
	    if ( z1 == 18 ) {
	      SingTripRatioR = 0.2701+0.003379*LET-4.7338e-5*pow(LET,2.)+8.1449e-6*pow(LET,3.); SingTripRatioX = SingTripRatioR;
	      if( LET < 3 ) {
		SingTripRatioX = G4RandGauss::shoot(0.36,0.06);
		SingTripRatioR = G4RandGauss::shoot(0.5,0.2); 
	      }
	    }
	  }//end electron recoil
	  else if ( fAlpha ) { //alpha particles
	    SingTripRatioR = G4RandGauss::shoot(2.3,0.51);
	    //currently based on Dawson 05 and Tey. 11 (arXiv:1103.3689)
	    //real ratio is likely a gentle function of LET
	    if (z1==18) SingTripRatioR = (-0.065492+1.9996*exp(-dE/MeV))/(1+0.082154/pow(dE/MeV,2.)) + 2.1811;
	    SingTripRatioX = SingTripRatioR;
	  }//end Alpha
	  else { //nuclear recoil
	    //based loosely on Hitachi et al. Phys. Rev. B 27 (1983) 5279
	    //with an eye to reproducing Akimov 2002 Fig. 9
	    SingTripRatioR = G4RandGauss::shoot(7.8,1.5);
	    if (z1==18) SingTripRatioR = 0.22218*pow(energ/keV,0.48211);
	    SingTripRatioX = SingTripRatioR;
	  }
	  // now, use binomial distributions to determine singlet and
	  // triplet states (and do separately for initially excited guys
	  // and recombining)
	  if ( k > NumExcitons ) {
	    //the recombination time is non-exponential, but approximates
	    //to exp at long timescales (see Kubota '79)
	    aSecondaryTime += tauR*(1./G4UniformRand()-1);
	    if(G4UniformRand()<SingTripRatioR/(1+SingTripRatioR))
	      aSecondaryTime -= tau1*log(G4UniformRand());
	    else aSecondaryTime -= tau3*log(G4UniformRand());
	  }
	  else {
	    if(G4UniformRand()<SingTripRatioX/(1+SingTripRatioX))
	      aSecondaryTime -= tau1*log(G4UniformRand());
	    else aSecondaryTime -= tau3*log(G4UniformRand());
	  }
	}

	// emission position distribution -- 
	// Generate the position of a new photon or electron, with NO
	// stochastic variation because that could lead to particles
	// being mistakenly generated outside of your active region by
	// Geant4, but real-life finite detector position resolution
	// wipes out any effects from here anyway...
	x0[0] = aMaterialPropertiesTable->GetConstProperty( xCoord );
	x0[1] = aMaterialPropertiesTable->GetConstProperty( yCoord );
	x0[2] = aMaterialPropertiesTable->GetConstProperty( zCoord );
	G4double radius = sqrt(pow(x0[0],2.)+pow(x0[1],2.));
	//re-scale radius to ensure no generation of quanta outside
	//the active volume of your simulation due to Geant4 rounding
	if ( radius >= R_TOL ) {
	  if (x0[0] == 0) x0[0] = 1*nm; 
	  if (x0[1] == 0) x0[1] = 1*nm;
	  radius -= R_TOL;
	  phi = atan ( x0[1] / x0[0] );
	  x0[0] = fabs(radius*cos(phi))*((fabs(x0[0]))/(x0[0]));
	  x0[1] = fabs(radius*sin(phi))*((fabs(x0[1]))/(x0[1]));
	}
	//position of the new secondary particle is ready for use
	G4ThreeVector aSecondaryPosition = x0;
	
	
	if ( k >= NumPhotons && diffusion && ElectricField > 0 ) {// simulate drift and e diffusion

	  //LAr
	  G4double D_T = 4.8*cm2/s; //ICARUS NIM A527 (2004) 329
	  G4double D_L = 18*cm2/s; //atrazhev  Timoshkim theory at 1kV/cm
	  if (ElectricField<100 && Phase == kStateLiquid) D_L = 8*cm2/s;
	  G4double vDrift=GetLArDriftVelocity(Temperature,ElectricField)*cm/s;
	  //G4cout<<" Dt "<<D_T<<" D_L "<<D_L<<" T "<<Temperature<<" vDrift "<<vDrift/(cm/s)<<" B "<<BORDER/cm<<" x0 "<<x0[2]/cm<<" dif "<<fabs(BORDER-x0[2])<<G4endl;
	  if ( BORDER == 0 ) x0[2] = 0;
	  G4double sigmaDT = sqrt(2*D_T*fabs(BORDER-x0[2])/vDrift);//cm
	  G4double sigmaDL = sqrt(2*D_L*fabs(BORDER-x0[2])/vDrift);//cm
	  G4double dr = fabs(G4RandGauss::shoot(0.,sigmaDT));
	  G4double tdrift=(fabs(BORDER-x0[2])/vDrift);
	  G4double dt = G4RandGauss::shoot(0.,(sigmaDL/vDrift));

	  phi = twopi * G4UniformRand();
	  aSecondaryPosition[0] += cos(phi) * dr;
	  aSecondaryPosition[1] += sin(phi) * dr;
	  aSecondaryPosition[2] = BORDER+1*um;
	  aSecondaryTime +=tdrift+dt;
	  radius = sqrt(pow(aSecondaryPosition[0],2.)+pow(aSecondaryPosition[1],2.));
	  
	  //kill thermal electrons according to purity...
	  double prob=exp(-(tdrift/ns)/(TaueLAr/ns));
	  double testvalue=G4UniformRand();
	  if(testvalue>prob)
	    continue;
	  else
	    if(!DSStorage::Get()->GetKillS2()) CreateS2(aSecondaryPosition,aSecondaryTime,&aParticleChange);
	  
	  continue;

	} //end of electron diffusion code
	
	
	// GEANT4 business: stuff you need to make a new track
	if ( aSecondaryTime < 0 ) aSecondaryTime = 0; //no neg. time
	G4Track * aSecondaryTrack = new G4Track(aQuantum,aSecondaryTime,aSecondaryPosition);
	if ( k < NumPhotons || radius < R_MAX )
	  aParticleChange.AddSecondary(aSecondaryTrack);
      }
	    
      //reset bunch of things when done with an interaction site
      aMaterialPropertiesTable->AddConstProperty( xCoord, 999*km );
      aMaterialPropertiesTable->AddConstProperty( yCoord, 999*km );
      aMaterialPropertiesTable->AddConstProperty( zCoord, 999*km );
      aMaterialPropertiesTable->AddConstProperty( trackL, 0*um );
      aMaterialPropertiesTable->AddConstProperty( energy, 0*eV );
      aMaterialPropertiesTable->AddConstProperty( time00, DBL_MAX );
      aMaterialPropertiesTable->AddConstProperty( time01, -1*ns );
      
      DSLog(debugging) << "\n Exiting from DSLight::DoIt -- "<< "NumberOfSecondaries = "<< aParticleChange.GetNumberOfSecondaries() << G4endl;
      
    } //end of interaction site loop

    //more things to reset...
    aMaterialPropertiesTable->AddConstProperty( "TOTALNUM_INT_SITES", 0 );
    aMaterialPropertiesTable->AddConstProperty( "ENERGY_DEPOSIT_TOT", 0*keV );
    aMaterialPropertiesTable->AddConstProperty( "ENERGY_DEPOSIT_GOL", 0*MeV );
    fExcitedNucleus = false;
    fAlpha = false;
  }
  
  //don't do anything when you're not ready to scintillate
  else {
    aParticleChange.SetNumberOfSecondaries(0);
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
	
  //the end (exiting)
  return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
}


//////////////////////////////
//  GetMeanFreePath
//////////////////////////////
G4double DSLight::GetMeanFreePath(const G4Track&,
				  G4double ,
				  G4ForceCondition* condition)
{
  *condition = StronglyForced;
  // what this does is enforce the DSLight physics process as always
  // happening, so in effect scintillation is a meta-process on top of
  // any and all other energy depositions which may occur, just like the
  // original G4Scintillation (disregard DBL_MAX, this function makes the
  // mean free path zero really, not infinite)
  
  return DBL_MAX; //a C-defined constant
}

//////////////////////////////
//  GetMeanLifeTime
//////////////////////////////
G4double DSLight::GetMeanLifeTime(const G4Track&,
				  G4ForceCondition* condition)
{
  *condition = Forced;
  // this function and this condition has the same effect as the above
  return DBL_MAX;
}

//////////////////////////////
//  GetLArDriftVelocity
//////////////////////////////
G4double DSLight::GetLArDriftVelocity(G4double tempinput, G4double efieldinput) {
  //from ATLAS internal note LARG-NO-058
  
  double p1=-0.016863;
  double p2 = -0.0083412;
  double p3 = 0.18088; 
  double p4 = 8.9751; 
  double p5 = 1.4614;
  double p6 = 0.32891;
  double T= tempinput;
  double T0 = 92.91;
  double E = efieldinput/1000.; //in kV/cm
  
  double vdrift = (p1*(T-T0)+1)*(p3*E*log(1+p4/E)+p5*pow(E,p6))+p2*(T-T0);
  
  //G4cout<<" test drift "<<vdrift<<G4endl;
  vdrift*=1e5;//cm/s
  return vdrift;
  
}


//////////////////////////////
//  CalculateElectronLET
//////////////////////////////
G4double DSLight::CalculateElectronLET ( G4double E) {
  G4double LET;  // in eV
  //LAr
  //use a spline fit to online ESTAR data
  //at energies <1 keV, use a different spline, determined manually by
  //generating sub-keV electrons in Geant4 and looking at their ranges, since
  //ESTAR does not go this low  

  if ( E >= 1 ) 
    LET = 116.70-162.97*log10(E)+99.361*pow(log10(E),2)-33.405*pow(log10(E),3)+6.5069*pow(log10(E),4)-0.69334*pow(log10(E),5)+.031563*pow(log10(E),6);
  else if ( E>0 && E<1 )
    LET = 100;
  else
    LET = 0;
  
  return LET;
}


/////////////////////////////
//        BinomFluct 
/////////////////////////////
G4int DSLight::BinomFluct ( G4int N0, G4double prob ) {
  G4double mean = N0*prob;
  G4double sigma = sqrt(N0*prob*(1-prob));
  G4int N1 = 0;
  if ( prob == 0.00 ) return N1;
  if ( prob == 1.00 ) return N0;
  
  if ( N0 < 10 ) {
    for(G4int i = 0; i < N0; i++) {
      if(G4UniformRand() < prob) N1++;
    }
  }
  else {
    N1 = G4int(floor(G4RandGauss::shoot(mean,sigma)+0.5));
  }
  if ( N1 > N0 ) N1 = N0;
  if ( N1 < 0 ) N1 = 0;
  return N1;
}

/////////////////////////////
///  InitMatPropValues
/////////////////////////////
void DSLight::InitMatPropValues ( G4MaterialPropertiesTable *nobleElementMat ) {
  char xCoord[80]; char yCoord[80]; char zCoord[80];
  char numExc[80]; char numIon[80]; char numPho[80]; char numEle[80];
  char trackL[80]; char time00[80]; char time01[80]; char energy[80];
  
  // for loop to initialize the interaction site mat'l properties
  for( G4int i=0; i<10000; i++ ) {
    sprintf(xCoord,"POS_X_%d",i); 
    sprintf(yCoord,"POS_Y_%d",i);
    sprintf(zCoord,"POS_Z_%d",i);
    nobleElementMat->AddConstProperty( xCoord, 999*km );
    nobleElementMat->AddConstProperty( yCoord, 999*km );
    nobleElementMat->AddConstProperty( zCoord, 999*km );
    sprintf(numExc,"N_EXC_%d",i); 
    sprintf(numIon,"N_ION_%d",i);
    sprintf(numPho,"N_PHO_%d",i); 
    sprintf(numEle,"N_ELE_%d",i);
    nobleElementMat->AddConstProperty( numExc, 0 );
    nobleElementMat->AddConstProperty( numIon, 0 );
    nobleElementMat->AddConstProperty( numPho, 0 );
    nobleElementMat->AddConstProperty( numEle, 0 );
    sprintf(trackL,"TRACK_%d",i); 
    sprintf(energy,"ENRGY_%d",i);
    sprintf(time00,"TIME0_%d",i); 
    sprintf(time01,"TIME1_%d",i);
    nobleElementMat->AddConstProperty( trackL, 0*um );
    nobleElementMat->AddConstProperty( energy, 0*eV );
    nobleElementMat->AddConstProperty( time00, DBL_MAX );
    nobleElementMat->AddConstProperty( time01,-1*ns );
  }
  
  // we initialize the total number of interaction sites, a variable for
  // updating the amount of energy deposited thus far in the medium, and a
  // variable for storing the amount of energy expected to be deposited
  nobleElementMat->AddConstProperty( "TOTALNUM_INT_SITES", 0 );
  nobleElementMat->AddConstProperty( "ENERGY_DEPOSIT_TOT", 0*keV );
  nobleElementMat->AddConstProperty( "ENERGY_DEPOSIT_GOL", 0*MeV );
  return;
}

//////////////////////////////
//       UnivScreenFunc
//////////////////////////////
G4double DSLight::UnivScreenFunc ( G4double E, G4double Z, G4double A ) {
  G4double a_0 = 5.29e-11*m; G4double a = 0.626*a_0*pow(Z,(-1./3.));
  G4double epsilon_0 = 8.854e-12*(farad/m);
  G4double epsilon = (a*E*2.*twopi*epsilon_0)/(2*pow(eplus,2.)*pow(Z,2.));
  G4double zeta_0 = pow(Z,(1./6.)); G4double m_N = A*1.66e-27*kg;
  G4double hbar = 6.582e-16*eV*s;
  if ( Z == 54 ) {
    epsilon *= 1.068; //zeta_0 = 1.63;
  } //special case for LXe from Bezrukov et al. 2011
  G4double s_n = log(1+1.1383*epsilon)/(2.*(epsilon +
                 0.01321*pow(epsilon,0.21226) +
                 0.19593*sqrt(epsilon)));
  G4double s_e = (a_0*zeta_0/a)*hbar*sqrt(8*epsilon*2.*twopi*epsilon_0/
                 (a*m_N*pow(eplus,2.)));
  return 1.38e5*0.5*(1+tanh(50*epsilon-0.25))*epsilon*(s_e/s_n);
}

//////////////////////////////
//       Create S2
//////////////////////////////
void DSLight::CreateS2 (G4ThreeVector position, G4double time, G4ParticleChange *apc){

  G4ThreeVector x1 = position;
  G4double      t1 = time;
  //set to gas
  const G4Material* aMaterial = DSMaterial::Get()->GetGaseousArgon();
  //G4double Density = aMaterial->GetDensity()/(g/cm3);
  G4double Pressure = aMaterial->GetPressure();
  //G4double nDensity = (Density/MolarMassLAr)*AVO;
  //G4int Phase = aMaterial->GetState();
  G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();
  
  if ( !YieldFactor || !YieldFactorS2 ){ 
    return;
  }

  G4double ElectricField;
  ElectricField = aMaterialPropertiesTable->GetConstProperty("ELECTRICFIELD");
  
   // in Ar at 3kV/cm all extracted //http://arxiv.org/pdf/astro-ph/0701286.pdf
  //otherwise a probability has to be set
  
  G4int numberofsecondary;//from http://arxiv.org/pdf/1207.2292.pdf
  G4double phgasA=0.0813; 
  G4double phgasB=139; 
  G4double phgasG=30.6; 
  
  G4double phpercm=phgasA*ElectricField/(volt/cm) - phgasB*Pressure/bar - phgasG;
  numberofsecondary=G4int(G4int(floor(YieldFactor))* G4RandGauss::shoot(phpercm*GASGAP/cm,sqrt(phpercm*GASGAP/cm)) + 0.5);
  numberofsecondary = G4int(numberofsecondary/DSStorage::Get()->GetScaleS2() + 0.5);
  
  DSLog(debugging)<< "\n Exiting from DS2Light::DoIt -- "<< "NumberOfSecondaries = "<<" "<<numberofsecondary<<endlog;
 
  for(int k = 0; k<numberofsecondary; k++) {	
    // start particle creation
    G4double sampledEnergy;
    G4DynamicParticle* aQuantum;
    
    // Generate random direction
    G4double cost = 1. - 2.*G4UniformRand();
    G4double sint = std::sqrt((1.-cost)*(1.+cost));
    G4double phi = twopi*G4UniformRand();
    G4double sinp = std::sin(phi);
    G4double cosp = std::cos(phi);
    G4double px = sint*cosp; G4double py = sint*sinp;
    G4double pz = cost;
    
    // Create momentum direction vector
    G4ParticleMomentum photonMomentum(px, py, pz);
    
    // Determine polarization of new photon
    G4double sx = cost*cosp; G4double sy = cost*sinp;
    G4double sz = -sint;
    G4ThreeVector photonPolarization(sx, sy, sz);
    G4ThreeVector perp = photonMomentum.cross(photonPolarization);
    phi = twopi*G4UniformRand();
    sinp = std::sin(phi); cosp = std::cos(phi);
    photonPolarization = cosp * photonPolarization + sinp * perp;
    photonPolarization = photonPolarization.unit();
    
    // Generate a new photon:
    sampledEnergy = G4RandGauss::shoot(E_eVLAr*eV,0.2*eV);
    aQuantum = new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(),photonMomentum);
    aQuantum->SetPolarization(photonPolarization.x(),
			      photonPolarization.y(),
			      photonPolarization.z());
    
    //assign energy to make particle real
    aQuantum->SetKineticEnergy(sampledEnergy);
    //G4cout << "sampledEnergy = " << sampledEnergy << G4endl;
    
    // electroluminesence emission time distribution
    G4double aSecondaryTime = t1, 
      SingTripRatio=.1; //guess: revisit!!!
    if(G4UniformRand()<SingTripRatio/(1+SingTripRatio))
      aSecondaryTime -= tau1LAr*log(G4UniformRand());
    else 
      aSecondaryTime -= tau3LAr*log(G4UniformRand());
    
    G4ThreeVector aSecondaryPosition = x1;
    aSecondaryPosition[2]+= (G4UniformRand()*GASGAP/mm) * mm;
    
    // GEANT4 business: stuff you need to make a new track
    G4Track* aSecondaryTrack = 
      new G4Track(aQuantum,aSecondaryTime,aSecondaryPosition);
    apc->AddSecondary(aSecondaryTrack);
    
  }
}

/*
 * $Log: DSLight.cc,v $
 * Revision 1.1  2014/05/07 12:21:03  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.10  2013/07/24 09:49:01  dfranco
 * Added S1 and S2 equivalent energies, the material index for each deposit, the command killS1S2 to kill photons and electrons generated by DSLight (after storing the equivalent energies)
 *
 * Revision 1.9  2013/07/05 09:09:08  meregaglia
 * tunded the DokeBirks[2] parameter and the light yield to match the 122 keV line in the DS10 paper
 *
 * Revision 1.8  2013/06/24 13:05:55  dfranco
 * TPC QE values were filled twice: once in the standard way, the second deconvoluting the reflections. The second has been commented
 *
 * Revision 1.7  2013/06/22 13:52:40  dfranco
 * Changed the Fano factor in DSLight from 0.07 (Doke 1976) to 0.20 (Doke 1984))
 *
 * Revision 1.6  2013/06/10 14:15:39  dfranco
 * Added two commands: /ds/physics/killS2 and /ds/physics/scaleS2 to kill or scale the S2 light
 *
 * Revision 1.5  2013/05/02 12:37:48  dfranco
 * Added ProcessType(fUserDefined) to DSLight to differentiate it from the Scintillation process
 *
 * Revision 1.4  2013/04/19 10:22:12  meregaglia
 * DSLight first major cleaning
 *
 * Revision 1.3  2013/04/18 13:55:18  dfranco
 * removed stupid anselmo messages
 *
 * Revision 1.2  2013/04/18 13:48:05  meregaglia
 * get rid of couts
 *
 * Revision 1.1  2013/04/18 12:55:43  meregaglia
 * S1S2 merged in DSLight
 *
 *
 */
