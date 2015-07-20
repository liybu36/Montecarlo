//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DSOpWLS.cc,v 1.1 2014/05/07 12:21:04 dfranco Exp $
// GEANT4 tag $Name:  $
//
////////////////////////////////////////////////////////////////////////
// Optical Photon WaveLength Shifting (WLS) Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpWLS.cc
// Description: Discrete Process -- Wavelength Shifting of Optical Photons
// Version:     1.0
// Created:     2003-05-13
// Author:      John Paul Archambault
//              (Adaptation of G4Scintillation and G4OpAbsorption)
// Updated:     2005-07-28 - add G4ProcessType to constructor
//              2006-05-07 - add G4VWLSTimeGeneratorProfile
// mail:        gum@triumf.ca
//              jparcham@phys.ualberta.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4OpProcessSubType.hh"

#include "DSOpWLS.hh"
#include "G4WLSTimeGeneratorProfileDelta.hh"
#include "G4WLSTimeGeneratorProfileExponential.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
/////////////////////////
// Class Implementation
/////////////////////////

/////////////////
// Constructors
/////////////////

DSOpWLS::DSOpWLS(const G4String& processName, G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
{
  SetProcessSubType(fOpWLS);

  theIntegralTable = 0;
 
  if (verboseLevel>0) {
    G4cout << GetProcessName() << " is created " << G4endl;
  }

  WLSTimeGeneratorProfile = 
       new G4WLSTimeGeneratorProfileDelta("WLSTimeGeneratorProfileDelta");

  BuildThePhysicsTable();
}

////////////////
// Destructors
////////////////

DSOpWLS::~DSOpWLS()
{
  if (theIntegralTable != 0) {
    theIntegralTable->clearAndDestroy();
    delete theIntegralTable;
  }
  delete WLSTimeGeneratorProfile;
}

////////////
// Methods
////////////

// PostStepDoIt
// -------------
//
G4VParticleChange*
DSOpWLS::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);
  
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  if (verboseLevel>0) {
    G4cout << "\n** Photon absorbed! **" << G4endl;
  }
  const G4Material* aMaterial = aTrack.GetMaterial();

  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
    
  G4MaterialPropertiesTable* aMaterialPropertiesTable =
    aMaterial->GetMaterialPropertiesTable();
  if (!aMaterialPropertiesTable)
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

  const G4MaterialPropertyVector* WLS_Intensity = 
    aMaterialPropertiesTable->GetProperty("WLSCOMPONENT"); 

  if (!WLS_Intensity)
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

  G4int NumPhotons = 1;

  if (aMaterialPropertiesTable->ConstPropertyExists("WLSMEANNUMBERPHOTONS")) {
     G4double MeanNumberOfPhotons = aMaterialPropertiesTable->
                                    GetConstProperty("WLSMEANNUMBERPHOTONS");

     
     ///////////////////////////////////////////////////////////
     // DF modification    
     // the WLS does not increase the photon. Either it absorbs or re-emit.
     // Hence, WLSMEANNUMBERPHOTONS becomes the efficiency in re-emitting the photon
     // next line is part of the original code
     //NumPhotons = G4int(G4Poisson(MeanNumberOfPhotons));
     
     if(MeanNumberOfPhotons > 1) MeanNumberOfPhotons = 1;
     NumPhotons = 0 ;
     if(G4UniformRand() < MeanNumberOfPhotons) NumPhotons = 1;
     
     ///////////////////////////////////////////////////////////
     
     if (NumPhotons <= 0) {

        // return unchanged particle and no secondaries

        aParticleChange.SetNumberOfSecondaries(0);

        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

     }
  }

  aParticleChange.SetNumberOfSecondaries(NumPhotons);

  G4int materialIndex = aMaterial->GetIndex();

  // Retrieve the WLS Integral for this material
  // new G4PhysicsOrderedFreeVector allocated to hold CII's

  G4double WLSTime = 0.*ns;
  G4PhysicsOrderedFreeVector* WLSIntegral = 0;

  WLSTime   = aMaterialPropertiesTable->
    GetConstProperty("WLSTIMECONSTANT");
  WLSIntegral =
    (G4PhysicsOrderedFreeVector*)((*theIntegralTable)(materialIndex));
   
  // Max WLS Integral
  
  G4double CIImax = WLSIntegral->GetMaxValue();
  
  for (G4int i = 0; i < NumPhotons; i++) {
    // First check the quantum yield of the WLS to see if this photon gets made
    G4double wlsefficiency = aMaterialPropertiesTable->GetConstProperty("WLSEFFICIENCY");
    if(!wlsefficiency) wlsefficiency = 1.;
    if(G4UniformRand() > wlsefficiency) continue;

    // Determine photon energy
    G4double oldEnergy = aTrack.GetKineticEnergy();
    G4double sampledEnergy = 0;
    G4double CIIvalue = 0;
    do
      {
        CIIvalue = G4UniformRand()*CIImax;
	sampledEnergy = 
	  WLSIntegral->GetEnergy(CIIvalue);
      }
    while(sampledEnergy > oldEnergy && false);
    
    //if (verboseLevel>1) {
    //  G4cout << "sampledEnergy = " << sampledEnergy << G4endl;
    //  G4cout << "CIIvalue =        " << CIIvalue << G4endl;
   // }
    
    // Generate random photon direction
    
    G4double cost = 1. - 2.*G4UniformRand();
    G4double sint = std::sqrt((1.-cost)*(1.+cost));

    G4double phi = twopi*G4UniformRand();
    G4double sinp = std::sin(phi);
    G4double cosp = std::cos(phi);
    
    G4double px = sint*cosp;
    G4double py = sint*sinp;
    G4double pz = cost;
    
    // Create photon momentum direction vector
    
    G4ParticleMomentum photonMomentum(px, py, pz);
    
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
    
    // Generate a new photon:
    
    G4DynamicParticle* aWLSPhoton =
      new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(),
			    photonMomentum);
    aWLSPhoton->SetPolarization
      (photonPolarization.x(),
       photonPolarization.y(),
       photonPolarization.z());
    
    aWLSPhoton->SetKineticEnergy(sampledEnergy);
    
    // Generate new G4Track object:
    
    // Must give position of WLS optical photon

    G4double TimeDelay = WLSTimeGeneratorProfile->GenerateTime(WLSTime);
    G4double aSecondaryTime = (pPostStepPoint->GetGlobalTime()) + TimeDelay;

    G4ThreeVector aSecondaryPosition = pPostStepPoint->GetPosition();

    G4Track* aSecondaryTrack = 
      new G4Track(aWLSPhoton,aSecondaryTime,aSecondaryPosition);
   
    aSecondaryTrack->SetTouchableHandle(aTrack.GetTouchableHandle()); 
    // aSecondaryTrack->SetTouchableHandle((G4VTouchable*)0);
    
    aSecondaryTrack->SetParentID(aTrack.GetTrackID());
    
    aParticleChange.AddSecondary(aSecondaryTrack);
  }
  if (verboseLevel>0) {
    G4cout << "\n Exiting from DSOpWLS::DoIt -- NumberOfSecondaries = " 
	   << aParticleChange.GetNumberOfSecondaries() << G4endl;  
  }
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// BuildThePhysicsTable for the wavelength shifting process
// --------------------------------------------------
//

void DSOpWLS::BuildThePhysicsTable()
{
  if (theIntegralTable) return;
  
  const G4MaterialTable* theMaterialTable = 
    G4Material::GetMaterialTable();
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();
  
  // create new physics table
  
  if(!theIntegralTable)theIntegralTable = new G4PhysicsTable(numOfMaterials);
  
  // loop for materials
  
  for (G4int i=0 ; i < numOfMaterials; i++)
    {
      G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
	new G4PhysicsOrderedFreeVector();
      
      // Retrieve vector of WLS wavelength intensity for
      // the material from the material's optical properties table.
      
      G4Material* aMaterial = (*theMaterialTable)[i];

      G4MaterialPropertiesTable* aMaterialPropertiesTable =
	aMaterial->GetMaterialPropertiesTable();

      if (aMaterialPropertiesTable) {

	G4MaterialPropertyVector* theWLSVector = 
	  aMaterialPropertiesTable->GetProperty("WLSCOMPONENT");

	if (theWLSVector) {
	  
	  // Retrieve the first intensity point in vector
	  // of (photon energy, intensity) pairs
	  
	  G4double currentIN = (*theWLSVector)[0];
	  
	  if (currentIN >= 0.0) {

	    // Create first (photon energy) 
	   
	    G4double currentPM = theWLSVector->Energy(0);
	    
	    G4double currentCII = 0.0;
	    
	    aPhysicsOrderedFreeVector->
	      InsertValues(currentPM , currentCII);
	    
	    // Set previous values to current ones prior to loop
	    
	    G4double prevPM  = currentPM;
	    G4double prevCII = currentCII;
	    G4double prevIN  = currentIN;
	    
	    // loop over all (photon energy, intensity)
	    // pairs stored for this material

            for (size_t i = 1;
                 i < theWLSVector->GetVectorLength();
                 i++)	    
	      {
		currentPM = theWLSVector->Energy(i);
		currentIN = (*theWLSVector)[i];
		
		currentCII = 0.5 * (prevIN + currentIN);
		
		currentCII = prevCII +
		  (currentPM - prevPM) * currentCII;
		
		aPhysicsOrderedFreeVector->
		  InsertValues(currentPM, currentCII);
		
		prevPM  = currentPM;
		prevCII = currentCII;
		prevIN  = currentIN;
	      }
	  }
	}
      }
	// The WLS integral for a given material
	// will be inserted in the table according to the
	// position of the material in the material table.

	theIntegralTable->insertAt(i,aPhysicsOrderedFreeVector);
    }
}

// GetMeanFreePath
// ---------------
//
G4double DSOpWLS::GetMeanFreePath(const G4Track& aTrack,
 				         G4double ,
				         G4ForceCondition* )
{
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4Material* aMaterial = aTrack.GetMaterial();

  G4double thePhotonEnergy = aParticle->GetTotalEnergy();

  G4MaterialPropertiesTable* aMaterialPropertyTable;
  G4MaterialPropertyVector* AttenuationLengthVector;
	
  G4double AttenuationLength = DBL_MAX;

  aMaterialPropertyTable = aMaterial->GetMaterialPropertiesTable();

  if ( aMaterialPropertyTable ) {
    AttenuationLengthVector = aMaterialPropertyTable->
      GetProperty("WLSABSLENGTH");
    if ( AttenuationLengthVector ){
      AttenuationLength = AttenuationLengthVector->
	Value(thePhotonEnergy);
    }
    else {
      //      G4cout << "No WLS absorption length vector specified for " << aMaterial->GetName() << G4endl;
    }
  }
  else {
    //    G4cout << "No WLS absortion length MPT specified" << G4endl;
  }
  return AttenuationLength;
}

void DSOpWLS::UseTimeProfile(const G4String name)
{
  if (name == "delta")
    {
      delete WLSTimeGeneratorProfile;
      WLSTimeGeneratorProfile = 
             new G4WLSTimeGeneratorProfileDelta("delta");
    }
  else if (name == "exponential")
    {
      delete WLSTimeGeneratorProfile;
      WLSTimeGeneratorProfile =
             new G4WLSTimeGeneratorProfileExponential("exponential");
    }
  else
    {
      G4Exception("DSOpWLS::UseTimeProfile", "em0202",
                  FatalException,
                  "generator does not exist");
    }
}
/*
 * $Log: DSOpWLS.cc,v $
 * Revision 1.1  2014/05/07 12:21:04  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.3  2013/06/22 11:43:15  dfranco
 * WLSMEANNUMBERPHOTONS defines now the efficiency in re-emitting the photon. 1 - WLSMEANNUMBERPHOTONS is the absorbtion probability
 *
 *
 */
