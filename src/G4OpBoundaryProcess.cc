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
////////////////////////////////////////////////////////////////////////
// Optical Photon Boundary Process Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpBoundaryProcess.cc
// Description: Discrete Process -- reflection/refraction at
//                                  optical interfaces
// Version:     1.1
// Created:     1997-06-18
// Modified:    1998-05-25 - Correct parallel component of polarization
//                           (thanks to: Stefano Magni + Giovanni Pieri)
//              1998-05-28 - NULL Rindex pointer before reuse
//                           (thanks to: Stefano Magni)
//              1998-06-11 - delete *sint1 in oblique reflection
//                           (thanks to: Giovanni Pieri)
//              1998-06-19 - move from GetLocalExitNormal() to the new 
//                           method: GetLocalExitNormal(&valid) to get
//                           the surface normal in all cases
//              1998-11-07 - NULL OpticalSurface pointer before use
//                           comparison not sharp for: std::abs(cost1) < 1.0
//                           remove sin1, sin2 in lines 556,567
//                           (thanks to Stefano Magni)
//              1999-10-10 - Accommodate changes done in DoAbsorption by
//                           changing logic in DielectricMetal
//              2001-10-18 - avoid Linux (gcc-2.95.2) warning about variables
//                           might be used uninitialized in this function
//                           moved E2_perp, E2_parl and E2_total out of 'if'
//              2003-11-27 - Modified line 168-9 to reflect changes made to
//                           G4OpticalSurface class ( by Fan Lei)
//              2004-02-02 - Set theStatus = Undefined at start of DoIt
//              2005-07-28 - add G4ProcessType to constructor
//              2006-11-04 - add capability of calculating the reflectivity
//                           off a metal surface by way of a complex index 
//                           of refraction - Thanks to Sehwook Lee and John 
//                           Hauptman (Dept. of Physics - Iowa State Univ.)
//              2009-11-10 - add capability of simulating surface reflections
//                           with Look-Up-Tables (LUT) containing measured
//                           optical reflectance for a variety of surface
//                           treatments - Thanks to Martin Janecek and
//                           William Moses (Lawrence Berkeley National Lab.)
//
// Author:      Peter Gumplinger
// 		adopted from work by Werner Keil - April 2/96
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4OpProcessSubType.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4GeometryTolerance.hh"
#include "G4SystemOfUnits.hh"
//// DS modification
#include "DSStorage.hh"
#include "DSParameters.hh"
#include <iostream>
#include "TMath.h"

using namespace std ;

/////////////////////////
// Class Implementation
/////////////////////////

        //////////////
        // Operators
        //////////////

// G4OpBoundaryProcess::operator=(const G4OpBoundaryProcess &right)
// {
// }

        /////////////////
        // Constructors
        /////////////////

G4OpBoundaryProcess::G4OpBoundaryProcess(const G4String& processName,
					 G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
{
  if ( verboseLevel > 0) {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
  
  SetProcessSubType(fOpBoundary);
  
  theStatus = Undefined;
  theModel = glisur;
  theFinish = polished;
  theReflectivity =  1.;
  theEfficiency   =  0.;
  theTransmittance = 0.;
  
  prob_sl = 0.;
  prob_ss = 0.;
  prob_bs = 0.;
  
  PropertyPointer  = NULL;
  PropertyPointer1 = NULL;
  PropertyPointer2 = NULL;
  
  Material1 = NULL;
  Material2 = NULL;
  
  OpticalSurface = NULL;
  
  kCarTolerance = G4GeometryTolerance::GetInstance()
    ->GetSurfaceTolerance();
  
  iTE = iTM = 0;
  thePhotonMomentum = 0.;
  Rindex1 = Rindex2 = cost1 = cost2 = sint1 = sint2 = 0.;

  isTransparent = false;
}

// G4OpBoundaryProcess::G4OpBoundaryProcess(const G4OpBoundaryProcess &right)
// {
// }

////////////////
// Destructors
////////////////

G4OpBoundaryProcess::~G4OpBoundaryProcess(){}

////////////
// Methods
////////////

// PostStepDoIt
// ------------
//
G4VParticleChange*
G4OpBoundaryProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  theStatus = Undefined;
  
  // DF: ITO variable   PA: grid variable
  G4bool IsITO = false ;
  G4bool IsGrid = false ;
  G4bool IsGridExit = false ;
  G4bool IsArTPB = false ;

  aParticleChange.Initialize(aTrack);
  aParticleChange.ProposeVelocity(aTrack.GetVelocity());
  
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

  if ( verboseLevel > 0 ) {
    G4cout << " Photon at Boundary! " << G4endl;
    G4VPhysicalVolume* thePrePV = pPreStepPoint->GetPhysicalVolume();
    G4VPhysicalVolume* thePostPV = pPostStepPoint->GetPhysicalVolume();
    if (thePrePV)  G4cout << " thePrePV:  " << thePrePV->GetName()  << G4endl;
    if (thePostPV) G4cout << " thePostPV: " << thePostPV->GetName() << G4endl;
  }
  
  if (pPostStepPoint->GetStepStatus() != fGeomBoundary){
    theStatus = NotAtBoundary;
    if ( verboseLevel > 0) BoundaryProcessVerbose();
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  if (aTrack.GetStepLength()<=kCarTolerance/2){
    theStatus = StepTooSmall;
    if ( verboseLevel > 0) BoundaryProcessVerbose();
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  
  Material1 = pPreStepPoint  -> GetMaterial();
  Material2 = pPostStepPoint -> GetMaterial();
  
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  
  thePhotonMomentum = aParticle->GetTotalMomentum();
  OldMomentum       = aParticle->GetMomentumDirection();
  OldPolarization   = aParticle->GetPolarization();
	
  if ( verboseLevel > 0 ) {
    G4cout << " Old Momentum Direction: " << OldMomentum     << G4endl;
    G4cout << " Old Polarization:       " << OldPolarization << G4endl;
  }
  
  G4ThreeVector theGlobalPoint = pPostStepPoint->GetPosition();
  
  G4Navigator* theNavigator =
    G4TransportationManager::GetTransportationManager()->
    GetNavigatorForTracking();
  
  G4bool valid;
  //  Use the new method for Exit Normal in global coordinates,
  //    which provides the normal more reliably. 
  theGlobalNormal = 
    theNavigator->GetGlobalExitNormal(theGlobalPoint,&valid);
  
  if (valid) {
    theGlobalNormal = -theGlobalNormal;
  }
  else 
    {
      G4ExceptionDescription ed;
      ed << " G4OpBoundaryProcess/PostStepDoIt(): "
	 << " The Navigator reports that it returned an invalid normal"
	 << G4endl;
      G4Exception("G4OpBoundaryProcess::PostStepDoIt", "OpBoun01",
		  EventMustBeAborted,ed,
		  "Invalid Surface Normal - Geometry must return valid surface normal");
    }
  
  if (OldMomentum * theGlobalNormal > 0.0) {
#ifdef G4OPTICAL_DEBUG
    G4ExceptionDescription ed;
    ed << " G4OpBoundaryProcess/PostStepDoIt(): "
       << " theGlobalNormal points in a wrong direction. "
       << G4endl;
    ed << "    The momentum of the photon arriving at interface (oldMomentum)"
       << " must exit the volume cross in the step. " << G4endl;
    ed << "  So it MUST have dot < 0 with the normal that Exits the new volume (globalNormal)." << G4endl;
    ed << "  >> The dot product of oldMomentum and global Normal is " << OldMomentum*theGlobalNormal << G4endl;
    ed << "     Old Momentum  (during step)     = " << OldMomentum << G4endl;
    ed << "     Global Normal (Exiting New Vol) = " << theGlobalNormal << G4endl;
    ed << G4endl;
    G4Exception("G4OpBoundaryProcess::PostStepDoIt", "OpBoun02",
		EventMustBeAborted,  // Or JustWarning to see if it happens repeatedbly on one ray
		ed,
		"Invalid Surface Normal - Geometry must return valid surface normal pointing in the right direction");
#else
    theGlobalNormal = -theGlobalNormal;
#endif
  }
  
  G4MaterialPropertiesTable* aMaterialPropertiesTable;
  G4MaterialPropertyVector* Rindex;

  // pdm
  //  verboseLevel = 1;
  
  aMaterialPropertiesTable = Material1->GetMaterialPropertiesTable();
  if (aMaterialPropertiesTable) {
    Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
  }
  else {
    theStatus = NoRINDEX;
    if ( verboseLevel > 0) BoundaryProcessVerbose();
    aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  
  if (Rindex) {
    Rindex1 = Rindex->Value(thePhotonMomentum);
  }
  else {
    theStatus = NoRINDEX;
    if ( verboseLevel > 0) BoundaryProcessVerbose();
    aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  
  theReflectivity =  1.;
  theEfficiency   =  0.;
  theTransmittance = 0.;
  
  theModel = glisur;
  theFinish = polished;
  
  G4SurfaceType type = dielectric_dielectric;
  
  Rindex = NULL;
  OpticalSurface = NULL;
  
  G4LogicalSurface* Surface = NULL;
  
  Surface = G4LogicalBorderSurface::GetSurface
    (pPreStepPoint ->GetPhysicalVolume(),
     pPostStepPoint->GetPhysicalVolume());
  
  if (Surface == NULL){
    G4bool enteredDaughter=(pPostStepPoint->GetPhysicalVolume()
			    ->GetMotherLogical() ==
			    pPreStepPoint->GetPhysicalVolume()
			    ->GetLogicalVolume());
    if(enteredDaughter){
      Surface = G4LogicalSkinSurface::GetSurface
	(pPostStepPoint->GetPhysicalVolume()->
	 GetLogicalVolume());
      if(Surface == NULL)
	Surface = G4LogicalSkinSurface::GetSurface
	  (pPreStepPoint->GetPhysicalVolume()->
	   GetLogicalVolume());
    }
    else {
      Surface = G4LogicalSkinSurface::GetSurface
	(pPreStepPoint->GetPhysicalVolume()->
	 GetLogicalVolume());
      if(Surface == NULL)
	Surface = G4LogicalSkinSurface::GetSurface
	  (pPostStepPoint->GetPhysicalVolume()->
	   GetLogicalVolume());
    }
  }
  
  if (Surface) 
    {
      OpticalSurface = 
	dynamic_cast <G4OpticalSurface*> (Surface->GetSurfaceProperty());
    }
  
  if (OpticalSurface) {
    
    type      = OpticalSurface->GetType();
    theModel  = OpticalSurface->GetModel();
    theFinish = OpticalSurface->GetFinish();
    
    aMaterialPropertiesTable = OpticalSurface->
      GetMaterialPropertiesTable();
    
    if (aMaterialPropertiesTable) {
      // SW: Allow boundaries to be transparent
      if(aMaterialPropertiesTable->ConstPropertyExists("DOTRANSMISSION"))
	isTransparent = true;
      else
	isTransparent = false;
      
      if (theFinish == polishedbackpainted ||
	  theFinish == groundbackpainted ) {
	Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
	if (Rindex) {
	  Rindex2 = Rindex->Value(thePhotonMomentum);
	}
	else {
	  theStatus = NoRINDEX;
	  if ( verboseLevel > 0) BoundaryProcessVerbose();
	  aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
	  aParticleChange.ProposeTrackStatus(fStopAndKill);
	  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}
      }
      
      PropertyPointer =
	aMaterialPropertiesTable->GetProperty("REFLECTIVITY");
      PropertyPointer1 =
	aMaterialPropertiesTable->GetProperty("REALRINDEX");
      PropertyPointer2 =
	aMaterialPropertiesTable->GetProperty("IMAGINARYRINDEX");
      
      iTE = 1;
      iTM = 1;
      
      if (PropertyPointer) {
	
                 theReflectivity =
		   PropertyPointer->Value(thePhotonMomentum);
		 
      } else if (PropertyPointer1 && PropertyPointer2) {
	
	CalculateReflectivity();
        
      }
      
      PropertyPointer =
	aMaterialPropertiesTable->GetProperty("EFFICIENCY");
      if (PropertyPointer) {
	theEfficiency =
	  PropertyPointer->Value(thePhotonMomentum);
      }
      
      PropertyPointer =
	aMaterialPropertiesTable->GetProperty("TRANSMITTANCE");
      if (PropertyPointer) {
	theTransmittance =
	  PropertyPointer->Value(thePhotonMomentum);
      }
      
      if ( theModel == unified ) {
	PropertyPointer =
	  aMaterialPropertiesTable->GetProperty("SPECULARLOBECONSTANT");
	if (PropertyPointer) {
	  prob_sl =
	    PropertyPointer->Value(thePhotonMomentum);
	} else {
	  prob_sl = 0.0;
	}
	
	PropertyPointer =
	  aMaterialPropertiesTable->GetProperty("SPECULARSPIKECONSTANT");
	if (PropertyPointer) {
	  prob_ss =
	    PropertyPointer->Value(thePhotonMomentum);
	} else {
	  prob_ss = 0.0;
	}
	
	PropertyPointer =
	  aMaterialPropertiesTable->GetProperty("BACKSCATTERCONSTANT");
	if (PropertyPointer) {
	  prob_bs =
	    PropertyPointer->Value(thePhotonMomentum);
	} else {
	  prob_bs = 0.0;
	}
      }
    }
    else if (theFinish == polishedbackpainted ||
	     theFinish == groundbackpainted ) {
      aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
      aParticleChange.ProposeTrackStatus(fStopAndKill);
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
  }
  
  if (type == dielectric_dielectric ) {
    if (theFinish == polished || theFinish == ground ) {
      
      if (Material1 == Material2){
	theStatus = SameMaterial;
	if ( verboseLevel > 0) BoundaryProcessVerbose();
	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
      }
      ///////////////////////////////////////////////////////////////////////
      /// Modfied by DF for adding ITO properties, assuming 15 nm thickness
      ///////////////////////////////////////////////////////////////////////
      if(aMaterialPropertiesTable->ConstPropertyExists("DOITO")) {

	//pdm       	std::cout << " Yo: DOITO  " << std::endl;

	G4double hc = 1239.84172; //eV*nm
	
       	theFacetNormal = theGlobalNormal;

	// PDM: print a table of ITO results for plotting.  Execute by running a single photon that hits ITO.

	//		G4double ITOth0 = 0.;
	//		do {
	//		  G4double *mytestITO = DSStorage::Get()->GetITOParameters(1.47, 1.233, 15., 420., ITOth0);
	//		  std::cout << " th0,R,T,A = " << ITOth0 <<" "<< mytestITO[0] <<" "<< mytestITO[1] <<" "<< mytestITO[2] <<std::endl;
	//		  ITOth0 += 0.01;
	//		} while (ITOth0 < TMath::Pi()/2.);

	// Single test value
	//	G4double *mytestITO = DSStorage::Get()->GetITOParameters(1.47, 1.233, 15., 420., 1.4);
	//	std::cout << " R,T,A = " << mytestITO[0] <<" "<< mytestITO[1] <<" "<< mytestITO[2] <<std::endl;

	G4double *myITO = DSStorage::Get()->GetITOParameters(
							     Material1->GetMaterialPropertiesTable()->GetProperty("RINDEX")->Value(thePhotonMomentum),
							     Material2->GetMaterialPropertiesTable()->GetProperty("RINDEX")->Value(thePhotonMomentum),
							     15.0, 
							     hc/(thePhotonMomentum/eV),
							     GetIncidentAngle()
							     );
	//pdm       	std::cout << " R,T,A = " << myITO[0] <<" "<< myITO[1] <<" "<< myITO[2] <<" sum= "<< myITO[0]+myITO[1]+myITO[2]<<std::endl;

	if(myITO[0] > 0) {
	  // myITO[2]: absorbance
	  // myITO[1]: transmittance
	  // myITO[0]: reflectivity	

	  theTransmittance =  myITO[1]/(myITO[0] + myITO[1]); 
	  theReflectivity  =  myITO[0]/(myITO[0] + myITO[1]);
	  theEfficiency    =  0;
	  
	  //theTransmittance =  myITO[1];
	  theReflectivity  =  myITO[0] + myITO[1];
	  
	  //if(G4BooleanRand(myITO[2]*2)) DoAbsorption();
	  IsITO = true ;
	  //return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}		
      }
      ////////////////////////////////////////////////////////////////////////////
      
      ///////////////////////////////////////////////////////////////////////
      /// Modfied by PA for adding new grid model (transmittance as function of incident angle)
      /// Updated by P. Meyers 5/21/2014
      ///////////////////////////////////////////////////////////////////////
      
      if(aMaterialPropertiesTable->ConstPropertyExists("DOGRID")) {	
	G4double hc = 1239.84172; //eV*nm
	//pdm       	std::cout << " New DOGRID  " << std::endl;

	// PDM 5/21/14: An update to the "New grid model."  The transparency is calculated in two parts:
	//  1) The transparency at normal incidence -- the hole/total area,
	//  2) An angle-dependent transparency of a hole due to the thickness of the grid.
	// Note that we violate usual G4 usage, where theReflectivity is really 1-absorption.
	// For the grid only, theReflectivity is the actual TRANSMISSION.  Sorry!

	G4double GridNormalTransp = DSParameters::Get()->GetGridNormalTransparency();

	// PDM 5/20/14: GetIncidentAngle uses theFacetNormal -- not yet defined.  For grid (with current model) it will be specular.
	//pdm	std::cout<<"***G4OBP:theGlobalNormal ="<<theGlobalNormal<<std::endl;

	theFacetNormal = theGlobalNormal;

	G4double myIncAngle = GetIncidentAngle();

       	theReflectivity = DSStorage::Get()->GetGridParameters(
							     Material1->GetMaterialPropertiesTable()->GetProperty("RINDEX")->Value(thePhotonMomentum),
							     Material2->GetMaterialPropertiesTable()->GetProperty("RINDEX")->Value(thePhotonMomentum),
							     20000.0, 
							     hc/(thePhotonMomentum/eV),
							     //							     GetIncidentAngle()
							     myIncAngle
							     );
	theReflectivity = theReflectivity * GridNormalTransp;
	IsGrid = true;

	//pdm       	std::cout<<"***G4OBP: thinc="<<myIncAngle*180./pi<<"  GridNormT="<< GridNormalTransp <<"  theReflectivity="<<theReflectivity<<std::endl;

	}
      if(aMaterialPropertiesTable->ConstPropertyExists("DOGRIDEXIT")) {	
	IsGridExit = true;
	//pdm       	std::cout << " DOGRIDEXIT  " << std::endl;   
 	}
           
      if(aMaterialPropertiesTable->ConstPropertyExists("DOArTPB")) {	
	IsArTPB = true;
	//pdm       	std::cout << " DOArTPB  " << std::endl;   
 	}
           
      aMaterialPropertiesTable =
	Material2->GetMaterialPropertiesTable();
      if (aMaterialPropertiesTable)
	Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
      if (Rindex) {
	Rindex2 = Rindex->Value(thePhotonMomentum);
      }
      else {
	theStatus = NoRINDEX;
	if ( verboseLevel > 0) BoundaryProcessVerbose();
	aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
	aParticleChange.ProposeTrackStatus(fStopAndKill);
	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
      }
    }
  }
  
  if (type == dielectric_metal) {
    
    DielectricMetal();
    
    // Uncomment the following lines if you wish to have 
    //         Transmission instead of Absorption
    if (isTransparent && theStatus == Absorption) {
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
    
  }
  else if (type == dielectric_LUT) {
    
    DielectricLUT();
    
  }
  else if (type == dielectric_dielectric) {
//    if(!IsITO) { // DF: ITO protection
    if(IsGridExit) {   // No action, all grid optics handled at grid entry.
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
      if ( theFinish == polishedbackpainted || theFinish == groundbackpainted ) {
	DielectricDielectric();
      }
      else {
	if (!G4BooleanRand(theReflectivity) ) {  // To G4, "!reflect" means "absorb". 
	  // For the grid, "theReflectivity" is the TRANSMISSION, and thus this branch does absorption and reflection.
	  // PDM 5/22/2014: only call if IsGrid; use LArGridVisRef/LArGridUVRef here, not GridReflection
	  if(IsGrid) {
	    //	    double refGrid = DSParameters::Get()->GetGridReflection();
	    //	    if (thePhotonMomentum < 8.3e-6 ) refGrid = .435;  
	    double refGrid = DSParameters::Get()->GetLArGridUVRef();
	    if (thePhotonMomentum < 8.3e-6 ) refGrid = DSParameters::Get()->GetLArGridVisRef();  
	    if (G4BooleanRand(refGrid)) {
	      //pdm	      std::cout << " GRID REFLECTION  " << std::endl;   
	      DoReflection();
	    } 
	    else {
	      //pdm	      std::cout << " GRID ABSORPTION  " << std::endl;   
	      DoAbsorption();
	    }

	  } 
	  else DoAbsorption();
	}
	else {
	  if ( theFinish == polishedfrontpainted ) {
                     DoReflection();
	  }
	  else if ( theFinish == groundfrontpainted ) {
	    theStatus = LambertianReflection;
	    DoReflection();
	  }
	  else if (IsGrid )  {   // PA trasnmittance trough the GRID
	    //pdm	    std::cout << " GRID TRANSMISSION  " << std::endl;   
	    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	  }
	  else if (IsArTPB )  {   // Taking over reflection/transmission at Ar (GAr or LAr) <--> TPB
	    if (!G4BooleanRand(theTransmittance) ) {  // Reflect
	      //pdm	      std::cout << " Ar-TPB reflection  " << std::endl;
	      theStatus = LambertianReflection;
	      DoReflection();
	    }
	    else {  // Model "Lambertian transmission".  Do Lambertian _reflection_, then reverse momentum and polarization.
	      theStatus = LambertianReflection;
	      DoReflection();
	      theStatus = LambertianTransmission;
	      NewMomentum = -NewMomentum;
	      NewPolarization = -NewPolarization;
	      //pdm	      std::cout << " Ar-TPB transmission  " << std::endl;
	    }
	  }
	  else DielectricDielectric();
	  }
	}
      
//    }
  }
  else {
    
    G4cerr << " Error: G4BoundaryProcess: illegal boundary type " << G4endl;
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    
  }
  
  NewMomentum = NewMomentum.unit();
  NewPolarization = NewPolarization.unit();
  
  if ( verboseLevel > 0) {
    G4cout << " New Momentum Direction: " << NewMomentum     << G4endl;
    G4cout << " New Polarization:       " << NewPolarization << G4endl;
    BoundaryProcessVerbose();
  }
  
  aParticleChange.ProposeMomentumDirection(NewMomentum);
  aParticleChange.ProposePolarization(NewPolarization);
  
  if ( theStatus == FresnelRefraction ) {
    G4MaterialPropertyVector* groupvel =
      Material2->GetMaterialPropertiesTable()->GetProperty("GROUPVEL");
    G4double finalVelocity = groupvel->Value(thePhotonMomentum);
    aParticleChange.ProposeVelocity(finalVelocity);
  }
  
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

void G4OpBoundaryProcess::BoundaryProcessVerbose() const
{
  if ( theStatus == Undefined )
    G4cout << " *** Undefined *** " << G4endl;
  if ( theStatus == FresnelRefraction )
    G4cout << " *** FresnelRefraction *** " << G4endl;
  if ( theStatus == FresnelReflection )
    G4cout << " *** FresnelReflection *** " << G4endl;
  if ( theStatus == TotalInternalReflection )
    G4cout << " *** TotalInternalReflection *** " << G4endl;
  if ( theStatus == LambertianReflection )
    G4cout << " *** LambertianReflection *** " << G4endl;
  if ( theStatus == LambertianTransmission )
    G4cout << " *** LambertianTransmission *** " << G4endl;
  if ( theStatus == LobeReflection )
    G4cout << " *** LobeReflection *** " << G4endl;
  if ( theStatus == SpikeReflection )
    G4cout << " *** SpikeReflection *** " << G4endl;
  if ( theStatus == BackScattering )
    G4cout << " *** BackScattering *** " << G4endl;
  if ( theStatus == PolishedLumirrorAirReflection )
    G4cout << " *** PolishedLumirrorAirReflection *** " << G4endl;
  if ( theStatus == PolishedLumirrorGlueReflection )
    G4cout << " *** PolishedLumirrorGlueReflection *** " << G4endl;
  if ( theStatus == PolishedAirReflection )
    G4cout << " *** PolishedAirReflection *** " << G4endl;
  if ( theStatus == PolishedTeflonAirReflection )
    G4cout << " *** PolishedTeflonAirReflection *** " << G4endl;
  if ( theStatus == PolishedTiOAirReflection )
    G4cout << " *** PolishedTiOAirReflection *** " << G4endl;
  if ( theStatus == PolishedTyvekAirReflection )
    G4cout << " *** PolishedTyvekAirReflection *** " << G4endl;
  if ( theStatus == PolishedVM2000AirReflection )
    G4cout << " *** PolishedVM2000AirReflection *** " << G4endl;
  if ( theStatus == PolishedVM2000GlueReflection )
    G4cout << " *** PolishedVM2000GlueReflection *** " << G4endl;
  if ( theStatus == EtchedLumirrorAirReflection )
    G4cout << " *** EtchedLumirrorAirReflection *** " << G4endl;
  if ( theStatus == EtchedLumirrorGlueReflection )
    G4cout << " *** EtchedLumirrorGlueReflection *** " << G4endl;
  if ( theStatus == EtchedAirReflection )
    G4cout << " *** EtchedAirReflection *** " << G4endl;
  if ( theStatus == EtchedTeflonAirReflection )
    G4cout << " *** EtchedTeflonAirReflection *** " << G4endl;
  if ( theStatus == EtchedTiOAirReflection )
    G4cout << " *** EtchedTiOAirReflection *** " << G4endl;
  if ( theStatus == EtchedTyvekAirReflection )
    G4cout << " *** EtchedTyvekAirReflection *** " << G4endl;
  if ( theStatus == EtchedVM2000AirReflection )
    G4cout << " *** EtchedVM2000AirReflection *** " << G4endl;
  if ( theStatus == EtchedVM2000GlueReflection )
    G4cout << " *** EtchedVM2000GlueReflection *** " << G4endl;
  if ( theStatus == GroundLumirrorAirReflection )
    G4cout << " *** GroundLumirrorAirReflection *** " << G4endl;
  if ( theStatus == GroundLumirrorGlueReflection )
    G4cout << " *** GroundLumirrorGlueReflection *** " << G4endl;
  if ( theStatus == GroundAirReflection )
    G4cout << " *** GroundAirReflection *** " << G4endl;
  if ( theStatus == GroundTeflonAirReflection )
    G4cout << " *** GroundTeflonAirReflection *** " << G4endl;
  if ( theStatus == GroundTiOAirReflection )
    G4cout << " *** GroundTiOAirReflection *** " << G4endl;
  if ( theStatus == GroundTyvekAirReflection )
    G4cout << " *** GroundTyvekAirReflection *** " << G4endl;
  if ( theStatus == GroundVM2000AirReflection )
    G4cout << " *** GroundVM2000AirReflection *** " << G4endl;
  if ( theStatus == GroundVM2000GlueReflection )
    G4cout << " *** GroundVM2000GlueReflection *** " << G4endl;
  if ( theStatus == Absorption )
    G4cout << " *** Absorption *** " << G4endl;
  if ( theStatus == Detection )
    G4cout << " *** Detection *** " << G4endl;
  if ( theStatus == NotAtBoundary )
    G4cout << " *** NotAtBoundary *** " << G4endl;
  if ( theStatus == SameMaterial )
    G4cout << " *** SameMaterial *** " << G4endl;
  if ( theStatus == StepTooSmall )
    G4cout << " *** StepTooSmall *** " << G4endl;
  if ( theStatus == NoRINDEX )
    G4cout << " *** NoRINDEX *** " << G4endl;
}

G4ThreeVector
G4OpBoundaryProcess::GetFacetNormal(const G4ThreeVector& Momentum,
			            const G4ThreeVector&  Normal ) const
{
  G4ThreeVector FacetNormal;
  
  if (theModel == unified || theModel == LUT) {
    
    /* This function code alpha to a random value taken from the
       distribution p(alpha) = g(alpha; 0, sigma_alpha)*std::sin(alpha),
       for alpha > 0 and alpha < 90, where g(alpha; 0, sigma_alpha)
       is a gaussian distribution with mean 0 and standard deviation
       sigma_alpha.  */
    
    G4double alpha;
    
    G4double sigma_alpha = 0.0;
    if (OpticalSurface) sigma_alpha = OpticalSurface->GetSigmaAlpha();
    
    if (sigma_alpha == 0.0) return FacetNormal = Normal;
    
    G4double f_max = std::min(1.0,4.*sigma_alpha);
    
    do {
      do {
	alpha = G4RandGauss::shoot(0.0,sigma_alpha);
      } while (G4UniformRand()*f_max > std::sin(alpha) || alpha >= halfpi );
      
      G4double phi = G4UniformRand()*twopi;
      
      G4double SinAlpha = std::sin(alpha);
      G4double CosAlpha = std::cos(alpha);
      G4double SinPhi = std::sin(phi);
      G4double CosPhi = std::cos(phi);
      
      G4double unit_x = SinAlpha * CosPhi;
      G4double unit_y = SinAlpha * SinPhi;
      G4double unit_z = CosAlpha;
      
      FacetNormal.setX(unit_x);
      FacetNormal.setY(unit_y);
      FacetNormal.setZ(unit_z);
      
      G4ThreeVector tmpNormal = Normal;
      
      FacetNormal.rotateUz(tmpNormal);
    } while (Momentum * FacetNormal >= 0.0);
  }
  else {
    
    G4double  polish = 1.0;
    if (OpticalSurface) polish = OpticalSurface->GetPolish();
    
    if (polish < 1.0) {
      do {
	G4ThreeVector smear;
	do {
	  smear.setX(2.*G4UniformRand()-1.0);
	  smear.setY(2.*G4UniformRand()-1.0);
	  smear.setZ(2.*G4UniformRand()-1.0);
	} while (smear.mag()>1.0);
	smear = (1.-polish) * smear;
	FacetNormal = Normal + smear;
              } while (Momentum * FacetNormal >= 0.0);
      FacetNormal = FacetNormal.unit();
    }
    else {
      FacetNormal = Normal;
    }
  }
  return FacetNormal;
}

void G4OpBoundaryProcess::DielectricMetal()
{
  G4int n = 0;
  
  do {
    
    n++;
    
    if( !G4BooleanRand(theReflectivity) && n == 1 ) {
      
      // Comment out DoAbsorption and uncomment theStatus = Absorption;
      // if you wish to have Transmission instead of Absorption

      if(isTransparent)
	theStatus = Absorption;
      else
	DoAbsorption();
      // theStatus = Absorption;
      break;
      
           }
    else {
      
      if (PropertyPointer1 && PropertyPointer2) {
	if ( n > 1 ) {
	  CalculateReflectivity();
	  if ( !G4BooleanRand(theReflectivity) ) {
	    DoAbsorption();
	    break;
	  }
	}
      }
      
      if ( theModel == glisur || theFinish == polished ) {
	
	DoReflection();
	
      } else {
	
	if ( n == 1 ) ChooseReflection();
        
	if ( theStatus == LambertianReflection ) {
	  DoReflection();
	}
	else if ( theStatus == BackScattering ) {
	  NewMomentum = -OldMomentum;
	  NewPolarization = -OldPolarization;
	}
	else {
	  
	  if(theStatus==LobeReflection){
	    if ( PropertyPointer1 && PropertyPointer2 ){
	    } else {
	      theFacetNormal =
		GetFacetNormal(OldMomentum,theGlobalNormal);
	    }
	  }
	  
	  G4double PdotN = OldMomentum * theFacetNormal;
	  NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
	  G4double EdotN = OldPolarization * theFacetNormal;
	  
	  G4ThreeVector A_trans, A_paral;
	  
	  if (sint1 > 0.0 ) {
	    A_trans = OldMomentum.cross(theFacetNormal);
	    A_trans = A_trans.unit();
	  } else {
	    A_trans  = OldPolarization;
	  }
	  A_paral   = NewMomentum.cross(A_trans);
	  A_paral   = A_paral.unit();
	  
	  if(iTE>0&&iTM>0) {
	    NewPolarization = 
	      -OldPolarization + (2.*EdotN)*theFacetNormal;
	  } else if (iTE>0) {
	    NewPolarization = -A_trans;
	  } else if (iTM>0) {
	    NewPolarization = -A_paral;
	  }
	  
	}
	
      }
      
      OldMomentum = NewMomentum;
      OldPolarization = NewPolarization;
      
    }

  } while (NewMomentum * theGlobalNormal < 0.0);
}

void G4OpBoundaryProcess::DielectricLUT()
{
  G4int thetaIndex, phiIndex;
  G4double AngularDistributionValue, thetaRad, phiRad, EdotN;
  G4ThreeVector PerpendicularVectorTheta, PerpendicularVectorPhi;
  
  theStatus = G4OpBoundaryProcessStatus(G4int(theFinish) + 
					(G4int(NoRINDEX)-G4int(groundbackpainted)));
  
  G4int thetaIndexMax = OpticalSurface->GetThetaIndexMax();
  G4int phiIndexMax   = OpticalSurface->GetPhiIndexMax();
  
  do {
           if ( !G4BooleanRand(theReflectivity) ) // Not reflected, so Absorbed
	     DoAbsorption();
           else {
	     // Calculate Angle between Normal and Photon Momentum
	     G4double anglePhotonToNormal = 
	       OldMomentum.angle(-theGlobalNormal);
	     // Round it to closest integer
	     G4int angleIncident = G4int(std::floor(180/pi*anglePhotonToNormal+0.5));

	     // Take random angles THETA and PHI, 
	     // and see if below Probability - if not - Redo
	     do {
	       thetaIndex = CLHEP::RandFlat::shootInt(thetaIndexMax-1);
	       phiIndex = CLHEP::RandFlat::shootInt(phiIndexMax-1);
	       // Find probability with the new indeces from LUT
	       AngularDistributionValue = OpticalSurface -> 
		 GetAngularDistributionValue(angleIncident,
					     thetaIndex,
					     phiIndex);
	     } while ( !G4BooleanRand(AngularDistributionValue) );
	     
	     thetaRad = (-90 + 4*thetaIndex)*pi/180;
	     phiRad = (-90 + 5*phiIndex)*pi/180;
	     // Rotate Photon Momentum in Theta, then in Phi
	     NewMomentum = -OldMomentum;
	     PerpendicularVectorTheta = NewMomentum.cross(theGlobalNormal);
	     if (PerpendicularVectorTheta.mag() > kCarTolerance ) {
	       PerpendicularVectorPhi = 
		 PerpendicularVectorTheta.cross(NewMomentum);
	     }
	     else {
	       PerpendicularVectorTheta = NewMomentum.orthogonal();
	       PerpendicularVectorPhi =
                                  PerpendicularVectorTheta.cross(NewMomentum);
	     }
	     NewMomentum =
	       NewMomentum.rotate(anglePhotonToNormal-thetaRad,
				  PerpendicularVectorTheta);
	     NewMomentum = NewMomentum.rotate(-phiRad,PerpendicularVectorPhi);
	     // Rotate Polarization too:
	     theFacetNormal = (NewMomentum - OldMomentum).unit();
	     EdotN = OldPolarization * theFacetNormal;
	     NewPolarization = -OldPolarization + (2.*EdotN)*theFacetNormal;
           }
  } while (NewMomentum * theGlobalNormal <= 0.0);
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

void G4OpBoundaryProcess::DielectricDielectric()
{
  G4bool Inside = false;
  G4bool Swap = false;

 leap:
  
  G4bool Through = false;
  G4bool Done = false;
  do {
    
    if (Through) {
      Swap = !Swap;
      Through = false;
      theGlobalNormal = -theGlobalNormal;
      G4SwapPtr(Material1,Material2);
      G4SwapObj(&Rindex1,&Rindex2);
    }
    
    if ( theFinish == polished ) {
      theFacetNormal = theGlobalNormal;
    }
    else {
      theFacetNormal =
	GetFacetNormal(OldMomentum,theGlobalNormal);
    }
    
    G4double PdotN = OldMomentum * theFacetNormal;
    G4double EdotN = OldPolarization * theFacetNormal;
    
    cost1 = - PdotN;
    if (std::abs(cost1) < 1.0-kCarTolerance){
      sint1 = std::sqrt(1.-cost1*cost1);
      sint2 = sint1*Rindex1/Rindex2;     // *** Snell's Law ***
    }
    else {
      sint1 = 0.0;
      sint2 = 0.0;
    }
    
    if (sint2 >= 1.0) {
      
      // Simulate total internal reflection
      
      if (Swap) Swap = !Swap;
      
      theStatus = TotalInternalReflection;
      
      if ( theModel == unified && theFinish != polished )
	ChooseReflection();
      
      if ( theStatus == LambertianReflection ) {
	DoReflection();
      }
      else if ( theStatus == BackScattering ) {
	NewMomentum = -OldMomentum;
	NewPolarization = -OldPolarization;
      }
      else {
	
	PdotN = OldMomentum * theFacetNormal;
	NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
	EdotN = OldPolarization * theFacetNormal;
	NewPolarization = -OldPolarization + (2.*EdotN)*theFacetNormal;
	
      }
    }
    else if (sint2 < 1.0) {
      
      // Calculate amplitude for transmission (Q = P x N)
      
      if (cost1 > 0.0) {
	cost2 =  std::sqrt(1.-sint2*sint2);
      }
      else {
	cost2 = -std::sqrt(1.-sint2*sint2);
      }
      
      G4ThreeVector A_trans, A_paral, E1pp, E1pl;
      G4double E1_perp, E1_parl;
      
      if (sint1 > 0.0) {
	A_trans = OldMomentum.cross(theFacetNormal);
	A_trans = A_trans.unit();
	E1_perp = OldPolarization * A_trans;
	E1pp    = E1_perp * A_trans;
	E1pl    = OldPolarization - E1pp;
	E1_parl = E1pl.mag();
      }
      else {
	A_trans  = OldPolarization;
	// Here we Follow Jackson's conventions and we set the
	// parallel component = 1 in case of a ray perpendicular
	// to the surface
	E1_perp  = 0.0;
	E1_parl  = 1.0;
      }
      
      G4double s1 = Rindex1*cost1;
      G4double E2_perp = 2.*s1*E1_perp/(Rindex1*cost1+Rindex2*cost2);
      G4double E2_parl = 2.*s1*E1_parl/(Rindex2*cost1+Rindex1*cost2);
      G4double E2_total = E2_perp*E2_perp + E2_parl*E2_parl;
      G4double s2 = Rindex2*cost2*E2_total;
      
      G4double TransCoeff;
      if (theTransmittance > 0 ) TransCoeff = theTransmittance;
      else if (cost1 != 0.0) TransCoeff = s2/s1;
      else TransCoeff = 0.0;
      
      G4double E2_abs, C_parl, C_perp;
      
      if ( !G4BooleanRand(TransCoeff) ) {
	
	// Simulate reflection
	if (Swap) Swap = !Swap;

	theStatus = FresnelReflection;

	if ( theModel == unified && theFinish != polished )
	  ChooseReflection();
	
	if ( theStatus == LambertianReflection ) {
	  DoReflection();
	}
	else if ( theStatus == BackScattering ) {
	  NewMomentum = -OldMomentum;
	  NewPolarization = -OldPolarization;
	}
	else {
	  
	  PdotN = OldMomentum * theFacetNormal;
	  NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
	  
	  if (sint1 > 0.0) {   // incident ray oblique

	    E2_parl   = Rindex2*E2_parl/Rindex1 - E1_parl;
	    E2_perp   = E2_perp - E1_perp;
	    E2_total  = E2_perp*E2_perp + E2_parl*E2_parl;
	    A_paral   = NewMomentum.cross(A_trans);
	    A_paral   = A_paral.unit();
	    E2_abs    = std::sqrt(E2_total);
	    C_parl    = E2_parl/E2_abs;
	    C_perp    = E2_perp/E2_abs;

	    NewPolarization = C_parl*A_paral + C_perp*A_trans;
	    
	    
	  }
	  
	  else {               // incident ray perpendicular
	    
	    if (Rindex2 > Rindex1) {
	      NewPolarization = - OldPolarization;
	               }
	    else {
	      NewPolarization =   OldPolarization;
	               }
	    
	  }
	}
      }
      
      else { // photon gets transmitted
	
	// Simulate transmission/refraction
	
	Inside = !Inside;
	Through = true;
	theStatus = FresnelRefraction;
	
	if (sint1 > 0.0) {      // incident ray oblique
	  
	  G4double alpha = cost1 - cost2*(Rindex2/Rindex1);
	  NewMomentum = OldMomentum + alpha*theFacetNormal;
	  NewMomentum = NewMomentum.unit();
	  PdotN = -cost2;
	  A_paral = NewMomentum.cross(A_trans);
	  A_paral = A_paral.unit();
	  E2_abs     = std::sqrt(E2_total);
	  C_parl     = E2_parl/E2_abs;
	  C_perp     = E2_perp/E2_abs;
	  
	  NewPolarization = C_parl*A_paral + C_perp*A_trans;
	  
	}
	else {                  // incident ray perpendicular
	  
	  NewMomentum = OldMomentum;
	  NewPolarization = OldPolarization;

	}
      }
    }
    
    OldMomentum = NewMomentum.unit();
    OldPolarization = NewPolarization.unit();
    
    if (theStatus == FresnelRefraction) {
      Done = (NewMomentum * theGlobalNormal <= 0.0);
    } 
    else {
      Done = (NewMomentum * theGlobalNormal >= 0.0);
    }
    
  } while (!Done);
  
  if (Inside && !Swap) {
    if( theFinish == polishedbackpainted ||
	theFinish == groundbackpainted ) {
      
      if( !G4BooleanRand(theReflectivity) ) {
	DoAbsorption();
      }
      else {
	if (theStatus != FresnelRefraction ) {
	  theGlobalNormal = -theGlobalNormal;
	}
	else {
	  Swap = !Swap;
	  G4SwapPtr(Material1,Material2);
	  G4SwapObj(&Rindex1,&Rindex2);
	}
	if ( theFinish == groundbackpainted )
	  theStatus = LambertianReflection;
	
	DoReflection();
	
	theGlobalNormal = -theGlobalNormal;
	OldMomentum = NewMomentum;
	
	goto leap;
      }
    }
  }
}

// GetMeanFreePath
// ---------------
//
G4double G4OpBoundaryProcess::GetMeanFreePath(const G4Track& ,
                                              G4double ,
                                              G4ForceCondition* condition)
{
  *condition = Forced;
  
  return DBL_MAX;
}

G4double G4OpBoundaryProcess::GetIncidentAngle() 
{
  G4double PdotN = OldMomentum * theFacetNormal;
  G4double magP= OldMomentum.mag();
  G4double magN= theFacetNormal.mag();
  G4double incidentangle = pi - std::acos(PdotN/(magP*magN));
  //pdm
  //  std::cout<<"G4OBP::GIA: OldMomentum"<<OldMomentum<<std::endl;
  //  std::cout<<"G4OBP::GIA: theFacetNormal"<<theFacetNormal<<std::endl;
  //  std::cout<<"G4OBP::GIA: PdotN, incidentangle "<<PdotN<<"\t"<<incidentangle<<std::endl;

  return incidentangle;
}

G4double G4OpBoundaryProcess::GetReflectivity(G4double E1_perp,
                                              G4double E1_parl,
                                              G4double incidentangle,
                                              G4double RealRindex,
                                              G4double ImaginaryRindex)
{
  
  G4complex Reflectivity, Reflectivity_TE, Reflectivity_TM;
  G4complex N(RealRindex, ImaginaryRindex);
  G4complex CosPhi;
  
  G4complex u(1,0);           //unit number 1
  
  G4complex numeratorTE;      // E1_perp=1 E1_parl=0 -> TE polarization
  G4complex numeratorTM;      // E1_parl=1 E1_perp=0 -> TM polarization
  G4complex denominatorTE, denominatorTM;
  G4complex rTM, rTE;
  
  // Following two equations, rTM and rTE, are from: "Introduction To Modern
  // Optics" written by Fowles
  
  CosPhi=std::sqrt(u-((std::sin(incidentangle)*std::sin(incidentangle))/(N*N)));
  
  numeratorTE   = std::cos(incidentangle) - N*CosPhi;
  denominatorTE = std::cos(incidentangle) + N*CosPhi;
  rTE = numeratorTE/denominatorTE;
  
  numeratorTM   = N*std::cos(incidentangle) - CosPhi;
  denominatorTM = N*std::cos(incidentangle) + CosPhi;
  rTM = numeratorTM/denominatorTM;
  
  // This is my calculaton for reflectivity on a metalic surface
  // depending on the fraction of TE and TM polarization
  // when TE polarization, E1_parl=0 and E1_perp=1, R=abs(rTE)^2 and
  // when TM polarization, E1_parl=1 and E1_perp=0, R=abs(rTM)^2
  
  Reflectivity_TE =  (rTE*conj(rTE))*(E1_perp*E1_perp)
    / (E1_perp*E1_perp + E1_parl*E1_parl);
  Reflectivity_TM =  (rTM*conj(rTM))*(E1_parl*E1_parl)
    / (E1_perp*E1_perp + E1_parl*E1_parl);
  Reflectivity    = Reflectivity_TE + Reflectivity_TM;
  
  do {
    if(G4UniformRand()*real(Reflectivity) > real(Reflectivity_TE))
      {iTE = -1;}else{iTE = 1;}
    if(G4UniformRand()*real(Reflectivity) > real(Reflectivity_TM))
      {iTM = -1;}else{iTM = 1;}
  } while(iTE<0&&iTM<0);
  
  return real(Reflectivity);
  
}

void G4OpBoundaryProcess::CalculateReflectivity()
{
  G4double RealRindex =
    PropertyPointer1->Value(thePhotonMomentum);
  G4double ImaginaryRindex =
    PropertyPointer2->Value(thePhotonMomentum);
  
  // calculate FacetNormal
  if ( theFinish == ground ) {
    theFacetNormal =
      GetFacetNormal(OldMomentum, theGlobalNormal);
  } else {
    theFacetNormal = theGlobalNormal;
  }
  
  G4double PdotN = OldMomentum * theFacetNormal;
  cost1 = -PdotN;
  
  if (std::abs(cost1) < 1.0 - kCarTolerance) {
    sint1 = std::sqrt(1. - cost1*cost1);
  } else {
    sint1 = 0.0;
  }
  
  G4ThreeVector A_trans, A_paral, E1pp, E1pl;
  G4double E1_perp, E1_parl;
  
  if (sint1 > 0.0 ) {
    A_trans = OldMomentum.cross(theFacetNormal);
    A_trans = A_trans.unit();
    E1_perp = OldPolarization * A_trans;
    E1pp    = E1_perp * A_trans;
    E1pl    = OldPolarization - E1pp;
    E1_parl = E1pl.mag();
  }
  else {
    A_trans  = OldPolarization;
    // Here we Follow Jackson's conventions and we set the
    // parallel component = 1 in case of a ray perpendicular
    // to the surface
    E1_perp  = 0.0;
    E1_parl  = 1.0;
  }
  
  //calculate incident angle
  G4double incidentangle = GetIncidentAngle();
  
  //calculate the reflectivity depending on incident angle,
  //polarization and complex refractive
  
  theReflectivity =
    GetReflectivity(E1_perp, E1_parl, incidentangle,
		    RealRindex, ImaginaryRindex);
}


/*
 * $Log: G4OpBoundaryProcess.cc,v $
 * Revision 1.4  2014/07/25 16:52:47  meyers
 * Add (commented-out) ITO diagnostics
 *
 * Revision 1.3  2014/07/02 17:07:21  meyers
 * Update model for LAr-TPB and GAr-TPB interface.  Now all Lambertian (reflection and transmission).
 *
 * Revision 1.2  2014/06/03 13:31:36  meyers
 * Migrate TPC grid and ITO optics updates to g4ds10
 *
 * Revision 1.7  2014/05/28 14:35:05  meyers
 * Set normal for ITO angle of incidence calculation
 *
 * Revision 1.6  2014/05/23 17:55:49  meyers
 * Update grid optical model
 *
 * Revision 1.5  2014/03/19 16:37:30  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.4  2014/03/11 09:51:23  dfranco
 * preliminary optical optimization
 *
 * Revision 1.3  2013/08/05 03:13:52  swesterd
 * some fine tuning of bscint and veto parameters
 *
 * Revision 1.2  2013/06/11 22:48:33  dfranco
 * Added ITO optical boundary. The Sernelius function is defined in DSStorage, and called by G4OpBoundaryProcess. New ITO bool variable added to DSDetectorDS50.cc as surface property (G4MaterialPropertyTable)
 *
 *
 */
