#ifndef DSLight_h
#define DSLight_h 1

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"


#define AVO 6.022e23 //Avogadro's number (#/mol)
#define EMASS 9.109e-31*kg
#define MilleDriftSpeed true
#define GASGAP 0.25*cm //S2 generation region
#define BORDER 27.467*cm //liquid-gas border z-coordinate
#define QE_EFF 1 //a base or maximum quantum efficiency


// this entire file is adapted from G4Scintillation.hh from Geant4.9.4
class DSLight : public G4VRestDiscreteProcess //class definition
{
  // Class inherits publicly from G4VRestDiscreteProcess
private:
public: // constructor and destructor
  
  DSLight(const G4String& processName = "DSLight",G4ProcessType type = fElectromagnetic);
  ~DSLight();
  
public: // methods, with descriptions
  G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
  // Returns true -> 'is applicable', for any particle type except for an
  // 'opticalphoton' and for short-lived particles
  
  G4double GetMeanFreePath(const G4Track& aTrack,G4double,G4ForceCondition* );
  // Returns infinity; i. e. the process does not limit the step, but 
  // sets the 'StronglyForced' condition for the DoIt to be invoked at
  // every step.
  
  G4double GetMeanLifeTime(const G4Track& aTrack,G4ForceCondition* );
  // Returns infinity; i. e. the process does not limit the time, but
  // sets the 'StronglyForced' condition for the DoIt to be invoked at
  // every step.
  
  // For in-flight particles losing energy (or those stopped)
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
				  const G4Step& aStep);
  G4VParticleChange* AtRestDoIt ( const G4Track& aTrack,
				  const G4Step& aStep);
  
  // These are the methods implementing the scintillation process.
  
  void SetTrackSecondariesFirst(const G4bool state);
  // If set, the primary particle tracking is interrupted and any
  // produced scintillation quanta are tracked next. When all have been
  // tracked, the tracking of the primary resumes.

  G4bool GetTrackSecondariesFirst() const;
  // Returns the boolean flag for tracking secondaries first.
  
  void SetScintillationYieldFactor(const G4double yieldfactor);
  // Called to set the scintillation quantum yield factor, useful for
  // shutting off scintillation entirely, or for producing a universal 
  // re-scaling to for example represent detector effects. Internally is
  // used for Lindhard yield factor for NR. Default should be user-set
  // to be 1 (for ER) in your simulation -- see NEST readme
  
  G4double GetScintillationYieldFactor() const;
  // Returns the quantum (photon/electron) yield factor. See above.

  G4double GetScintillationYieldFactorS2() const;
  void SetScintillationYieldFactorS2(const G4double yieldfactor);
  
  void SetScintillationExcitationRatio(const G4double excitationratio);
  // Called to set the scintillation exciton-to-ion ratio, useful for
  // when NumExcitons/NumIons is different for different types of parent
  // particles. This overwrites the ExcitationRatio obtained from the
  // G4MaterialPropertiesTable (e.g., 0.06 for LXe).

  G4double GetScintillationExcitationRatio() const;
  // Returns the ratio of the number of excitons to ions. Read above.


  void CreateS2(G4ThreeVector position, G4double time, G4ParticleChange* apc);
  
  G4double GetLArDriftVelocity(double T, double F);
  G4double CalculateElectronLET ( G4double E);

  G4double UnivScreenFunc ( G4double E, G4double Z, G4double A );

  G4int BinomFluct(G4int N0, G4double prob); //function for doing fluctuations

  void InitMatPropValues ( G4MaterialPropertiesTable* nobleElementMat );







protected:
  G4bool fTrackSecondariesFirst; // see above
  //bools for tracking some special particle cases
  G4bool fExcitedNucleus, fAlpha, fVeryHighEnergy, fMultipleScattering;
  G4double YieldFactor; // turns scint. on/off
  G4double YieldFactorS2; // turns S2 scint. on/off
  G4double ExcitationRatio; // N_ex/N_i, the dimensionless ratio of
  //initial excitons to ions
private:
  
};

////////////////////
// Inline methods
////////////////////
inline 
G4bool DSLight::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  if (aParticleType.GetParticleName() == "opticalphoton") return false;
  if (aParticleType.IsShortLived()) return false;
  //if(abs(aParticleType.GetPDGEncoding())==2112 || //neutron (no E-dep.)
  if(abs(aParticleType.GetPDGEncoding())==12 || //neutrinos (ditto) 
     abs(aParticleType.GetPDGEncoding())==14 ||
     abs(aParticleType.GetPDGEncoding())==16) return false;
  
  return true;
}

inline 
void DSLight::SetTrackSecondariesFirst(const G4bool state) 
{ 
  fTrackSecondariesFirst = state;
}

inline
G4bool DSLight::GetTrackSecondariesFirst() const
{
  return fTrackSecondariesFirst;
}

inline
void DSLight::SetScintillationYieldFactor(const G4double yieldfactor)
{
  YieldFactor = yieldfactor;
}

inline
G4double DSLight::GetScintillationYieldFactor() const
{
  return YieldFactor;
}

inline
void DSLight::SetScintillationYieldFactorS2(const G4double yieldfactor)
{
  YieldFactorS2 = yieldfactor;
}

inline
G4double DSLight::GetScintillationYieldFactorS2() const
{
  return YieldFactorS2;
}

inline
void DSLight::SetScintillationExcitationRatio(const G4double excitationratio)
{
  ExcitationRatio = excitationratio;
}

inline
G4double DSLight::GetScintillationExcitationRatio() const
{
  return ExcitationRatio;
}

#endif /* DSLight_h */

/*
 * $Log: DSLight.hh,v $
 * Revision 1.1  2014/05/07 12:20:53  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/04/19 10:22:12  meregaglia
 * DSLight first major cleaning
 *
 * Revision 1.1  2013/04/18 12:55:43  meregaglia
 * S1S2 merged in DSLight
 *
 * Revision 1.2  2013/04/11 14:13:48  meregaglia
 * fixed wrong gas z coordinate
 *
 * Revision 1.1  2013/04/10 21:06:06  meregaglia
 * added S1S2 from NEST and possibility to switch on or off this physics. Work in progress on the physics itsel
 *
 *
 */
