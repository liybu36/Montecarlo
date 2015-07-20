#ifndef DSLight3_h
#define DSLight3_h 1

#include <vector>
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

using namespace std ;



struct EnDepStructure {

  long int      pdg ;
  G4int         z;
  G4double      t0 ;
  G4ThreeVector X0 ;
  G4ThreeVector X1 ;
  G4double      depEne ;
  G4double      kEne ;
  G4double      dEdx ;

  G4int         nExcitons;
  G4int         nPhotons;
  G4int         nElectrons;

};



class DSLight3 : public G4VRestDiscreteProcess {
private:
public: 
  
  DSLight3(const G4String& processName = "DSLight3",G4ProcessType type = fElectromagnetic);
  ~DSLight3();
  
public: 
  G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
  G4double GetMeanFreePath(const G4Track& aTrack,G4double,G4ForceCondition* );
  G4double GetMeanLifeTime(const G4Track& aTrack,G4ForceCondition* );
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
  G4VParticleChange* AtRestDoIt ( const G4Track& aTrack, const G4Step& aStep);  
  void CreateS2 ( G4ThreeVector, G4double, G4ParticleChange* );
  G4double GetLArDriftVelocity( G4double, G4double );
  

  void SetTrackSecondariesFirst(G4bool val)                 { fTrackSecondariesFirst = val;  }
  G4bool GetTrackSecondariesFirst()                         { return fTrackSecondariesFirst; }
  
/*
  void SetScintillationYieldFactor(G4double val)            { fYieldFactorS1 =  val; } 
  G4double GetScintillationYieldFactor()                    { return fYieldFactorS1; }

  void SetScintillationYieldFactorS2(G4double val)          { fYieldFactorS2 =  val; } 
  G4double GetScintillationYieldFactorS2()                  { return fYieldFactorS2; } 

  void SetScintillationExcitationRatio(G4double val)        { fExcitationRatio = val; }
  G4double GetScintillationExcitationRatio()                { return fExcitationRatio; }
*/

protected:

  int  fTracker[10000];
private:
  
  
  void ReadData();
  G4int BinomFluct ( G4int, G4double);
  
  G4double Interpolator(double, vector<float>&, vector<float>&);
  
  G4double GetRecoProb(double, double) ;
  
  G4bool   fTrackSecondariesFirst; 
  
  G4bool   isFirstEnDepInLAr;

  G4double fRecoProb;
  G4double fQuenchingFactor;
  G4double fLightYield; 
  G4double fS2Yield ; 
  G4double fS1ScaleFactor ;
  G4double fMeanQuantaEnergy ; 
  
  G4double fExcitationRatioER;
  G4double fExcitationRatioNR;
  G4double fPhotonEneMean;
  G4double fPhotonEneWidth;

  G4double fLArGArBoundaryZ;
  G4double fDriftField;

  G4double fD_T;
  G4double fD_L;
  G4double fTauFast;
  G4double fTauSlow;
  G4double fTauReco;

  G4double fTaueLAr;

    
  vector<float> fArgonEnergy;
  vector<float> fArgonELoss;
  
  vector<float> fElectronEnergy;
  vector<float> fElectronELoss;
  
  vector<float> fAlphaEnergy;
  vector<float> fAlphaELoss;
  
  vector<float> fRecombinationEnergy;
  vector<float> fRecombinationProbability;


};



#endif 

/*
 * $Log: DSLight3.hh,v $
 * Revision 1.4  2015/01/20 09:22:05  dfranco
 * do not produce light for neutral particles in DSLight3 and new model for f90
 *
 * Revision 1.3  2014/12/22 14:40:49  dfranco
 * added the option to activate the recombination probability at 200 V/cm (/ds/physics/tunedS1); this option is by default true; selecting a specific drift field automatically switch off the tunedS1 option
 *
 * Revision 1.2  2014/07/23 14:55:39  pagnes
 * S2 first tuning
 *
 * Revision 1.1  2014/06/27 14:54:32  perassos
 * Implementation of a new scintillation model
 *
 * Revision 1.1  2014/05/07 12:20:53  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.7  2014/01/07 14:14:19  perassos
 * Bug fixed in the storage of the energy deposits and major changes in the application of the TI and DB models
 *
 * Revision 1.6  2013/08/16 15:33:51  perassos
 * Added S1 and S2 to DSLight3
 *
 * Revision 1.5  2013/08/06 13:58:18  dfranco
 * Added 3 variables to store the deposited energy in LAr, scintillator, and water. The last two are not yet implemented. g4rooter has been updated with 3 new variables: tpcene, vetoene, and muene
 *
 * Revision 1.4  2013/08/02 15:46:00  dfranco
 * Further development on DSLight3
 *
 * Revision 1.3  2013/08/02 12:25:24  dfranco
 * Development of new DSLight class
 *
 * Revision 1.2  2013/08/01 14:28:09  dfranco
 * added energy loss data from SRIM/TRIM
 *
 * Revision 1.1  2013/07/25 09:55:59  dfranco
 * Added a second version of NEST, still in development, and actually not working. The default version is DSLight
 *
 *
 */
