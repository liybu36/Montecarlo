#ifndef DSLight2_h
#define DSLight2_h 1

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


struct EnDepCluster {

  G4double       genPartPDG ;   // PDG of the particle generating the energy deposit
  G4double       genPartZ ;     // Z of the particle generating the energy deposit
  G4double       Length ; 
  G4double       Radius ;
  G4double       Energy ;       
  G4double       kinEne ;       // kinetic energy of the particle generating the energy deposit
  G4double       dEdx ;
  G4double       T0 ;
  G4ThreeVector  X0 ;
  G4ThreeVector  X1 ;
  G4ThreeVector  Baricenter ;
  G4int          nDeposits;

  G4int          nElectrons ;
  G4int          nPhotons ;
  G4int          nExcitons ;

};

struct S1S2Structure {

  long int      pdg ;
  double        t0 ;
  G4ThreeVector X0 ;
  G4ThreeVector X1 ;
  G4double      numCluster;
  G4double      depEne ;
  G4double      kEne ;
  G4double      dEdx ;
  G4int         z;
  G4bool        IsInTheBox ;
  G4int         nExcitons;
  G4int         nPhotons;
  G4int         nElectrons; 
  
};

class DSLight2 : public G4VRestDiscreteProcess {
private:
public: 
  
  DSLight2(const G4String& processName = "DSLight2",G4ProcessType type = fElectromagnetic);
  ~DSLight2();
  
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
  
  void SetScintillationYieldFactor(G4double val)            { fYieldFactorS1 =  val; } 
  G4double GetScintillationYieldFactor()                    { return fYieldFactorS1; }

  void SetScintillationYieldFactorS2(G4double val)          { fYieldFactorS2 =  val; } 
  G4double GetScintillationYieldFactorS2()                  { return fYieldFactorS2; } 
 
  void SetScintillationExcitationRatio(G4double val)        { fExcitationRatio = val; }
  G4double GetScintillationExcitationRatio()                { return fExcitationRatio; }


protected:

  int  fTracker[10000];
private:
  
  void ClearDeposits()                { fDeposits.clear();            }
  void SetDeposits(S1S2Structure val) { fDeposits.push_back(val);     }
  void SetDeposits()                  { fDeposits.push_back(fLocalDep);}
  
  void ReadData();
  G4int BinomFluct ( G4int, G4double);
  G4double CalculateElectronLET(G4double);
  
  G4double Interpolator(double, vector<float>&, vector<float>&);
  
  S1S2Structure              fLocalDep;
  EnDepCluster               fLocalClus;
  vector<S1S2Structure>      fDeposits;
  vector<EnDepCluster>       fDepClusters;

  G4double tauFast;
  G4double tauSlow;
  
  G4bool   fTrackSecondariesFirst; 
  G4double fYieldFactorS1; // turns scint. on/off
  G4double fYieldFactorS2; // turns S2 scint. on/off
  G4double fExcitationRatio; // N_ex/N_i, the dimensionless ratio of
  
  G4TrackVector*  trackVector;

  bool fEndOfEvent ;
  bool fStartSecondaries ;

  double fEne ;
    
  vector<float> fArgonEnergy;
  vector<float> fArgonELoss;
  
  vector<float> fElectronEnergy;
  vector<float> fElectronELoss;
  
  vector<float> fAlphaEnergy;
  vector<float> fAlphaELoss;
  


};



#endif 

/*
 * $Log: DSLight2.hh,v $
 * Revision 1.1  2014/05/07 12:20:53  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.7  2014/01/07 14:14:19  perassos
 * Bug fixed in the storage of the energy deposits and major changes in the application of the TI and DB models
 *
 * Revision 1.6  2013/08/16 15:33:51  perassos
 * Added S1 and S2 to DSLight2
 *
 * Revision 1.5  2013/08/06 13:58:18  dfranco
 * Added 3 variables to store the deposited energy in LAr, scintillator, and water. The last two are not yet implemented. g4rooter has been updated with 3 new variables: tpcene, vetoene, and muene
 *
 * Revision 1.4  2013/08/02 15:46:00  dfranco
 * Further development on DSLight2
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
