//---------------------------------------------------------------------------//
//    Adapted by davide.franco@mi.infn.it from the MaGe code:
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//                                                                           //
//                                                                           //
//                         MaGe Simulation                                   //
//                                                                           //
//      This code implementation is the intellectual property of the         //
//      MAJORANA and Gerda Collaborations. It is based on Geant4, an         //
//      intellectual property of the RD44 GEANT4 collaboration.              //
//                                                                           //
//                        *********************                              //
//                                                                           //
//    Neither the authors of this software system, nor their employing       //
//    institutes, nor the agencies providing financial support for this      //
//    work  make  any representation or  warranty, express or implied,       //
//    regarding this software system or assume any liability for its use.    //
//    By copying, distributing or modifying the Program (or any work based   //
//    on on the Program) you indicate your acceptance of this statement,     //
//    and all its terms.                                                     //
//                                                                           //
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//---------------------------------------------------------------------------//
#ifndef DSGeneratorCosmicRayMuons_h
#define DSGeneratorCosmicRayMuons_h 1

#include "DSVGenerator.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4Navigator.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"
#include "G4DataVector.hh"
#include "globals.hh"
#include "Randomize.hh"

#include "DSGeneratorCosmicRayMuonsMessenger.hh"

class DSGeneratorCosmicRayMuons : public DSVGenerator {

public:
  DSGeneratorCosmicRayMuons (); 
  virtual ~DSGeneratorCosmicRayMuons ();
  virtual void DSGeneratePrimaries(G4Event *event);

public:

  // position distribution  
  void SetHalfZ(G4double zhalf) {halfz = zhalf;};
  void SetRadius(G4double rad0) {Radius = rad0;};
  void SetSpectralIndex(G4double ind) {spectralIndex = ind;};
  void SetRockDepth(G4double rock) {rockDepth = rock;};
  void SetSamplerInitialized(G4bool samp) {samplerInitialized = samp;};
  void SetFileName(G4String stri) {fileName = stri;};
  void SetEnergySup(G4double ene) {energysup = ene;};
  void SetEnergyInf(G4double ene) {energyinf = ene;};
  G4String GetFileName() {return fileName;};

private: 
  // angular distribution
  void GenerateAngularSpectrum();

  // energy distribution 
  void ShootEnergy();

  // initial position
  void SampleInitialPosition();

public:
  //Get methods
  inline G4double GetParticleEnergy() {return particle_energy;}
  inline G4String GetParticleName() {return particle_definition->GetParticleName();};
  
  //Set methods
  void SetParticleDefinition(G4ParticleDefinition * aParticleDefinition);
  inline void SetParticleCharge(G4double aCharge)
  { particle_charge = aCharge; }
  
private:
  // position distribution
  G4double halfz;
  G4double Radius;

   // energy sampling
  static const G4int nbin=1000;
  G4double userpdf[nbin]; //nbin
  G4double spectralIndex;
  G4double rockDepth;
  G4bool samplerInitialized; 
  G4double energysup;
  G4double energyinf;

  //angular sampling
  G4DataVector* angularpdf;
  G4DataVector* costheta_pdf;
  G4DataVector* phi_pdf;
  G4int nx;
  G4int ny;
  G4String fileName; 
  G4bool fileFound;

  // particle properties 
  G4int                  NumberOfParticlesToBeGenerated;
  G4ParticleDefinition*  particle_definition;
  G4ParticleMomentum     particle_momentum_direction;
  G4double               particle_energy;
  G4double               particle_charge;
  G4ThreeVector          particle_position;
  G4double               particle_time;
  G4ThreeVector          particle_polarization;


  G4double MuonSpectrum(G4double ene);

private:
  DSGeneratorCosmicRayMuonsMessenger *fTheMessenger;
  G4Navigator *gNavigator;
  
};


#endif

/*
 * Revision 1.2  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 */
