#ifndef _BXGENERATORAMBESOURCE_HH
#define _BXGENERATORAMBESOURCE_HH

//---------------------------------------------------------------------------//

#include "DSGeneratorAmBeSourceMessenger.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "DSVGenerator.hh"
#include "G4Event.hh"

//---------------------------------------------------------------------------//

class DSGeneratorAmBeSource : public DSVGenerator {
public:

  //default constructor
  DSGeneratorAmBeSource();

  //destructor
  virtual ~DSGeneratorAmBeSource();

  //public interface
  virtual void DSGeneratePrimaries(G4Event *event);
  
   void     SetNeutrinoType(G4int k)  { fNeutrinoType = k ;}
   G4int    GetNeutrinoType()         { return fNeutrinoType ;}

   void     SetDisableParticleType(G4int k)  { fDisableParticleType = k ;}
   G4int    GetDisableParticleType()         { return fDisableParticleType ;}

  
  //protected members
  enum Neutrinos
  {AllPN=0,NeutronNoG=1,Neutron1G=2,Neutron2G=3};


private:

  G4double ShootEnergyNeutron0G()  ;
  G4double ShootEnergyNeutron1G()  ;
  G4double ShootEnergyNeutron2G()  ;


  //fills the probability array and the energy bin array, called in the constructor
  void InitNeutron0G()  ;
  void InitNeutron1G()  ;
  void InitNeutron2G()  ;
  
  G4double Neutron0G_Probability [58];
  G4double Neutron0G_EnergyBin [58];
  G4double Neutron1G_Probability [58];
  G4double Neutron1G_EnergyBin [58];
  G4double Neutron2G_Probability [58];
  G4double Neutron2G_EnergyBin [58];

 
  DSGeneratorAmBeSourceMessenger*   fTheMessenger;

  G4ParticleTable*             fParticleTable;
  G4ParticleDefinition*        fParticle;
  G4ThreeVector                fPosition;
  G4ThreeVector                fDirection;

  G4int                        fNeutrinoType;
  G4int                        fDisableParticleType;
  G4bool                       isFirstTime;
};
#endif
