#ifndef DSPhysicsList_h
#define DSPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"
#include "DSParameters.hh"

class DSPhysicsListMessenger;

class DSPhysicsList: public G4VUserPhysicsList
{
public:
  DSPhysicsList();
  ~DSPhysicsList();
  virtual void SetCuts();

  void SetHPRangeCuts(G4bool value)  { isHPRangeCuts = value; }
  G4bool GetHPRangeCuts()            { return isHPRangeCuts; }

  void SetDepositCuts(G4bool value)  { isDepositCuts = value; }
  G4bool GetDepositCuts()            { return isDepositCuts; }

  void SetS1S2(G4bool value)         { IsS1S2On = value;}
  G4bool GetS1S2()                   { return IsS1S2On; }

  void  SetHadronicList(G4int value) { fHadronic = value;}
  G4int GetHadronicList()            { return fHadronic; }

  void  SetEMList(G4int value)       { fEM = value;}
  G4int GetEMList()                  { return fEM; }

  void SetOptics(G4int  value)       { fOptics = value;}
  G4int  GetOptics()                 { return fOptics; }

 
protected:
  //! \brief Construct particle
  virtual void ConstructParticle();
  //! \brief Construct physics process
  virtual void ConstructProcess();
    
  // these methods Construct physics processes and register them
  //! \brief Construct general physics process and register it
  virtual void ConstructGeneral();
  virtual void ConstructOp();
  virtual void ConstructOp2();
  virtual void ConstructOp3();
  virtual void AddTransportation();
  
  

private:
  
  void ConstructEMStandard();
  void EnableAugerEffect();
  void ConstructEMLivermore();
  void ConstructIonPhysics();
  void DumpProcessTable();
  void ConstructHad_Custom() ;
  void ConstructHad_QGSP_BERT_HP() ;
  void ConstructHad_QSGP_BIC_HP() ;
  void ConstructHad_FTF_BIC_HP() ;
  void ConstructHad_FTFP_BERT_HP(); 
  void ConstructHad_Shielding(); 
  
  G4int VerboseLevel;
  G4int OpVerbLevel;

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  G4double cutForProton;
  G4double cutForAlpha;
  G4double cutForGenericIon;

  DSPhysicsListMessenger *fMessenger;

  G4bool isHPRangeCuts;
  G4bool isDepositCuts;

  G4bool IsS1S2On;
  G4int  fOptics;
  
  G4int  fHadronic;
  G4int  fEM;

  // these methods Construct particles 
  //! \brief these Construct bosons 
  void ConstructMyBosons();
  //! \brief these Construct leptons 
  void ConstructMyLeptons();
  //! \brief these Construct hadrons 
  void ConstructMyHadrons();
  //! \brief these Construct short-lived particles 
  void ConstructMyShortLiveds();
};

#endif
/*
 * $Log: DSPhysicsList.hh,v $
 * Revision 1.3  2015/04/29 14:49:17  dfranco
 * Fixed a bug in the DSGeneratorEnergyDeposit
 *
 * Revision 1.2  2014/06/27 14:54:33  perassos
 * Implementation of a new scintillation model
 *
 * Revision 1.1  2014/05/07 12:20:54  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.13  2014/01/21 10:50:53  perassos
 * Range cuts set to 1um in LAr and to 1 mm elsewhere
 *
 * Revision 1.12  2013/10/01 06:26:23  swesterd
 * added waveform averager
 *
 * Revision 1.11  2013/08/27 07:13:12  swesterd
 * add visible energy for the neutron veto
 *
 * Revision 1.10  2013/07/25 09:55:59  dfranco
 * Added a second version of NEST, still in development, and actually not working. The default version is DSLight
 *
 * Revision 1.9  2013/07/16 08:52:19  dfranco
 * added high precision (HP) neutron physics list for low energies
 *
 * Revision 1.8  2013/06/03 08:37:10  dfranco
 * added cvs logger
 *
 *
 */
