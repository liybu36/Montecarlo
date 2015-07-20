#include "DSLogger.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "DSPhysicsList.hh"
#include "DSPhysicsListMessenger.hh"

#include "globals.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4HadronCaptureProcess.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

#include "G4ios.hh"
#include <iomanip>

#include "G4UserLimits.hh"
#include "G4Version.hh"

// Constructor /////////////////////////////////////////////////////////////
DSPhysicsList::DSPhysicsList() : G4VUserPhysicsList() {

  
  defaultCutValue     =0.05*mm;    //0.5*mm; 
  cutForGamma	      = defaultCutValue;
  cutForElectron      = defaultCutValue; //1.0*nanometer;
  cutForPositron      = defaultCutValue;

  VerboseLevel    = 1;
  OpVerbLevel     = 0;

  fHadronic     =  1 ;
  fEM           =  2 ;
  fOptics       =  3 ;
  
  
  SetVerboseLevel(VerboseLevel);
  isHPRangeCuts =  true;
  isDepositCuts =  false;
  IsS1S2On      =  true; //default value of S1S2 scintillation physics OFF
   
  defaultCutValue     = 1.0*mm; //
  
  fMessenger = new DSPhysicsListMessenger(this);



}


// Destructor //////////////////////////////////////////////////////////////
DSPhysicsList::~DSPhysicsList()
{;}


// Construct Particles /////////////////////////////////////////////////////
#include "G4MesonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Geantino.hh"
#include "G4ChargedGeantino.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"


void DSPhysicsList::ConstructParticle()
{

  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

   
  //  mesons
  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();

  //  baryons
  G4BaryonConstructor bConstructor;
  bConstructor.ConstructParticle();

  //  ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();

  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  //  leptons
  G4LeptonConstructor lConstructor;
  lConstructor.ConstructParticle();
 
  //  bosons
    // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
  
  
  //G4BosonConstructor boConstructor;
  //boConstructor.ConstructParticle();
  
  // short lived
  if(fHadronic) {
    G4ShortLivedConstructor  sConstructor;
    sConstructor.ConstructParticle();
  }
  
}


// Construct Processes //////////////////////////////////////////////////////
#include "G4ProcessVector.hh"
void DSPhysicsList::ConstructProcess()
{

  AddTransportation();
       
       if(fEM == 1)           ConstructEMStandard();
  else if(fEM == 2)           ConstructEMLivermore();

       if(fOptics == 1)       ConstructOp();
  else if(fOptics == 2)       ConstructOp2();
  else if(fOptics == 3)       ConstructOp3();

       if(fHadronic == 0)     ConstructHad_Custom();  
  else if(fHadronic == 1)     ConstructHad_QGSP_BERT_HP();
  else if(fHadronic == 2)     ConstructHad_QSGP_BIC_HP();
  else if(fHadronic == 3)     ConstructHad_FTF_BIC_HP();
  else if(fHadronic == 4)     ConstructHad_FTFP_BERT_HP();
  else if(fHadronic == 5)     ConstructHad_Shielding();

  ConstructGeneral();
  
  ConstructIonPhysics();

  //  EnableAugerEffect();
  
  DumpProcessTable();
}  



// Transportation ///////////////////////////////////////////////////////////
void DSPhysicsList::AddTransportation() {
  G4VUserPhysicsList::AddTransportation();
}

// Standard EM  /////////////////////////////////////////////////////////////
#include "G4EmStandardPhysics.hh"
void DSPhysicsList::ConstructEMStandard()  {
  DSLog(routine) << "Standard EM Physics List Active"      << endlog;
  G4EmStandardPhysics *phys = new G4EmStandardPhysics(1);
  phys->ConstructProcess();
  
}
// Livermore EM  /////////////////////////////////////////////////////////////
#include "G4EmLivermorePhysics.hh"
void DSPhysicsList::ConstructEMLivermore()  {
  DSLog(routine) << "Livermore EM Physics List Active"      << endlog;
  G4EmLivermorePhysics *phys = new G4EmLivermorePhysics(1);
  phys->ConstructProcess();
}

// Enable Auger Effect //////////////////////////////////////////////////////
#include "G4EmProcessOptions.hh" 
void DSPhysicsList::EnableAugerEffect() {
  DSLog(routine) << "Enabling Auger Effect" << endlog;
  G4EmProcessOptions opt;
  opt.SetMscStepLimitation(fMinimal);
  opt.SetFluo(true);
  opt.SetAuger(true);
  opt.SetPIXE(true);
}







// IonPhysics  /////////////////////////////////////////////////////////////
#include "G4IonPhysics.hh"
void DSPhysicsList::ConstructIonPhysics()  {
  DSLog(routine) << "Ion Physics List Active"      << endlog;
  G4IonPhysics *phys = new G4IonPhysics ;
  phys->ConstructProcess();
}


// Optical Processes ////////////////////////////////////////////////////////
#include "G4Scintillation.hh"
#include "DSLight.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "DSOpWLS.hh"
#include "G4OpWLS.hh"
//#include "DSScintillation.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4EmSaturation.hh"
#include "G4LossTableManager.hh"

void DSPhysicsList::ConstructOp() {


   DSLight* theScintProcessDefS1S2 = new DSLight();
   theScintProcessDefS1S2->SetTrackSecondariesFirst(true);
   theScintProcessDefS1S2->SetScintillationYieldFactor(1.0); //
   theScintProcessDefS1S2->SetScintillationYieldFactorS2(1.0); //
   theScintProcessDefS1S2->SetVerboseLevel(OpVerbLevel);

   // scintillation process for alpha:
   DSLight* theScintProcessAlphaS1S2 = new DSLight();
   theScintProcessAlphaS1S2->SetTrackSecondariesFirst(true);
   theScintProcessAlphaS1S2->SetScintillationYieldFactor(1.0);
   theScintProcessAlphaS1S2->SetScintillationYieldFactorS2(1.0);
   theScintProcessAlphaS1S2->SetVerboseLevel(OpVerbLevel);

   // scintillation process for heavy nuclei
   DSLight* theScintProcessNucS1S2 = new DSLight();
   theScintProcessNucS1S2->SetTrackSecondariesFirst(true);
   theScintProcessNucS1S2->SetScintillationYieldFactor(1);
   theScintProcessNucS1S2->SetScintillationYieldFactorS2(1);
   theScintProcessNucS1S2->SetVerboseLevel(OpVerbLevel);


   // default scintillation process
   G4Scintillation* theScintProcessDef = new G4Scintillation("Scintillation");
   theScintProcessDef->SetTrackSecondariesFirst(true);
   theScintProcessDef->SetScintillationYieldFactor(1.0); //
   theScintProcessDef->SetScintillationExcitationRatio(1.0); //
   theScintProcessDef->SetVerboseLevel(OpVerbLevel);
   // Use Birks's correction in the scintillation process
   G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
   theScintProcessDef->AddSaturation(emSaturation);
   DSParameters::Get()->AddSaturation(emSaturation);
   //theScintProcessDef->DumpPhysicsTable();

   // optical processes
   G4OpAbsorption      *theAbsorptionProcess         = new G4OpAbsorption();
   G4OpRayleigh        *theRayleighScatteringProcess = new G4OpRayleigh();
   G4OpBoundaryProcess *theBoundaryProcess           = new G4OpBoundaryProcess();
   DSOpWLS             *theWLSProcess                = new DSOpWLS();
   //  theAbsorptionProcess->DumpPhysicsTable();
   //  theRayleighScatteringProcess->DumpPhysicsTable();
   theAbsorptionProcess->SetVerboseLevel(OpVerbLevel);
   // theRayleighScatteringProcess->SetVerboseLevel(OpVerbLevel);
   //theBoundaryProcess->SetVerboseLevel(OpVerbLevel);
   //G4OpticalSurfaceModel themodel = unified;
   //theBoundaryProcess->SetModel(themodel);
   theWLSProcess->UseTimeProfile("exponential");

   theParticleIterator->reset();
   while( (*theParticleIterator)() )
     {
       G4ParticleDefinition* particle = theParticleIterator->value();
       G4ProcessManager* pmanager = particle->GetProcessManager();
       G4String particleName = particle->GetParticleName();
       
       
       if (theScintProcessDefS1S2->IsApplicable(*particle)) {
	 //      if(particle->GetPDGMass() > 5.0*GeV)
	 if(particle->GetParticleName() == "GenericIon") {
	   pmanager->AddProcess(theScintProcessNucS1S2,ordDefault,ordInActive,ordDefault); // AtRestDiscrete
	   pmanager->SetProcessOrderingToLast(theScintProcessNucS1S2,idxAtRest);
	   pmanager->SetProcessOrderingToLast(theScintProcessNucS1S2,idxPostStep);
	 }
	 else if(particle->GetParticleName() == "alpha") {
	   pmanager->AddProcess(theScintProcessAlphaS1S2,ordDefault,ordInActive,ordDefault);
	   pmanager->SetProcessOrderingToLast(theScintProcessAlphaS1S2,idxAtRest);
	   pmanager->SetProcessOrderingToLast(theScintProcessAlphaS1S2,idxPostStep);
	 }
	 else {
	   pmanager->AddProcess(theScintProcessDefS1S2,ordDefault,ordInActive,ordDefault);
	   pmanager->SetProcessOrderingToLast(theScintProcessDefS1S2,idxAtRest);
	   pmanager->SetProcessOrderingToLast(theScintProcessDefS1S2,idxPostStep);
	 }
       }
       
       if (theScintProcessDef->IsApplicable(*particle)) {
	 
	 pmanager->AddProcess(theScintProcessDef,ordDefault,ordInActive,ordDefault);
	 pmanager->SetProcessOrderingToLast(theScintProcessDef,idxAtRest);
	 pmanager->SetProcessOrderingToLast(theScintProcessDef,idxPostStep);
	 
       }
       
       if (particleName == "opticalphoton") {
	 G4bool longRangeProcesses = true;
	 if(longRangeProcesses)
	   {
	     pmanager->AddDiscreteProcess(theAbsorptionProcess);
	     pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
	     pmanager->AddDiscreteProcess(theWLSProcess);
	   }
	 pmanager->AddDiscreteProcess(theBoundaryProcess);
       }
     }
     //theWLSProcess->DumpPhysicsTable();
}


// Optical Processes 2////////////////////////////////////////////////////////
#include "G4Scintillation.hh"
#include "DSLight2.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "DSOpWLS.hh"
#include "G4OpWLS.hh"
//#include "DSScintillation.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4EmSaturation.hh"
#include "G4LossTableManager.hh"

void DSPhysicsList::ConstructOp2() {


   DSLight2* theScintProcessDefS1S2 = new DSLight2();

   // scintillation process for alpha:
   DSLight2* theScintProcessAlphaS1S2 = new DSLight2();

   // scintillation process for heavy nuclei
   DSLight2* theScintProcessNucS1S2 = new DSLight2();


   // default scintillation process
   G4Scintillation* theScintProcessDef = new G4Scintillation("Scintillation");
   theScintProcessDef->SetTrackSecondariesFirst(true);
   theScintProcessDef->SetScintillationYieldFactor(1.0); //
   theScintProcessDef->SetScintillationExcitationRatio(1.0); //
   theScintProcessDef->SetVerboseLevel(OpVerbLevel);
   // Use Birks's correction in the scintillation process
   G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
   theScintProcessDef->AddSaturation(emSaturation);
   //theScintProcessDef->DumpPhysicsTable();

   // optical processes
   G4OpAbsorption      *theAbsorptionProcess         = new G4OpAbsorption();
   G4OpRayleigh        *theRayleighScatteringProcess = new G4OpRayleigh();
   G4OpBoundaryProcess *theBoundaryProcess           = new G4OpBoundaryProcess();
   DSOpWLS             *theWLSProcess                = new DSOpWLS();
   //  theAbsorptionProcess->DumpPhysicsTable();
   //  theRayleighScatteringProcess->DumpPhysicsTable();
   theAbsorptionProcess->SetVerboseLevel(OpVerbLevel);
   // theRayleighScatteringProcess->SetVerboseLevel(OpVerbLevel);
   //theBoundaryProcess->SetVerboseLevel(OpVerbLevel);
   //G4OpticalSurfaceModel themodel = unified;
   //theBoundaryProcess->SetModel(themodel);
   theWLSProcess->UseTimeProfile("exponential");

   theParticleIterator->reset();
   while( (*theParticleIterator)() )
     {
       G4ParticleDefinition* particle = theParticleIterator->value();
       G4ProcessManager* pmanager = particle->GetProcessManager();
       G4String particleName = particle->GetParticleName();
       
       if (theScintProcessDefS1S2->IsApplicable(*particle)) {
	 //      if(particle->GetPDGMass() > 5.0*GeV)
	 if(particle->GetParticleName() == "GenericIon") {
	   pmanager->AddProcess(theScintProcessNucS1S2,ordDefault,ordInActive,ordDefault); // AtRestDiscrete
	   pmanager->SetProcessOrderingToLast(theScintProcessNucS1S2,idxAtRest);
	   pmanager->SetProcessOrderingToLast(theScintProcessNucS1S2,idxPostStep);
	 }
	 else if(particle->GetParticleName() == "alpha") {
	   pmanager->AddProcess(theScintProcessAlphaS1S2,ordDefault,ordInActive,ordDefault);
	   pmanager->SetProcessOrderingToLast(theScintProcessAlphaS1S2,idxAtRest);
	   pmanager->SetProcessOrderingToLast(theScintProcessAlphaS1S2,idxPostStep);
	 }
	 else {
	   pmanager->AddProcess(theScintProcessDefS1S2,ordDefault,ordInActive,ordDefault);
	   pmanager->SetProcessOrderingToLast(theScintProcessDefS1S2,idxAtRest);
	   pmanager->SetProcessOrderingToLast(theScintProcessDefS1S2,idxPostStep);
	 }
       }
       if (theScintProcessDef->IsApplicable(*particle)) {
	 //      if(particle->GetPDGMass() > 5.0*GeV)
	 pmanager->AddProcess(theScintProcessDef,ordDefault,ordInActive,ordDefault);
	 pmanager->SetProcessOrderingToLast(theScintProcessDef,idxAtRest);
	 pmanager->SetProcessOrderingToLast(theScintProcessDef,idxPostStep);
       }
       
       if (particleName == "opticalphoton") {
	 pmanager->AddDiscreteProcess(theAbsorptionProcess);
	 pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
	 pmanager->AddDiscreteProcess(theWLSProcess);
	 pmanager->AddDiscreteProcess(theBoundaryProcess);
       }
     }
     //theWLSProcess->DumpPhysicsTable();
}



// Optical Processes 3////////////////////////////////////////////////////////
#include "G4Scintillation.hh"
#include "DSLight3.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "DSOpWLS.hh"
#include "G4OpWLS.hh"
//#include "DSScintillation.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4EmSaturation.hh"
#include "G4LossTableManager.hh"

void DSPhysicsList::ConstructOp3() {


   DSLight3* theScintProcessDefS1S2 = new DSLight3();

   // scintillation process for alpha:
   DSLight3* theScintProcessAlphaS1S2 = new DSLight3();

   // scintillation process for heavy nuclei
   DSLight3* theScintProcessNucS1S2 = new DSLight3();


   // default scintillation process
   G4Scintillation* theScintProcessDef = new G4Scintillation("Scintillation");
   theScintProcessDef->SetTrackSecondariesFirst(true);
//   theScintProcessDef->SetScintillationYieldFactor(1.0); //
//   theScintProcessDef->SetScintillationExcitationRatio(1.0); //
   theScintProcessDef->SetVerboseLevel(OpVerbLevel);
   // Use Birks's correction in the scintillation process
   G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
   theScintProcessDef->AddSaturation(emSaturation);
   //theScintProcessDef->DumpPhysicsTable();

   // optical processes
   G4OpAbsorption      *theAbsorptionProcess         = new G4OpAbsorption();
   G4OpRayleigh        *theRayleighScatteringProcess = new G4OpRayleigh();
   G4OpBoundaryProcess *theBoundaryProcess           = new G4OpBoundaryProcess();
   DSOpWLS             *theWLSProcess                = new DSOpWLS();
   //  theAbsorptionProcess->DumpPhysicsTable();
   //  theRayleighScatteringProcess->DumpPhysicsTable();
   theAbsorptionProcess->SetVerboseLevel(OpVerbLevel);
   // theRayleighScatteringProcess->SetVerboseLevel(OpVerbLevel);
   //theBoundaryProcess->SetVerboseLevel(OpVerbLevel);
   //G4OpticalSurfaceModel themodel = unified;
   //theBoundaryProcess->SetModel(themodel);
   theWLSProcess->UseTimeProfile("exponential");

   theParticleIterator->reset();
   while( (*theParticleIterator)() )
     {
       G4ParticleDefinition* particle = theParticleIterator->value();
       G4ProcessManager* pmanager = particle->GetProcessManager();
       G4String particleName = particle->GetParticleName();
       
       if (theScintProcessDefS1S2->IsApplicable(*particle)) {
	 //      if(particle->GetPDGMass() > 5.0*GeV)
	 if(particle->GetParticleName() == "GenericIon") {
	   pmanager->AddProcess(theScintProcessNucS1S2,ordDefault,ordInActive,ordDefault); // AtRestDiscrete
	   pmanager->SetProcessOrderingToLast(theScintProcessNucS1S2,idxAtRest);
	   pmanager->SetProcessOrderingToLast(theScintProcessNucS1S2,idxPostStep);
	 }
	 else if(particle->GetParticleName() == "alpha") {
	   pmanager->AddProcess(theScintProcessAlphaS1S2,ordDefault,ordInActive,ordDefault);
	   pmanager->SetProcessOrderingToLast(theScintProcessAlphaS1S2,idxAtRest);
	   pmanager->SetProcessOrderingToLast(theScintProcessAlphaS1S2,idxPostStep);
	 }
	 else {
	   pmanager->AddProcess(theScintProcessDefS1S2,ordDefault,ordInActive,ordDefault);
	   pmanager->SetProcessOrderingToLast(theScintProcessDefS1S2,idxAtRest);
	   pmanager->SetProcessOrderingToLast(theScintProcessDefS1S2,idxPostStep);
	 }
       }
       if (theScintProcessDef->IsApplicable(*particle)) {
	 //      if(particle->GetPDGMass() > 5.0*GeV)
	 pmanager->AddProcess(theScintProcessDef,ordDefault,ordInActive,ordDefault);
	 pmanager->SetProcessOrderingToLast(theScintProcessDef,idxAtRest);
	 pmanager->SetProcessOrderingToLast(theScintProcessDef,idxPostStep);
       }
       
       if (particleName == "opticalphoton") {
	 pmanager->AddDiscreteProcess(theAbsorptionProcess);
	 pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
	 pmanager->AddDiscreteProcess(theWLSProcess);
	 pmanager->AddDiscreteProcess(theBoundaryProcess);
       }
     }
     //theWLSProcess->DumpPhysicsTable();
}




//////////////////////////////////////////////////////////////////////////////
// Hadronic processes ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


//Underground shielding
#include "G4HadronPhysicsShielding.hh"
void DSPhysicsList::ConstructHad_Shielding() {
  DSLog(routine) << "Hadronic Physics Active: Underground Shielding model" << endlog;
  G4HadronPhysicsShielding *hadPhysicsList = new G4HadronPhysicsShielding("Shielding",true);
  //HadronPhysicsShielding::useLEND("ss");
  hadPhysicsList->ConstructProcess();
   
}

//QGSP BERTINI cascade NEUTRONHP model
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
void DSPhysicsList::ConstructHad_QGSP_BERT_HP() {
  DSLog(routine) << "Hadronic Physics Active: QGSP_BERT_HP model" << endlog;
  G4HadronPhysicsQGSP_BERT_HP *hadPhysicsList = new G4HadronPhysicsQGSP_BERT_HP;
  hadPhysicsList->ConstructProcess();
   
}


//QGSP BINARY cascade NEUTRONHP model
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
void DSPhysicsList::ConstructHad_QSGP_BIC_HP() {
  DSLog(routine) << "Hadronic Physics Active: QSGP_BIC_HP model" << endlog;
  G4HadronPhysicsQGSP_BIC_HP *hadPhysicsList = new G4HadronPhysicsQGSP_BIC_HP;
  hadPhysicsList->ConstructProcess();
   
}

//FTF BINARY cascade NEUTRONHP model
#include "G4HadronPhysicsFTF_BIC.hh"
void DSPhysicsList::ConstructHad_FTF_BIC_HP() {
  DSLog(routine) << "Hadronic Physics Active: FTF_BIC_HP model" << endlog;
  G4HadronPhysicsFTF_BIC *hadPhysicsList = new G4HadronPhysicsFTF_BIC;
  hadPhysicsList->ConstructProcess();
   
}

//FTFP BERTINI cascade NEUTRONHP model
#include "G4HadronPhysicsFTFP_BERT.hh"
void DSPhysicsList::ConstructHad_FTFP_BERT_HP() {
  DSLog(routine) << "Hadronic Physics Active: FTFP_BERT_HP model" << endlog;
  G4HadronPhysicsFTFP_BERT *hadPhysicsList = new G4HadronPhysicsFTFP_BERT;
  hadPhysicsList->ConstructProcess();  
}



//Custom model
// Elastic processes:
#include "G4HadronElasticProcess.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"

// Inelastic processes:
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

//migration
// High energy FTFP model and Bertini cascade
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4CascadeInterface.hh"

// Cross sections
#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4CrossSectionElastic.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4AntiNuclElastic.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4GGNuclNuclCrossSection.hh"
#include "G4HadronElastic.hh"
#include "G4HadronCaptureProcess.hh"

// Neutron high-precision models: <20 MeV
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"
//migration
//#include "G4LCapture.hh"

// Stopping processes
// migration
#include "G4PiMinusAbsorptionBertini.hh"
#include "G4KaonMinusAbsorptionBertini.hh"
#include "G4AntiProtonAbsorptionFritiof.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"


// ConstructHad()
// Makes discrete physics processes for the hadrons, at present limited
// to those particles with GHEISHA interactions (INTRC > 0).
// The processes are: Elastic scattering and Inelastic scattering.
// F.W.Jones  09-JUL-1998
void DSPhysicsList::ConstructHad_Custom() 
{
//migration
//  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
//  G4LElastic* theElasticModel = new G4LElastic;
//  theElasticProcess->RegisterMe(theElasticModel);
  
    //Elastic models 
  const G4double elastic_elimitPi = 1.0*GeV;

  G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
  G4HadronElastic* elastic_lhep1 = new G4HadronElastic();
  elastic_lhep1->SetMaxEnergy( elastic_elimitPi );
  G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE(); 
  elastic_he->SetMinEnergy( elastic_elimitPi );

  
  // Inelastic scattering
  const G4double theFTFMin0 =    0.0*GeV;
  const G4double theFTFMin1 =    4.0*GeV;
  const G4double theFTFMax =   100.0*TeV;
  const G4double theBERTMin0 =   0.0*GeV;
  const G4double theBERTMin1 =  19.0*MeV;
  const G4double theBERTMax =    5.0*GeV;
  const G4double theHPMin =      0.0*GeV;
  const G4double theHPMax =     20.0*MeV;

  G4FTFModel * theStringModel = new G4FTFModel;
  G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay( new G4LundStringFragmentation );
  theStringModel->SetFragmentationModel( theStringDecay );
  G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel( new G4ExcitationHandler );
  G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface( thePreEquilib );

  G4TheoFSGenerator * theFTFModel0 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel0->SetHighEnergyGenerator( theStringModel );
  theFTFModel0->SetTransport( theCascade );
  theFTFModel0->SetMinEnergy( theFTFMin0 );
  theFTFModel0->SetMaxEnergy( theFTFMax );

  G4TheoFSGenerator * theFTFModel1 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel1->SetHighEnergyGenerator( theStringModel );
  theFTFModel1->SetTransport( theCascade );
  theFTFModel1->SetMinEnergy( theFTFMin1 );
  theFTFModel1->SetMaxEnergy( theFTFMax );

  G4CascadeInterface * theBERTModel0 = new G4CascadeInterface;
  theBERTModel0->SetMinEnergy( theBERTMin0 );
  theBERTModel0->SetMaxEnergy( theBERTMax );

  G4CascadeInterface * theBERTModel1 = new G4CascadeInterface;
  theBERTModel1->SetMinEnergy( theBERTMin1 );
  theBERTModel1->SetMaxEnergy( theBERTMax );

  G4VCrossSectionDataSet * thePiData = new G4CrossSectionPairGG( new G4PiNuclearCrossSection, 91*GeV );
  G4VCrossSectionDataSet * theAntiNucleonData = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
  G4VCrossSectionDataSet * theGGNuclNuclData = G4CrossSectionDataSetRegistry::Instance()->
    GetCrossSectionDataSet(G4GGNuclNuclCrossSection::Default_Name());
  
  
  theParticleIterator->reset();
  while ((*theParticleIterator)()) 
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "pi+") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_lhep1 );
          theElasticProcess->RegisterMe( elastic_he );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  //Inelastic scattering
	  G4PionPlusInelasticProcess* theInelasticProcess = 
	    new G4PionPlusInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( thePiData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	} 

      else if (particleName == "pi-") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_lhep1 );
          theElasticProcess->RegisterMe( elastic_he );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  //Inelastic scattering
	  G4PionMinusInelasticProcess* theInelasticProcess = 
	    new G4PionMinusInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( thePiData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
	  //Absorption
	  pmanager->AddRestProcess(new G4PiMinusAbsorptionBertini, ordDefault);
	}
      
      else if (particleName == "kaon+") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering	
	  G4KaonPlusInelasticProcess* theInelasticProcess = 
	    new G4KaonPlusInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      
      else if (particleName == "kaon0S") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering	 
	  G4KaonZeroSInelasticProcess* theInelasticProcess = 
	    new G4KaonZeroSInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
	}

      else if (particleName == "kaon0L") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  // Inelastic scattering
	  G4KaonZeroLInelasticProcess* theInelasticProcess = 
	    new G4KaonZeroLInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 ); 
	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
	}

      else if (particleName == "kaon-") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4KaonMinusInelasticProcess* theInelasticProcess = 
	    new G4KaonMinusInelasticProcess("inelastic");	
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
					   GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	  pmanager->AddRestProcess(new G4KaonMinusAbsorptionBertini, ordDefault);
	}

      else if (particleName == "proton") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->
					GetCrossSectionDataSet(G4ChipsProtonElasticXS::Default_Name()));
          theElasticProcess->RegisterMe( elastic_chip );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  // Inelastic scattering
	  G4ProtonInelasticProcess* theInelasticProcess = 
	    new G4ProtonInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      else if (particleName == "anti_proton") 
	{
	  // Elastic scattering
          const G4double elastic_elimitAntiNuc = 100.0*CLHEP::MeV;
          G4AntiNuclElastic* elastic_anuc = new G4AntiNuclElastic();
          elastic_anuc->SetMinEnergy( elastic_elimitAntiNuc );
          G4CrossSectionElastic* elastic_anucxs = new G4CrossSectionElastic( elastic_anuc->GetComponentCrossSection() );
          G4HadronElastic* elastic_lhep2 = new G4HadronElastic();
          elastic_lhep2->SetMaxEnergy( elastic_elimitAntiNuc );
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( elastic_anucxs );
          theElasticProcess->RegisterMe( elastic_lhep2 );
          theElasticProcess->RegisterMe( elastic_anuc );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  // Inelastic scattering
	  G4AntiProtonInelasticProcess* theInelasticProcess = 
	    new G4AntiProtonInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( theAntiNucleonData );
	  theInelasticProcess->RegisterMe( theFTFModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	  // Absorption
	  pmanager->AddRestProcess(new G4AntiProtonAbsorptionFritiof, ordDefault);
	}

      else if (particleName == "neutron") {
	// elastic scattering
	G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
        theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsNeutronElasticXS::Default_Name()));
        G4HadronElastic* elastic_neutronChipsModel = new G4ChipsElasticModel();
	elastic_neutronChipsModel->SetMinEnergy( 19.0*CLHEP::MeV );
        theElasticProcess->RegisterMe( elastic_neutronChipsModel );
	G4NeutronHPElastic * theElasticNeutronHP = new G4NeutronHPElastic;
        theElasticNeutronHP->SetMinEnergy( theHPMin );
        theElasticNeutronHP->SetMaxEnergy( theHPMax );
	theElasticProcess->RegisterMe( theElasticNeutronHP );
	theElasticProcess->AddDataSet( new G4NeutronHPElasticData );
	pmanager->AddDiscreteProcess( theElasticProcess );
	// inelastic scattering		
	G4NeutronInelasticProcess* theInelasticProcess =
	  new G4NeutronInelasticProcess("inelastic");
	theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Neutron::Neutron() ) );
	theInelasticProcess->RegisterMe( theFTFModel1 );
        theInelasticProcess->RegisterMe( theBERTModel1 );
	G4NeutronHPInelastic * theNeutronInelasticHPModel = new G4NeutronHPInelastic;
        theNeutronInelasticHPModel->SetMinEnergy( theHPMin );
        theNeutronInelasticHPModel->SetMaxEnergy( theHPMax );
	theInelasticProcess->RegisterMe( theNeutronInelasticHPModel );
	theInelasticProcess->AddDataSet( new G4NeutronHPInelasticData );
	pmanager->AddDiscreteProcess(theInelasticProcess);
	// capture
	G4HadronCaptureProcess* theCaptureProcess =
	  new G4HadronCaptureProcess;
	G4NeutronHPCapture * theLENeutronCaptureModel = new G4NeutronHPCapture;
	theLENeutronCaptureModel->SetMinEnergy(theHPMin);
	theLENeutronCaptureModel->SetMaxEnergy(theHPMax);
	theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
	theCaptureProcess->AddDataSet( new G4NeutronHPCaptureData);
	pmanager->AddDiscreteProcess(theCaptureProcess);

      }
      else if (particleName == "anti_neutron") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering (include annihilation on-fly)
	  G4AntiNeutronInelasticProcess* theInelasticProcess = 
	    new G4AntiNeutronInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( theAntiNucleonData );
	  theInelasticProcess->RegisterMe( theFTFModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
	}

      else if (particleName == "deuteron") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4DeuteronInelasticProcess* theInelasticProcess = 
	    new G4DeuteronInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      
      else if (particleName == "triton") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4TritonInelasticProcess* theInelasticProcess = 
	    new G4TritonInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      else if (particleName == "alpha") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4AlphaInelasticProcess* theInelasticProcess = 
	    new G4AlphaInelasticProcess("inelastic");	 
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
    }
}


// Decays ///////////////////////////////////////////////////////////////////
#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"

void DSPhysicsList::ConstructGeneral() {

  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();

      if (theDecayProcess->IsApplicable(*particle) && !particle->IsShortLived())
        {
          pmanager ->AddProcess(theDecayProcess);
          // set ordering for PostStepDoIt and AtRestDoIt
          pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
          pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
        }
    }

  // Declare radioactive decay to the GenericIon in the IonTable.
  const G4IonTable *theIonTable =
    G4ParticleTable::GetParticleTable()->GetIonTable();
  G4RadioactiveDecay *theRadioactiveDecay = new G4RadioactiveDecay();
  theRadioactiveDecay->SetICM(true);                //Internal Conversion
  theRadioactiveDecay->SetARM(true);               //Atomic Rearangement

  for (G4int i=0; i<theIonTable->Entries(); i++)
    {
      G4String particleName = theIonTable->GetParticle(i)->GetParticleName();
      G4String particleType = theIonTable->GetParticle(i)->GetParticleType();

      if (particleName == "GenericIon")
        {
          G4ProcessManager* pmanager =
            theIonTable->GetParticle(i)->GetProcessManager();
          pmanager->SetVerboseLevel(VerboseLevel);
          pmanager ->AddProcess(theRadioactiveDecay);
          pmanager ->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
          pmanager ->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
        }
    }
}



// Cuts /////////////////////////////////////////////////////////////////////
void DSPhysicsList::SetCuts() {
  
  if(isDepositCuts)  { 
    DSLog(routine) << "Default Physics Cut for the DSGenneratorEnergyDeposit: 10 cm" << endlog ;
    defaultCutValue     = 10.0*cm; //
  }



  cutForGamma	      = defaultCutValue;
  cutForElectron      = defaultCutValue; //1.0*nanometer;
  cutForPositron      = defaultCutValue;

  
  DSLog(routine) << " Default Physics Cut Set to " << defaultCutValue/cm << " cm " << endlog ;

  //special for low energy physics
  G4double lowlimit=250*eV;  
  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  
  SetCutValue(1*cm, "opticalphoton");
  SetCutValue(defaultCutValue, "neutron");
  SetCutValue(1*micrometer,"GenericIon");
  
  if(isDepositCuts) SetCutValue(10*cm,"GenericIon");

  //  G4double lowlimit=1*keV;  
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit,100.*GeV);

  // Set Range Cut to 1 mm for gammas, electrons and positrons
  SetCutsWithDefault();


  // If enabled, set the range cuts to 1 micron in the active LAr volume
  if ( isHPRangeCuts ){
    double HPRangeCut = 1*um; 
    if(isDepositCuts)  HPRangeCut = 10*cm; 

    G4Region* activeLArRegion = G4RegionStore::GetInstance()->GetRegion("LAr_Logic");
    G4ProductionCuts* cutsForLAr = new G4ProductionCuts;
    cutsForLAr->SetProductionCut( HPRangeCut, "gamma" );
    cutsForLAr->SetProductionCut( HPRangeCut, "e-" );
    cutsForLAr->SetProductionCut( HPRangeCut, "e+" );
    cutsForLAr->SetProductionCut( HPRangeCut, "proton" );
    cutsForLAr->SetProductionCut( HPRangeCut, "neutron" );
    activeLArRegion->SetProductionCuts( cutsForLAr );
  }

  if (verboseLevel>0) DumpCutValuesTable();
}


void DSPhysicsList::DumpProcessTable() {
  
  int processmap[10000];
  for(int i=0;i<10000;i++) {
    processmap[i] = 0;
  }
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    G4ProcessVector *temp =  pmanager->GetProcessList();
    for(int i=0;i<pmanager->GetProcessListLength ();i++) {
      int index = 1000*(*temp)[i]->GetProcessType () + (*temp)[i]->GetProcessSubType () ;
      if(!processmap[index]) {
        processmap[index] = 1;
	DSLog(routine) << " Process Map Table - Index: "
	     << index << " ; Process Name: "
             << (*temp)[i]->GetProcessName() 
	     << endlog; 
	
       }
    }
  }

}

/*
 * $Log: DSPhysicsList.cc,v $
 * Revision 1.8  2015/04/29 14:49:13  dfranco
 * Fixed a bug in the DSGeneratorEnergyDeposit
 *
 * Revision 1.7  2015/01/30 15:38:38  dfranco
 * add neutrons to the list of high precision cuts in LAr
 *
 * Revision 1.6  2015/01/21 10:21:06  pagnes
 * DSLight3 set as default optics, cout removed
 *
 * Revision 1.5  2014/10/13 18:43:57  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.4  2014/07/23 14:45:54  meyers
 * Change WLS time distribution to exponential (was fixed offset)
 *
 * Revision 1.3  2014/06/27 14:54:27  perassos
 * Implementation of a new scintillation model
 *
 * Revision 1.2  2014/05/21 16:15:52  dfranco
 * added internal conversions to the RDM generator
 *
 * Revision 1.1  2014/05/07 12:21:04  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.36  2014/01/21 10:50:49  perassos
 * Range cuts set to 1um in LAr and to 1 mm elsewhere
 *
 * Revision 1.35  2013/10/20 16:30:14  swesterd
 * updated the veto scintillator optical properties based on Aldos measurements
 *
 * Revision 1.34  2013/10/01 06:26:24  swesterd
 * added waveform averager
 *
 * Revision 1.33  2013/08/27 07:13:14  swesterd
 * add visible energy for the neutron veto
 *
 * Revision 1.32  2013/08/27 04:07:01  swesterd
 * some fine tuning of the boron scintillator kB and scint yield, and some modifications to the DSG2 geometry
 *
 * Revision 1.31  2013/08/06 12:53:05  dfranco
 * Physics cuts changed
 *
 * Revision 1.30  2013/08/05 03:13:52  swesterd
 * some fine tuning of bscint and veto parameters
 *
 * Revision 1.29  2013/07/25 09:55:54  dfranco
 * Added a second version of NEST, still in development, and actually not working. The default version is DSLight
 *
 * Revision 1.28  2013/07/16 08:52:22  dfranco
 * added high precision (HP) neutron physics list for low energies
 *
 * Revision 1.27  2013/06/16 05:28:25  swesterd
 * fixed the makefile so that things work on argus again
 *
 * Revision 1.26  2013/06/06 23:19:35  swesterd
 * added quenching to BScint, currently commented out in DSPhysicsList. Have not compared with physical measurements...but might be giving reasonable results? Hard to tell....
 *
 * Revision 1.25  2013/06/03 08:37:08  dfranco
 * added cvs logger
 *
 * Revision 1.24  2013/06/03 08:35:16  dfranco
 * added command to kill the optical processes
 *
 * Revision 1.23  2013/05/28 14:48:58  dfranco
 * DSOpBoundaryProcess removed and subsstituted with G4OpBoundaryProcess.cc
 *
 * Revision 1.22  2013/05/27 23:59:02  swesterd
 * added a (currently commented out) Lumirror sheath to the cryostat and introduced DSOpBoundaryProcess to try to figure out why the boundaries are being screwy, with some edits so that it can handle constant and vector properties with freaking out
 *
 * Revision 1.21  2013/05/27 11:28:20  dfranco
 * Segmentation fault bug fixed, for the TPC case only. Useless files removed. Fixed the RINDEX, which was not defined in the whole range, for the scintillator.
 *
 *
 */
