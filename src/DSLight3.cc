////////////////////////////////////////////////////////////////////////
//
// Author: D. Franco (dfranco@in2p3.fr)
//
// S1 S2  Scintillation Light Class Implementation
////////////////////////////////////////////////////////////////////////

#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4EmProcessSubType.hh" 
#include "DSMaterial.hh"
#include "DSStorage.hh"
#include "DSLight3.hh"
#include "DSLogger.hh"
#include "DSEventHandler.hh"
#include "G4EventManager.hh"
#include "G4SteppingManager.hh"
#include "G4ParticleTable.hh"
#include "G4TrackVector.hh"
#include "DSParameters.hh"
#include "DSParameters.hh"
#include "DSEventHandler.hh"
#include "G4EventManager.hh"

using namespace std ;

/* 
* some important parameters: 

* LAr scintillation properties are read from a specific file (see ReadData() ) 
* The drift time calculation is based on the fLArGArBoundaryZ value (set in DSDetectorConstruction.cc ) 
*/


///////////////////////////////////////////////////////////////////
DSLight3::DSLight3(const G4String& processName,G4ProcessType type):G4VRestDiscreteProcess(processName, type) {

  SetProcessType(fUserDefined);  

  fTrackSecondariesFirst = false; 

  //isFirstEnDepInLAr = true;

  fRecoProb = 0;  
  fLArGArBoundaryZ = DSStorage::Get()->GetLArGArBoundaryPosZ() ;
  fDriftField = DSStorage::Get()->GetDriftField() ;
  
  if(DSStorage::Get()->GetTunedS1At200V()) fDriftField = 200.0*volt/cm ;
  
  if (DSStorage::Get()->GetKillS1S2() ) {
    DSStorage::Get()->SetKillS1(true) ; 
    DSStorage::Get()->SetKillS2(true) ; 
  }

  fArgonEnergy.clear();
  fArgonELoss.clear();

  fElectronEnergy.clear();
  fElectronELoss.clear();

  fAlphaEnergy.clear();
  fAlphaELoss.clear();  

  fRecombinationEnergy.clear();
  fRecombinationProbability.clear();

  ReadData();

  DSLog(debugging) << GetProcessName() << " is created "<< endlog ;
}


DSLight3::~DSLight3(){} 


G4VParticleChange* DSLight3::AtRestDoIt(const G4Track& aTrack, const G4Step& aStep){
  return DSLight3::PostStepDoIt(aTrack, aStep);
}



///////////////////////////////////////////////////////////////////
G4VParticleChange* DSLight3::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {


  aParticleChange.Initialize(aTrack);
  //cout << G4EventManager::GetEventManager()->GetStackManager()->GetNTotalTrack () << endl ;

  // At the event begin
  //if(aTrack.GetTrackID() == 1 && aTrack.GetCurrentStepNumber () == 1)  { 
  //  isFirstEnDepInLAr = true;
  //  DSStorage::Get()->SetDSLightTrackSecondaries(1);
  //}

  const G4ParticleDefinition* myParticle = aTrack.GetParticleDefinition();    
  const G4Material* myMaterial = aStep.GetPreStepPoint()->GetMaterial();  

      // if optical photon...
  if( myParticle->GetPDGEncoding() == 50 ||
      // ... or if not in LAr...
      DSStorage::Get()->GetLiquidArgonIndex() != (G4int) myMaterial->GetIndex() ||
      // ... or if deposited energy is null...
      aStep.GetTotalEnergyDeposit () <= 0  ||
      // ... or if the charge is zero 
      ( myParticle->GetPDGCharge() == 0  && myParticle->GetPDGEncoding() != 22  ) )  
      // ... return unchanged track.
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);




  G4double      myPDG    = myParticle->GetPDGEncoding() ;
  G4double      myZ      = myParticle->GetAtomicNumber() ;
  G4double      myT0     = aStep.GetPreStepPoint()->GetGlobalTime() ;
  G4ThreeVector myX0     = aStep.GetPreStepPoint()->GetPosition();
  G4double      myT1     = aStep.GetPostStepPoint()->GetGlobalTime() ;
  G4ThreeVector myX1     = aStep.GetPostStepPoint()->GetPosition();
  G4double      myDepEne = aStep.GetTotalEnergyDeposit ();
  //G4double      myKinEne = aTrack.GetKineticEnergy () + aStep.GetTotalEnergyDeposit ();
  G4double      myKinEne = aStep.GetPreStepPoint()->GetKineticEnergy();

  vector<float> myNrgVect   = fArgonEnergy;
  vector<float> myELossVect = fArgonELoss;
  if( fabs( myPDG ) == 11 || myPDG == 22) {  // electron/positron/ gamma
    myNrgVect   = fElectronEnergy;
    myELossVect = fElectronELoss;
  } else if( myPDG == 1000020040) {         // alpha
    myNrgVect   = fAlphaEnergy;
    myELossVect = fAlphaELoss;
  }
  G4double myDEDx  = Interpolator( myKinEne, myNrgVect, myELossVect );


  // G4double myFanoFactor = 0.2;
  G4double myExcitationRatio = myZ == 18 ? fExcitationRatioNR : fExcitationRatioER;
           



  //////////////////////////
  //                      //
  //   Quenching Factor   //
  //                      //
  //////////////////////////

  fQuenchingFactor = myZ == 18 ? 0.25 : 1;


  ///////////////////////////////////
  //                               //
  //   Recombination Probability   //
  //                               //
  ///////////////////////////////////

  G4double InitialKinEne = fQuenchingFactor * aTrack.GetVertexKineticEnergy();


  fRecoProb  = GetRecoProb(InitialKinEne,fDriftField);


  //DSEventHandler::Get()->SetUserDouble( fRecoProb );
  //DSEventHandler::Get()->SetUsers();

  G4double myNumQuanta   = fQuenchingFactor * myDepEne /fMeanQuantaEnergy ;
  G4double myNumIons     = myNumQuanta / ( 1 + myExcitationRatio ) ;
  G4double myNumExcitons = myNumQuanta - myNumIons;

  
  G4double myphotons   = myNumExcitons + myNumIons*fRecoProb ;
  G4double myelectrons = myNumIons*(1-fRecoProb) ;

  //compute the visible energy as an input to the model of f90
  G4double myVisEne = fQuenchingFactor * InitialKinEne *fLightYield;

  G4int myNumPhotons   = 0;
  G4int myNumElectrons = 0;
  
  if(myphotons < 20)     myNumPhotons = int(G4Poisson(myphotons) + 0.5);
  else                   myNumPhotons = int(G4RandGauss::shoot(myphotons, sqrt(myphotons))+0.5);
  if(myNumPhotons < 0)   myNumPhotons = 0 ;
  
  if(myelectrons < 20)   myNumElectrons = int(G4Poisson(myelectrons) + 0.5);
  else                   myNumElectrons = int(G4RandGauss::shoot(myelectrons, sqrt(myelectrons))+0.5);
  if(myNumElectrons < 0) myNumElectrons = 0 ;
  

  // G4int myNumPhotons   = int ( myNumExcitons + myNumIons*fRecoProb + 0.5);
  // G4int myNumElectrons = int ( myNumQuanta + 0.5 ) - myNumPhotons;        // + 0.5 for the rounding
  bool isInActiveLAr = aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName() == "ActiveLAr" ? 1 : 0;  
  if( isInActiveLAr ){
    DSEventHandler::Get()->SetS1Energy(
      DSEventHandler::Get()->GetS1Energy() + (float( myNumPhotons ) * fMeanQuantaEnergy )/keV
    );
    DSEventHandler::Get()->SetS2Energy(
      DSEventHandler::Get()->GetS2Energy() + (float( myNumElectrons )* fMeanQuantaEnergy  )/keV
    );
  }
  
  // scale by a factor --> LArScintillationProperties ; 
  myNumPhotons *= fS1ScaleFactor ;  
  // Kill S1 and S2 if required   (see the constructor - we want to save electrons!) 
  //if( DSStorage::Get()->GetKillS1S2() ) 
  //  return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);


  aParticleChange.SetNumberOfSecondaries( myNumPhotons + 1000 * myNumElectrons );
  if( fTrackSecondariesFirst ){
    if( aTrack.GetTrackStatus() == fAlive )
      aParticleChange.ProposeTrackStatus(fSuspend);
  }

  if (DSStorage::Get()->GetFastSimulation() )  myNumPhotons*=DSParameters::Get()->GetRealPmtMaxQe() ;
  
  // if KillS1 or KillS1S2, No S1 photons
  if (DSStorage::Get()->GetKillS1())  myNumPhotons = 0 ;

  // Generate S1
  for( int j = 0; j < (int) myNumPhotons; j++) {

    // Momentum
    G4double cost = 1 - 2*G4UniformRand();
    G4double sint = sqrt( 1 - cost*cost );

    G4double phi  = 2*M_PI*G4UniformRand();
    G4double cosp = cos( phi );
    G4double sinp = sin( phi );

    G4double px = sint*cosp;
    G4double py = sint*sinp;
    G4double pz = cost;

    G4ThreeVector myPhotonMomentum ( px, py, pz ); 

    // Polarization
    G4double sx = cost*cosp;
    G4double sy = cost*sinp;
    G4double sz = -sint;

    G4ThreeVector myPhotonPolarization ( sx, sy, sz );
    G4ThreeVector myPerpendicular = myPhotonMomentum.cross( myPhotonPolarization );

    phi  = 2*M_PI*G4UniformRand() ;
    cosp = cos( phi ) ;
    sinp = sin( phi ) ;

    myPhotonPolarization = cosp * myPhotonPolarization + sinp * myPhotonMomentum ;
    myPhotonPolarization = myPhotonPolarization.unit();

    // Energy
    G4double mySampledEnergy = 0;
    while( mySampledEnergy <= 0 ) 
      mySampledEnergy = G4RandGauss::shoot( fPhotonEneMean, fPhotonEneWidth );

    // Position
    G4ThreeVector myPhotonPosition = myX0 + G4UniformRand()*( myX1 - myX0 );
    
    // Time
    //   photonTime = Global Time ( + Recombination Time ) + Singlet/Triplet Time
    G4double myPhotonTime = myT0 + G4UniformRand()*(myT1 - myT0);

    //G4double mySingletTripletRatioEx   = 0.;        // Singlet to triplet ratio for excitons
    //G4double mySingletTripletRatioReco = 0.;        // Singlet to triplet ratio for recombining electrons
    G4double mySingletTripletRatio = 0.;
    G4double ratio_p1 = 0.;
    G4double ratio_p2 = 0.;
    G4double ratio_p3 = 0.;
    G4double ratio_p4 = 0.;



    // -- Electron Recoil
    if( myZ != 18 && myZ != 2 ){

      if( myVisEne < 200 ){
	ratio_p1 = 2.87719e-01;
	ratio_p2 = 1.54526e-01;
	ratio_p3 = 1.98113e-02;
	ratio_p4 = -7.05507e-05;
      }
      else if( myVisEne < 450 ){
	ratio_p1 = 2.98454e-01;
	ratio_p2 = 1.72632e+01;
	ratio_p3 = 4.63456e-02;
	ratio_p4 = -2.00378e-05;
      }
      else if( myVisEne < 1000 ){
	ratio_p1 = 2.99995e-01;
	ratio_p2 = -1.30386e+06;
	ratio_p3 = 2.35527e-01;
	ratio_p4 = -1.80290e-05;
      }
      else {
	ratio_p1 = 3.14439e-01;
	ratio_p2 = -1.35561e+06;
	ratio_p3 = 2.39579e-01;
	ratio_p4 = -1.19891e-06;
      }

      mySingletTripletRatio = (ratio_p1 - ratio_p4*myVisEne + ratio_p2*exp(-ratio_p3*myVisEne)) / 1.13 ; ///1.13 from the difference between the mean of f90 and the fraction of singlet and triplet 

    }
    // -- Alphas
    else if( myZ == 2 ){
      mySingletTripletRatio = (-0.065492+1.9996*exp(-myDepEne/MeV))/(1+0.082154/pow(myDepEne/MeV,2.)) + 2.1811;   // check
    }
    // -- Nuclear Recoils
    else{   
     ratio_p1 = -2.28001e-01 ;
     ratio_p2 = -8.14668e-03 ;
     ratio_p3 = 2.23409e-05 ;
     ratio_p4 = 8.03127e-01 ;
     
     if(myVisEne>2000) //constant f90 for >2000 Photons
       mySingletTripletRatio = (ratio_p1*exp(ratio_p2*2000) + ratio_p3*2000 + ratio_p4) * 1.12;
     else
       mySingletTripletRatio = (ratio_p1*exp(ratio_p2*myVisEne) + ratio_p3*myVisEne + ratio_p4) * 1.12;
     
    }
 
    G4UniformRand() < mySingletTripletRatio
      ? myPhotonTime -= fTauFast * log( G4UniformRand() ) 
      : myPhotonTime -= fTauSlow * log( G4UniformRand() ) ;

    /* Original model for f90 
    // -- Electron Recoil
    if( myZ != 18 && myZ != 2 && myPDG != 2112 ){
      mySingletTripletRatioEx   = 0.36;          // Kubota 1979
      mySingletTripletRatioReco = 0.5;           // Kubota 1979

      G4double myLET = myDEDx / 1.4;
      if( myLET > 3 ){
        mySingletTripletRatioReco = 0.2701+0.003379*myLET - 4.7338e-5*pow(myLET,2.)+8.1449e-6*pow(myLET,3.); //check
        mySingletTripletRatioEx   = mySingletTripletRatioReco;
      } 
    }
    // -- Alphas
    else if( myZ == 2 ){
      mySingletTripletRatioReco = (-0.065492+1.9996*exp(-myDepEne/MeV))/(1+0.082154/pow(myDepEne/MeV,2.)) + 2.1811;   // check
      mySingletTripletRatioEx   = mySingletTripletRatioReco;
    }
    // -- Nuclear Recoils
    else{   
      mySingletTripletRatioReco = 0.22218*pow( myKinEne/keV,0.48211 );  // check
      mySingletTripletRatioEx   = mySingletTripletRatioReco;
    }


    G4double mySingletFractionEx   = mySingletTripletRatioEx  /(1 + mySingletTripletRatioEx) ;
    G4double mySingletFractionReco = mySingletTripletRatioReco/(1 + mySingletTripletRatioReco) ;    
    


    if( j < myNumExcitons ){  // -- Excitons

      G4UniformRand() < mySingletFractionEx 
        ? myPhotonTime -= fTauFast * log( G4UniformRand() ) 
        : myPhotonTime -= fTauSlow * log( G4UniformRand() ) ;
  
    }
    else{  // -- Recombination Electrons 

      myPhotonTime -= fTauReco * log( G4UniformRand() ) ;   // It should be non-exponential

      G4UniformRand() < mySingletFractionReco 
        ? myPhotonTime -= fTauFast * log( G4UniformRand() ) 
        : myPhotonTime -= fTauSlow * log( G4UniformRand() ) ;

    }
    */

   
    

    // Create a new Photon
    G4DynamicParticle* myScintillationPhoton = new G4DynamicParticle( G4OpticalPhoton::OpticalPhoton(), myPhotonMomentum ); 
    myScintillationPhoton->SetPolarization( myPhotonPolarization.x(), myPhotonPolarization.y(), myPhotonPolarization.z() );
    myScintillationPhoton->SetKineticEnergy( mySampledEnergy );
    
    // Create a new track
    G4Track* aSecondaryTrack = new G4Track( myScintillationPhoton, myPhotonTime, myPhotonPosition );
    aParticleChange.AddSecondary( aSecondaryTrack );      

  }

  // Generate S2
  if( isInActiveLAr && fDriftField > 0 && fLArGArBoundaryZ != -100*m ){

   for( int j = 0; j < myNumElectrons; j++) {
      // Position
      G4ThreeVector aSecondaryPosition = myX0 + G4UniformRand()*(myX1 - myX0);
      
      
      G4double myEnDepZ = aSecondaryPosition.z();
      G4double myVDrift = GetLArDriftVelocity( myMaterial->GetTemperature(), fDriftField);
     
      G4double mySigmaDT = sqrt( 2 * fD_T * fabs( fLArGArBoundaryZ - myEnDepZ ) / myVDrift );//cm
      G4double dr = fabs( G4RandGauss::shoot(0., mySigmaDT) );
      G4double myPhi = 2 * M_PI * G4UniformRand();

      aSecondaryPosition[0] += cos(myPhi) * dr;
      aSecondaryPosition[1] += sin(myPhi) * dr;
      G4ThreeVector aSecondaryPosition_old = aSecondaryPosition ; 
      aSecondaryPosition[2]  = fLArGArBoundaryZ + 1*um;

      // Time
      if ( fDriftField < 100*volt/cm && myMaterial->GetState() == kStateLiquid ) fD_L = 8*cm2/s;
      
      G4double mySigmaDL      = sqrt( 2 * fD_L * fabs( fLArGArBoundaryZ - myEnDepZ ) / myVDrift );//cm
      G4double dt             = G4RandGauss::shoot(0., ( mySigmaDL/myVDrift ));
      G4double myTDrift       = ( fLArGArBoundaryZ - myEnDepZ ) / myVDrift;     // fabs shouldn't be there
      
      G4double aSecondaryTime = myT0 + G4UniformRand()*(myT1 - myT0) + myTDrift + dt;
      
      // the aSecondaryPosition  spread may cause some electrons 
      // to be generated inside the gas pocket (for boundary evts).
      // we will these electrons
      if (aSecondaryTime < 0.0 )  continue ; 

      // Kill Thermal Electrons according to purity...
      G4double prob           = exp(-(myTDrift/ns) / (fTaueLAr/ns));
      G4double myTestVal      = G4UniformRand();
      if( myTestVal > prob )	continue;
      
      else  {
        if(DSStorage::Get()->GetKillS2()==false) {
          CreateS2( aSecondaryPosition, aSecondaryTime, &aParticleChange);
        }
        if(DSStorage::Get()->GetWriteThermalElectrons() ) {
          DSEventHandler::Get()->SetDepPID(-1);
          DSEventHandler::Get()->SetDepVolume(0);
          DSEventHandler::Get()->SetDepEnergy(0);
          DSEventHandler::Get()->SetDepStep(0);
          DSEventHandler::Get()->SetDepTime(aSecondaryTime/ns);
          DSEventHandler::Get()->SetDepPosition(aSecondaryPosition_old/cm);
          DSEventHandler::Get()->SetDeposits();
        }
      }
    }
  }
  // ** --------------------


  return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
}




//////////////////////////////////////////////
//  GetMeanFreePath
//////////////////////////////////////////////
G4double DSLight3::GetMeanFreePath(const G4Track&, G4double , G4ForceCondition* condition) {
  *condition = StronglyForced;
  
  return DBL_MAX;
}


G4double DSLight3::GetMeanLifeTime(const G4Track&, G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}



/////////////////////////////////////////////
//        BinomFluct 
/////////////////////////////////////////////
G4int DSLight3::BinomFluct ( G4int N0, G4double prob ) {
  G4int N1 = 0;
  if ( prob == 0.00 ) return N1;
  if ( prob == 1.00 ) return N0;
  
  if ( N0 < 10 ) {
    for(G4int i = 0; i < N0; i++) {
      if(G4UniformRand() < prob) N1++;
    }
  }
  else {
    G4double mean = N0*prob;
    G4double sigma = sqrt(N0*prob*(1-prob));
    N1 = G4int(G4RandGauss::shoot(mean,sigma)+0.5);
    if ( N1 > N0 ) N1 = N0;
    if ( N1 < 0 ) N1 = 0;
  }
  return N1;
}





//////////////////////////////////////////////
//  Read Energy Loss Data (form SRIM/TRIM)
//////////////////////////////////////////////
void DSLight3::ReadData() {
  float ene, val ;
  ifstream ff;
  
  ff.open("../data/physics/electron_eloss.dat");
  while(!ff.eof()) {
    ff >> ene >> val ;
    if(ff.eof()) break ;
    fElectronEnergy.push_back(ene*keV);
    fElectronELoss.push_back(val*10);
  }
  ff.close();

  ff.open("../data/physics/argon_eloss.dat");
  while(!ff.eof()) {
    ff >> ene >> val ;
    if(ff.eof()) break ;
    fArgonEnergy.push_back(ene*keV);
    fArgonELoss.push_back(val*10);
  }
  ff.close();

  ff.open("../data/physics/alpha_eloss.dat");
  while(!ff.eof()) {
    ff >> ene >> val ;
    if(ff.eof()) break ;
    fAlphaEnergy.push_back(ene*keV);
    fAlphaELoss.push_back(val*keV/micrometer);
  }
  ff.close();

  ff.open("../data/physics/recombination_probability_S1.dat");
  while( ff >> ene >> val ){
    fRecombinationEnergy.push_back( ene );
    fRecombinationProbability.push_back( val );
  }
  ff.close();
  
  //READ LAr Scintillation Properties

  //Default values, Real ones are Read in ReadData() --> ../data/physics/LArScintillationProperties.txt
  fMeanQuantaEnergy =  (19.5*eV);
  fExcitationRatioER = 0.21;
  fExcitationRatioNR = 1.;
  fPhotonEneMean  = 9.69*eV; // 9.81*eV;    lambda = 126.8 nm
  fPhotonEneWidth = 0.22*eV; // 0.60*eV;    FWHM = 7.8 nm
  fS1ScaleFactor  = 3.5 ; 
  fD_T = 4.8*cm2/s;   // Transverse Diffusion from ICARUS NIM A527 (2004) 329 
  fD_L = 18*cm2/s;    // Longitudinal Diffusion from atrazhev  Timoshkim theory at 1kV/cm
  fTauFast = 6*ns;
  fTauSlow = 1600*ns;
  fTauReco = 0.8*ns;  // Kubota Phys Rev B 20 (1979) 
  fTaueLAr = 15.8*ms;
  fLightYield = 7000 ; //to be tuned!!!! ph/MeV
  // end of default values 
  
  ff.open("../data/physics/LArScintillatioProperties.txt");
  string mys ; 
  while (getline (ff, mys ) ) {
    if (!mys.find("MeanQuantaEnergy")) {
     fMeanQuantaEnergy= atof(mys.substr(mys.find("=")+1).c_str()) * eV;
     DSLog(routine) << " MeanQuantaEnergy: " << fMeanQuantaEnergy << endlog; 
    }
    if (!mys.find("S1ScaleFactor")) {
     fS1ScaleFactor= atof(mys.substr(mys.find("=")+1).c_str());
     DSLog(routine) << " S1ScaleFactor: " << fS1ScaleFactor << endlog; 
    }
    if (!mys.find("S2Yield")) {
     fS2Yield= atof(mys.substr(mys.find("=")+1).c_str());
     DSLog(routine) << " S2Yield: " << fS2Yield << endlog; 
    }
    if (!mys.find("ExcitationRatioER")) {
     fExcitationRatioER= atof(mys.substr(mys.find("=")+1).c_str());
     DSLog(routine) << " ExcitationRatioER: " << fExcitationRatioER << endlog; 
    }
    if (!mys.find("ExcitationRatioNR")) {
     fExcitationRatioNR= atof(mys.substr(mys.find("=")+1).c_str());
     DSLog(routine) << " ExcitationRatioNR: " << fExcitationRatioNR << endlog; 
    }
    if (!mys.find("PhotonEneMean")) {
     fPhotonEneMean= atof(mys.substr(mys.find("=")+1).c_str())  * eV;
     DSLog(routine) << " PhotonEneMean: " << fPhotonEneMean << endlog; 
    }
    if (!mys.find("PhotonEneWidth")) {
     fPhotonEneWidth= atof(mys.substr(mys.find("=")+1).c_str())  * eV;
     DSLog(routine) << " PhotonEneWidth: " << fPhotonEneWidth << endlog; 
    }
    if (!mys.find("D_T")) {
     fD_T= atof(mys.substr(mys.find("=")+1).c_str())  * cm2/s;
     DSLog(routine) << " D_T: " << fD_T << endlog; 
    }
    if (!mys.find("D_L")) {
     fD_L= atof(mys.substr(mys.find("=")+1).c_str())  * cm2/s;
     DSLog(routine) << " D_L: " << fD_L << endlog; 
    }
    if (!mys.find("TauFast")) {
     fTauFast= atof(mys.substr(mys.find("=")+1).c_str())  * ns;
     DSLog(routine) << " TauFast: " << fTauFast << endlog; 
    }
    if (!mys.find("TauSlow")) {
     fTauSlow= atof(mys.substr(mys.find("=")+1).c_str())  * ns;
     DSLog(routine) << " TauSlow: " << fTauSlow << endlog; 
    }
    if (!mys.find("TauReco")) {
     fTauReco= atof(mys.substr(mys.find("=")+1).c_str())  * ns;
     DSLog(routine) << " TauReco: " << fTauReco << endlog; 
    }
    if (!mys.find("TaueLAr")) {
     fTaueLAr= atof(mys.substr(mys.find("=")+1).c_str())  * ms;
     DSLog(routine) << " TaueLAr: " << fTaueLAr << endlog; 
    }
  }
  
}



//////////////////////////////////////////////
//  Energy Loss Interpolator
//////////////////////////////////////////////
G4double DSLight3::Interpolator(double ene, vector<float>&v1, vector<float>&v2) {
  
  
  if(ene < v1.front()) return v2.front() ;
  if(ene > v1.back())  return v2.back() ;
  
  
  int ndim = int(v1.size());
  
  for(int i=0;i<ndim;i++) {
    if(v1[i]> ene) {
      double _mm = (v2[i-1] - v2[i])/(v1[i-1]-v1[i]);
      double q = v2[i] - _mm*v1[i];
      return (_mm*ene + q)  ;
    }
  }
  return 0. ;  

}



G4bool DSLight3::IsApplicable(const G4ParticleDefinition& aParticleType) {
  // DF:
  // return always true
  // this is crucial to understand the origin and the end of the event
  // DO NOT MODIFY IT
  
  return true;
}


//////////////////////////////
//       Create S2
//////////////////////////////
void DSLight3::CreateS2 (G4ThreeVector position, G4double time, G4ParticleChange *apc){

  G4ThreeVector x1 = position;
  G4double      myT1 = time;

  //set to gas
  /*
  
  const G4Material* aMaterial = DSMaterial::Get()->GetGaseousArgon();
  G4double Pressure = aMaterial->GetPressure();
  G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();
  G4double ExtractionField = DSStorage::Get()->GetExtractionField();
  
  // in Ar at 3kV/cm all extracted //http://arxiv.org/pdf/astro-ph/0701286.pdf
  //otherwise a probability has to be set
  
  G4int numberofsecondary;//from http://arxiv.org/pdf/1207.2292.pdf
  G4double phgasA = 0.0813; 
  G4double phgasB = 139;
  G4double phgasG = 30.6;
  G4double GASGAP = 0.25*cm;
  G4double phpercm  = phgasA * ExtractionField/(volt/cm) - phgasB*Pressure/bar - phgasG;
  numberofsecondary = G4int( G4RandGauss::shoot(phpercm*GASGAP/cm,sqrt(phpercm*GASGAP/cm)) + 0.5 );  // floor(YieldFactor) removed because always = 1                                                                                             // to check 
  numberofsecondary = G4int(numberofsecondary/DSStorage::Get()->GetScaleS2() + 0.5) ;
  */
  
  //float GASGAP = 0.25*cm;
  float GASGAP = DSParameters::Get()->GetGasPocketThickness();
  float meannumberofphotons = fS2Yield  ;
  float s2scaling_ = DSStorage::Get()->GetScaleS2() ; 
  if (s2scaling_ < 1.0  ) {
    DSLog(fatal) << "Fatal: ScaleS2 must be > 1.0 " << endlog ;
    s2scaling_ = 1. ; 
  }
    
  G4int numberofsecondary =    ( G4RandGauss::shoot(meannumberofphotons, sqrt(meannumberofphotons))+0.5 ) /  s2scaling_; 
  if(numberofsecondary < 0) numberofsecondary = 0 ;
  
  if (DSStorage::Get()->GetFastSimulation() )  numberofsecondary*=DSParameters::Get()->GetRealPmtMaxQe() ;
  
  DSLog(debugging)<< "\n Exiting from DS2Light::DoIt -- "<< "NumberOfSecondaries = "<<" "<<numberofsecondary<<endlog;

  numberofsecondary = int(numberofsecondary) ; 
  
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
    G4double E_eVLAr = 9.7;
    sampledEnergy = G4RandGauss::shoot(E_eVLAr*eV,0.2*eV);
    aQuantum = new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(),photonMomentum);
    aQuantum->SetPolarization(photonPolarization.x(),
            photonPolarization.y(),
            photonPolarization.z());

    //assign energy to make particle real
    aQuantum->SetKineticEnergy(sampledEnergy);
    //G4cout << "sampledEnergy = " << sampledEnergy << G4endl;

    // electroluminesence emission time distribution
    G4double aSecondaryTime = myT1,
      SingTripRatio=.1; //guess: revisit!!!
    if(G4UniformRand()<SingTripRatio/(1+SingTripRatio))
      aSecondaryTime -= fTauFast*log(G4UniformRand());
    else
      aSecondaryTime -= fTauSlow*log(G4UniformRand());

    G4ThreeVector aSecondaryPosition = x1;
    aSecondaryPosition[2]+= (G4UniformRand()*GASGAP/mm) * mm;

    // GEANT4 business: stuff you need to make a new track
    G4Track* aSecondaryTrack =
      new G4Track(aQuantum,aSecondaryTime,aSecondaryPosition);
    apc->AddSecondary(aSecondaryTrack);

  }
}


//////////////////////////////
//  GetLArDriftVelocity
//////////////////////////////
G4double DSLight3::GetLArDriftVelocity(G4double tempinput, G4double efieldinput) {
  //from ATLAS internal note LARG-NO-058

  double p1 = -0.016863;
  double p2 = -0.0083412;
  double p3 = 0.18088;
  double p4 = 8.9751;
  double p5 = 1.4614;
  double p6 = 0.32891;
  double T  = tempinput;
  double T0 = 92.91;
  double E  = efieldinput/(kilovolt/cm); //in kV/cm

  double vdrift = (p1*(T-T0)+1)*(p3*E*log(1+p4/E)+p5*pow(E,p6))+p2*(T-T0);

  vdrift*=1e5;//cm/s
  return vdrift*cm/s;

}

//////////////////////////////
//  Recombination Probability
//////////////////////////////
G4double DSLight3::GetRecoProb(double InitialKinEne, double DriftField) {
  //  - Extracted from Ar39 spectrum in DS10 data
  //  - 9 parameters (more detailed modelization of the field dependence)
  //  - Improved agreement with DS10 data in the region of the recombination peak
  double myRecoProb = 0 ;
  
  if(DSStorage::Get()->GetTunedS1At200V() == false ) {
  
    double p0 =  0.581234;
    double p1 =  0.720611;
    double p2 = -0.122622;
    double p3 =  0.000117115;
    double p4 = -0.226248;
    double p5 =  0.000311352;
    double p6 =  0.286737;
    double p7 = -0.000170164;
    double p8 =  0.567157;

    myRecoProb = p0 * (1 - p1*exp( p2 * InitialKinEne/keV )) * 
                exp( (p3 * DriftField*cm/volt + p4) * pow( InitialKinEne/keV, (p5*DriftField*cm/volt  + p6) ) ) + 
                p7*DriftField*cm/volt  + p8 ;
  
  } else {

  
    double   p0=    5.79472e-01  ;
    double   p1=    8.77165e-01  ;
    double   p2=   -5.14651e-02  ;
    double   p3=   -2.48278e-01  ;
    double   p4=    3.73392e-01  ;
    double   p5=    5.64454e-01;
  
    
    myRecoProb = p0 * (1 - p1*exp( p2 * InitialKinEne/keV )) * 
              exp( (p3) * pow( InitialKinEne/keV, (p4) ) ) + 
              p5 ;
    
    myRecoProb *= 19.5/23.6 ;
  }

  return myRecoProb;
}

/*
 * $Log: DSLight3.cc,v $
 * Revision 1.20  2015/04/20 08:28:53  pagnes
 * add light production to gamma deposits too
 *
 * Revision 1.19  2015/02/17 09:41:30  dfranco
 * removed a useless variable
 *
 * Revision 1.18  2015/01/30 13:57:43  dfranco
 * improve f90 modeling for neutrons
 *
 * Revision 1.17  2015/01/20 09:19:08  dfranco
 * do not produce light for neutral particles in DSLight3 and new model for f90
 *
 * Revision 1.16  2014/12/22 14:40:43  dfranco
 * added the option to activate the recombination probability at 200 V/cm (/ds/physics/tunedS1); this option is by default true; selecting a specific drift field automatically switch off the tunedS1 option
 *
 * Revision 1.15  2014/11/13 15:44:07  dfranco
 * remove printout of EventManager
 *
 * Revision 1.14  2014/11/13 15:26:21  dfranco
 * bug fixed on the recombination probability
 *
 * Revision 1.13  2014/11/12 13:24:50  dfranco
 * comment removed
 *
 * Revision 1.12  2014/11/12 13:18:40  dfranco
 * added units to drift field
 *
 * Revision 1.11  2014/10/28 13:56:22  perassos
 * Generation of S2 and update of s1ene and s2ene are done only if the deposit is in the ActiveLAr volume
 *
 * Revision 1.10  2014/07/25 14:11:07  perassos
 * Introduced the improved recombination probability model from DS10 data
 *
 * Revision 1.9  2014/07/23 14:55:34  pagnes
 * S2 first tuning
 *
 * Revision 1.8  2014/07/18 13:54:49  dfranco
 * Added a new position generator in the Gas Pocket; calibrated the numnber of photons generated per electron in S2; added a new command to generate multiple particles in the same position
 *
 * Revision 1.7  2014/07/17 13:26:11  dfranco
 * Added the CVS automatic comment
 *
 *
 */
