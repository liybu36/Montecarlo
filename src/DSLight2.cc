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
#include "DSLight2.hh"
#include "DSLogger.hh"
#include "DSEventHandler.hh"
#include "G4EventManager.hh"
#include "G4SteppingManager.hh"
#include "G4ParticleTable.hh"
#include "G4TrackVector.hh"
#include <iostream>
#include <fstream>
  

using namespace std ;



///////////////////////////////////////////////////////////////////
DSLight2::DSLight2(const G4String& processName,G4ProcessType type):G4VRestDiscreteProcess(processName, type) {

  SetProcessType(fUserDefined);  
  fTrackSecondariesFirst = false ;
  
  fArgonEnergy.clear();
  fArgonELoss.clear();
  fElectronEnergy.clear();
  fElectronELoss.clear();
  fAlphaEnergy.clear();
  fAlphaELoss.clear();  

  ReadData();
  DSLog(debugging) << GetProcessName() << " is created "<< endlog ;

  tauFast = 6*ns;
  tauSlow = 1600*ns;
}

DSLight2::~DSLight2(){} 

G4VParticleChange* DSLight2::AtRestDoIt(const G4Track& aTrack, const G4Step& aStep){
  return DSLight2::PostStepDoIt(aTrack, aStep);
}

///////////////////////////////////////////////////////////////////
G4VParticleChange* DSLight2::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {

  // begin of the event
  if(aTrack.GetTrackID() == 1 && aTrack.GetCurrentStepNumber () == 1)  { 
    fEndOfEvent = false ;
    fEne = 0;
    DSStorage::Get()->SetEndOfEvent(0);
    DSStorage::Get()->SetDSLightTrackSecondaries(0);
    ClearDeposits() ;
  }
  
  const G4ParticleDefinition *particle = aTrack.GetParticleDefinition();
 
  trackVector =  G4EventManager::GetEventManager()->GetTrackingManager ()->GimmeSecondaries () ;
  int NTracks      = G4EventManager::GetEventManager()->GetStackManager()->GetNTotalTrack ();
  int TrackSize    = int(trackVector->size());
  int TrackStatus  = aTrack.GetTrackStatus ();
  int NStep        = aTrack.GetCurrentStepNumber ();
  



  // Look for the end of the event
  if(  DSStorage::Get()->GetEndOfEvent() == 0 
    && NTracks == 0  
    && TrackSize == 0 
    && TrackStatus > 0  
    && NStep > 1) {

    DSStorage::Get()->SetEndOfEvent(1) ;  
  } 

  fEne += aStep.GetTotalEnergyDeposit ();
  DSEventHandler::Get()->SetUserDouble(fEne/MeV);
  
    
  aParticleChange.Initialize(aTrack);
  const G4Material* aMaterial = aStep.GetPreStepPoint()->GetMaterial();  
  
  G4double dEdx = 0.;    // Units: MeV/mm

  
  if(aStep.GetTotalEnergyDeposit () > 0) { 
  //if(DSStorage::Get()->GetEndOfEvent() == 0) {
  
    // if not in LAr, return unchanged track
    if(DSStorage::Get()->GetLiquidArgonIndex() != (G4int) aMaterial->GetIndex()) 
                     return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
    
    // 
    if(fabs(particle->GetPDGEncoding()) == 11 || particle->GetPDGEncoding() == 22)           // electron/positron/ gamma
      dEdx = Interpolator(aTrack.GetKineticEnergy () + aStep.GetTotalEnergyDeposit (), fElectronEnergy, fElectronELoss);
    else if(particle->GetPDGEncoding() == 1000020040)    // alpha
      dEdx = Interpolator(aTrack.GetKineticEnergy () + aStep.GetTotalEnergyDeposit (), fAlphaEnergy, fAlphaELoss);
    else                                                 // argon or others
      dEdx = Interpolator(aTrack.GetKineticEnergy () + aStep.GetTotalEnergyDeposit (), fArgonEnergy, fArgonELoss);
  
    // save info about the energy deposition, and stack them into a vector
    fLocalDep.pdg        = particle->GetPDGEncoding() ;
    fLocalDep.t0         = aTrack.GetGlobalTime() ;
    fLocalDep.X0         = aStep.GetPreStepPoint()->GetPosition();
    fLocalDep.X1         = aStep.GetPostStepPoint()->GetPosition();
    fLocalDep.depEne     = aStep.GetTotalEnergyDeposit ();
    fLocalDep.kEne       = aTrack.GetKineticEnergy () + aStep.GetTotalEnergyDeposit ();
    fLocalDep.dEdx       = dEdx;
    fLocalDep.z          = particle->GetAtomicNumber() ;
    fLocalDep.numCluster = -1000;

    fDeposits.push_back(fLocalDep);

  }
  if(aStep.GetTotalEnergyDeposit () == 0 && DSStorage::Get()->GetEndOfEvent() == 0 )  
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  


  //const G4ElementVector* theElementVector1 = aMaterial->GetElementVector();
  //G4Element *ElementA = (*theElementVector1)[0];
  //G4double a1 = ElementA->GetA();

  if(DSStorage::Get()->GetEndOfEvent() == 1 && DSStorage::Get()->GetDSLightTrackSecondaries() == 0) {
    
  
    DSStorage::Get()->SetDSLightTrackSecondaries(1);
    DSStorage::Get()->SetEndOfEvent(0);
    
    // if no depositions, return unchanged track
    if(int(fDeposits.size()) == 0) {
      ClearDeposits() ;
      return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }

    G4double ScintillationYield = 4. / (19.5*eV); // scale factor for increasing the yield
    G4double ExcitationRatio = 0.21; // Doke 1992
    //G4double FanoFactor = 0.20; //Doke 1984
    G4double R0 =  1.568*um; //Mozumder 1995
    //G4double R0 =  0*um;  // Only Doke Birks
    //G4double R0 =  5*m; // Only Thomas Imel
    G4double DriftField = DSStorage::Get()->GetDriftField() ;
    

    G4double recombProb        = 0 ;
 
    G4double NumQuanta         = 0;
    G4double NumExcitons       = 0;
    G4double NumIons           = 0;
    G4double trackLength       = 0 ;
    G4int    NumElectrons      = 0;
    G4int    NumPhotons        = 0;
    
    G4double TotEne = 0 ;

    G4double       dep_dist = 0;
    G4double       dist_max = 1*mm;
    G4bool         isToBeAdded = true;
    G4ThreeVector  stepMeanPoint;
    G4double       stepLength = 0 ;

    //cout << "N deposits:  " << fDeposits.size() << endl;

    G4ParticleTable* myG4Table = G4ParticleTable::GetParticleTable();


    for( int i = 0; i < int( fDeposits.size() ); i++ ){

      isToBeAdded = true;

      if( myG4Table->FindParticle( fDeposits[i].pdg )->GetPDGCharge() ) 
        stepLength = (fDeposits[i].X1-fDeposits[i].X0).mag() ;
      else
        stepLength = fDeposits[i].depEne/fDeposits[i].dEdx ;  // mm




      // Clustering
      for( int j = 0; j < int( fDepClusters.size() ); j++) {

        dep_dist = sqrt( pow( fDepClusters[j].X1.x() - fDeposits[i].X1.x(), 2 ) + 
                         pow( fDepClusters[j].X1.y() - fDeposits[i].X1.y(), 2 ) + 
                         pow( fDepClusters[j].X1.z() - fDeposits[i].X1.z(), 2 ) ) ;

        if( dep_dist < dist_max ) {

          fDeposits[i].numCluster = j ;

          stepMeanPoint = ( fDeposits[i].X0 + fDeposits[i].X1 )/2. ;
          fDepClusters[j].Baricenter = fDepClusters[j].Energy * fDepClusters[j].Baricenter 
                                     + fDeposits[i].depEne    * stepMeanPoint ; 


          if( fDeposits[i].t0 < fDepClusters[j].T0 )  fDepClusters[j].T0 = fDeposits[i].t0 ; 
          fDepClusters[j].Length     += stepLength ;
          fDepClusters[j].Energy     += fDeposits[i].depEne ;
          fDepClusters[j].Baricenter /= fDepClusters[j].Energy ;
          fDepClusters[j].nDeposits++;


          isToBeAdded = false;
          break;
        }  

      }       


      // Otherwise, new cluster
      if( isToBeAdded ){
        fDeposits[i].numCluster = int( fDepClusters.size() ) ;

        //fLocalClus.Length     =  fDeposits[i].depEne/fDeposits[i].dEdx ;
        fLocalClus.Length     =  stepLength ; 
        fLocalClus.genPartPDG =  fDeposits[i].pdg ;
        fLocalClus.genPartZ   =  fDeposits[i].z ;
        fLocalClus.Energy     =  fDeposits[i].depEne ;
        fLocalClus.kinEne     =  fDeposits[i].kEne ;
        fLocalClus.dEdx       =  fDeposits[i].dEdx ;
        fLocalClus.T0         =  fDeposits[i].t0 ;
        fLocalClus.X0         =  fDeposits[i].X0 ;
        fLocalClus.X1         =  fDeposits[i].X1 ;
        fLocalClus.Baricenter =  (fDeposits[i].X1 + fDeposits[i].X0)/2. ;
        fLocalClus.nDeposits  =  1;


        fDepClusters.push_back( fLocalClus ) ;
      }

    } 




    // Compute the cluster's radius
    G4double  myRadius = 0.;
    //cout << "\n---------------------\n Number of Clusters:  " << fDepClusters.size() << "     Number of Deposits:  " << fDeposits.size() << endl;
    for( int i = 0; i < int( fDepClusters.size() ); i++ ){
   
      for( int j = 0; j < int( fDeposits.size() ); j++ ){
 
        if( fDeposits[j].numCluster != i ) continue ;
        
        myRadius = ( fDeposits[j].X1 - fDepClusters[i].Baricenter ).mag() ; 
        if( myRadius > fDepClusters[i].Radius )
          fDepClusters[i].Radius = myRadius ; 
        
      }

      //DSEventHandler::Get()->SetUserInt1( DSEventHandler::Get()->GetUserInt1() + fDepClusters[i].nDeposits );
      //DSEventHandler::Get()->SetUsers();

    }


    
    
    int TI = 0 ;
    int DB = 0 ;

    // loop over the en-dep clusters
    for( int i = 0; i < int( fDepClusters.size() ); ++i) {

      G4double Leff = 1 ;

      if( fabs( fDepClusters[i].genPartPDG ) != 11 && fDepClusters[i].genPartPDG != 22 ){
        G4double a1 = 39.948 ;
        G4double epsilon = 11.5 * fDepClusters[i].Energy/keV * pow( fDepClusters[i].genPartZ, (-7./3.) ) ;
        G4double gamma   = 3. * pow( epsilon, 0.15 ) + 0.7 * pow( epsilon, 0.6 ) + epsilon ;
        G4double kappa   = 0.133 * pow( fDepClusters[i].genPartZ, 2./3. ) * pow( a1, (-1./2.) ) * (2./3.) ;
        if( fDepClusters[i].genPartZ == 18 ) Leff = 0.23 * ( 1 + exp( -5 * epsilon ) ) ;
        else                                 Leff = (kappa*gamma)/(1+kappa*gamma) ;
      }
      

      G4double Quantas   = fDepClusters[i].Energy*ScintillationYield*Leff ;
      G4double Excitons  = Quantas*ExcitationRatio/(1+ExcitationRatio) ;
      G4double Ions      = Quantas - Excitons ;
      G4double Electrons = 0 ;


      if( fDepClusters[i].Radius < R0 ){ 
        TI++ ;
        //Thomas Imel
        G4double ThomasImel = 0 ;
        if(DriftField == 0)  ThomasImel = DSStorage::Get()->GetThomasImelNullField() ;
        else                 ThomasImel = DSStorage::Get()->GetThomasImelEp0() * 
                                            pow( DriftField/(volt/cm), DSStorage::Get()->GetThomasImelEp1() );

        G4double xi = Ions/4.*ThomasImel;
        recombProb = 1.0-log(1+xi)/xi;
        //cout << "Thomas-Imel" << endl;
	
      }
      else {
        DB++ ;
  
        // Doke Birks
        G4double DokeBirks1 = DSStorage::Get()->GetDokeBirksNFp1() ;
        G4double DokeBirks3 = DSStorage::Get()->GetDokeBirksNFp3() ;

        if( DriftField > 0) {
          DokeBirks1 = DSStorage::Get()->GetDokeBirksEp1()
                       *pow((DriftField/(kilovolt/cm)), DSStorage::Get()->GetDokeBirksEp2() );
          DokeBirks3 = DSStorage::Get()->GetDokeBirksEp3();
        }
        G4double DokeBirks2 = DokeBirks1/(1-DokeBirks3);

        recombProb = (DokeBirks1*fDepClusters[i].dEdx/1.4)
	                   / ( 1 + DokeBirks2*fDepClusters[i].dEdx/1.4 ) + DokeBirks3;
        //cout << "Doke-Birks" << endl;
      }

      Electrons = ( 1 - recombProb ) * Ions ;

      NumElectrons  += int( Electrons );
      NumPhotons    += int( Quantas - Electrons );
      
      fDepClusters[i].nElectrons = int( Electrons );
      fDepClusters[i].nPhotons   = int( Quantas - Electrons ); 
      fDepClusters[i].nExcitons  = int( Excitons );


    }



    DSEventHandler::Get()->SetUserFloat1( 100.*TI/(TI+DB) );
    //DSEventHandler::Get()->SetUserInt1( TI );
    //DSEventHandler::Get()->SetUserInt2( DB );
    DSEventHandler::Get()->SetUsers();

    DSEventHandler::Get()->SetS1Energy(
      DSEventHandler::Get()->GetS1Energy() + (float(NumPhotons) / ScintillationYield)/keV
    );

    DSEventHandler::Get()->SetS2Energy(
      DSEventHandler::Get()->GetS2Energy() + (float(NumElectrons) / ScintillationYield)/keV
      //DSEventHandler::Get()->GetS2Energy() + NumElectrons
    );



    G4int nDepClusters = int( fDepClusters.size() );

    // Kill S1 and S2 if required
    if(DSStorage::Get()->GetKillS1S2()) {
      NumPhotons = 0 ;
      nDepClusters = 0;
    }


    aParticleChange.SetNumberOfSecondaries( NumPhotons + 1000*NumElectrons );


    // Generate S1 and S2
    for( int i =0; i < nDepClusters; ++i) {


      // Generate S1
      for( int j = 0; j < fDepClusters[i].nPhotons; j++) {

        // Momentum
        G4double cost = 1 - 2*G4UniformRand();
        G4double sint = sqrt( 1 - cost*cost );

        G4double phi  = 2*M_PI*G4UniformRand();
        G4double cosp = cos( phi );
        G4double sinp = sin( phi );

        G4double px = sint*cosp;
        G4double py = sint*sinp;
        G4double pz = cost;

        G4ThreeVector photonMomentum ( px, py, pz ); 

        
        // Polarization
        G4double sx = cost*cosp;
        G4double sy = cost*sinp;
        G4double sz = -sint;

        G4ThreeVector photonPolarization ( sx, sy, sz );
        G4ThreeVector perpendicular = photonMomentum.cross( photonPolarization );

        phi  = 2*M_PI*G4UniformRand() ;
        cosp = cos( phi ) ;
        sinp = sin( phi ) ;

        photonPolarization = cosp * photonPolarization + sinp * photonMomentum ;
        photonPolarization = photonPolarization.unit();


        // Energy
        G4double photonMean  = 9.69*eV; // 9.81*eV;    lambda = 126.8 nm
        G4double photonWidth = 0.22*eV; // 0.60*eV;    FWHM = 7.8 nm
        G4double sampledEnergy = G4RandGauss::shoot( photonMean, photonWidth );


        // Position
        G4ThreeVector photonPosition = fDepClusters[i].X0 + G4UniformRand()*(fDepClusters[i].X1 - fDepClusters[i].X0) ;
        

        // Time
        //         photonTime = Global Time ( + Recombination Time ) + Singlet/Triplet Time

        G4double tauR                    = 0.8*ns;                // Recombination time

        G4double singletTripletRatioEx   = 0.;                    // Singlet to triplet ratio for excitons
        G4double singletTripletRatioReco = 0.;                    // Singlet to triplet ratio for recombining electrons
                                                                  // The two ratios are computed for ER, NR or alphas separately

        // -- Electron Recoil
        if( fDepClusters[i].genPartZ != 18 && fDepClusters[i].genPartZ != 2 && fDepClusters[i].genPartPDG != 2112 ){
          singletTripletRatioEx   = 0.36;          // Kubota 1979
          singletTripletRatioReco = 0.5;           // Kubota 1979
    
          G4double let = fDepClusters[i].dEdx/1.4;
          if( let > 3 ){
            singletTripletRatioReco = 0.2701+0.003379*let-4.7338e-5*pow(let,2.)+8.1449e-6*pow(let,3.); // check
            singletTripletRatioEx   = singletTripletRatioReco;
          }
        }
        // -- Alphas
        else if( fDepClusters[i].genPartZ == 2 ){
          singletTripletRatioReco = (-0.065492+1.9996*exp(-fDepClusters[i].Energy/MeV))/(1+0.082154/pow(fDepClusters[i].Energy/MeV,2.)) + 2.1811;   // check
          singletTripletRatioEx   = singletTripletRatioReco;
        }
        // -- Nuclear Recoils
        else{   
          singletTripletRatioReco = 0.22218*pow(fDepClusters[i].kinEne/keV,0.48211);  // check
          singletTripletRatioEx   = singletTripletRatioReco;
        }



        G4double singletFractionEx   = singletTripletRatioEx  /(1 + singletTripletRatioEx) ;
        G4double singletFractionReco = singletTripletRatioReco/(1 + singletTripletRatioReco) ;
        
        
        G4double photonTime = fDepClusters[i].T0;
        // -- Excitons
        if( j < fDepClusters[i].nExcitons ){ 

          if( G4UniformRand() < singletFractionEx ) 
            photonTime -= tauFast * log( G4UniformRand() ) ;
          else 
            photonTime -= tauSlow * log( G4UniformRand() ) ;
      
        }
        // -- Recombination electrons
        else{  

          photonTime -= tauR * log( G4UniformRand() ) ;   // It should be non-exponential

          if( G4UniformRand() < singletFractionReco )  
            photonTime -= tauFast * log( G4UniformRand() );
          else     
            photonTime -= tauSlow * log( G4UniformRand() );

        }
     
        

        // Create a new Photon
        G4DynamicParticle* scintillationPhoton = new G4DynamicParticle( G4OpticalPhoton::OpticalPhoton(), photonMomentum ); 
        scintillationPhoton->SetPolarization( photonPolarization.x(), photonPolarization.y(), photonPolarization.z() );
        scintillationPhoton->SetKineticEnergy( sampledEnergy );
        
        // Create a new track
        G4Track* aSecondaryTrack = new G4Track( scintillationPhoton, photonTime, photonPosition );
        aParticleChange.AddSecondary( aSecondaryTrack );      

      }
      
 

      // Generate S2
      G4int    nnn = 0;
      G4double LArGArBoundaryZ = DSStorage::Get()->GetLArGArBoundaryPosZ()/cm ;

      if( DriftField > 0 && LArGArBoundaryZ != -100*m && !DSStorage::Get()->GetKillS2() ){

        for( int j = 0; j < fDepClusters[i].nElectrons; j++) {

          // Position

          G4double D_T = 4.8*cm2/s; //ICARUS NIM A527 (2004) 329
          G4double D_L = 18*cm2/s; //atrazhev  Timoshkim theory at 1kV/cm
          if ( DriftField < 100 && aMaterial->GetState() == kStateLiquid )   D_L = 8*cm2/s;

          G4double vDrift = GetLArDriftVelocity( aMaterial->GetTemperature(), DriftField)*cm/s;

          G4double enDepZ = fDepClusters[i].X0.z() ;

          G4double sigmaDT = sqrt( 2 * D_T * fabs( LArGArBoundaryZ - enDepZ ) / vDrift );//cm
          G4double sigmaDL = sqrt( 2 * D_L * fabs( LArGArBoundaryZ - enDepZ ) / vDrift );//cm
          G4double dr = fabs( G4RandGauss::shoot(0., sigmaDT) );


          G4ThreeVector aSecondaryPosition = fDepClusters[i].X0 + G4UniformRand()*(fDepClusters[i].X1 - fDepClusters[i].X0) ;

          G4double phi = 2 * M_PI * G4UniformRand();
          aSecondaryPosition[0] += cos(phi) * dr;
          aSecondaryPosition[1] += sin(phi) * dr;
          aSecondaryPosition[2]  = LArGArBoundaryZ + 1*um;
          //radius = sqrt(pow(aSecondaryPosition[0],2.)+pow(aSecondaryPosition[1],2.));


          // Time

          G4double tDrift         = fabs( LArGArBoundaryZ - enDepZ ) / vDrift;
          G4double dt             = G4RandGauss::shoot(0.,(sigmaDL/vDrift));
          G4double aSecondaryTime = fDepClusters[i].T0 + tDrift + dt;

          //kill thermal electrons according to purity...
          G4double TaueLAr   = 15.8*ms;                              //electron lifetime due to impurities //http://arxiv.org/pdf/0910.5087.pdf 
          G4double prob      = exp(-(tDrift/ns) / (TaueLAr/ns)); 
          G4double testvalue = G4UniformRand();
          if(testvalue>prob)
            continue;
          else{
            CreateS2( aSecondaryPosition, aSecondaryTime, &aParticleChange);
          }
          continue;

        }
      }

    }  
      


    fDepClusters.clear() ;
    ClearDeposits() ;
  }
  
  





  return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
}


//////////////////////////////////////////////
//  GetMeanFreePath
//////////////////////////////////////////////
G4double DSLight2::GetMeanFreePath(const G4Track&, G4double , G4ForceCondition* condition) {
  *condition = StronglyForced;
  
  return DBL_MAX;
}

G4double DSLight2::GetMeanLifeTime(const G4Track&, G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}

/////////////////////////////////////////////
//        BinomFluct 
/////////////////////////////////////////////
G4int DSLight2::BinomFluct ( G4int N0, G4double prob ) {
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
//  CalculateElectronLET (from NEST)
//////////////////////////////////////////////
G4double DSLight2::CalculateElectronLET ( G4double E) {
  G4double LET;  // in keV
  //LAr
  //use a spline fit to online ESTAR data
  //at energies <1 keV, use a different spline, determined manually by
  //generating sub-keV electrons in Geant4 and looking at their ranges, since
  //ESTAR does not go this low  

  if ( E >= 1 ) 
    LET = 116.70-162.97*log10(E)+99.361*pow(log10(E),2)-33.405*pow(log10(E),3)+6.5069*pow(log10(E),4)-0.69334*pow(log10(E),5)+.031563*pow(log10(E),6);
  else if ( E>0 && E<1 )
    LET = 100;
  else
    LET = 0;
  
  return LET;
}


//////////////////////////////////////////////
//  Read Energy Loss Data (form SRIM/TRIM)
//////////////////////////////////////////////
void DSLight2::ReadData() {
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

}
//////////////////////////////////////////////
//  Energy Loss Interpolator
//////////////////////////////////////////////
G4double DSLight2::Interpolator(double ene, vector<float>&v1, vector<float>&v2) {
  
  
  
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

G4bool DSLight2::IsApplicable(const G4ParticleDefinition& aParticleType) {
  // DF:
  // return always true
  // this is crucial to understand the origin and the end of the event
  // DO NOT MODIFY IT
  
  return true;
}


//////////////////////////////
//       Create S2
//////////////////////////////
void DSLight2::CreateS2 (G4ThreeVector position, G4double time, G4ParticleChange *apc){

  G4ThreeVector x1 = position;
  G4double      t1 = time;

  //set to gas
  const G4Material* aMaterial = DSMaterial::Get()->GetGaseousArgon();
  //G4double Density = aMaterial->GetDensity()/(g/cm3);
  G4double Pressure = aMaterial->GetPressure();
  //G4double nDensity = (Density/MolarMassLAr)*AVO;
  //G4int Phase = aMaterial->GetState();
  G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();
  
  //if ( !YieldFactor || !YieldFactorS2 ){
  //  return;
  //}

  //G4double ElectricField = aMaterialPropertiesTable->GetConstProperty("ELECTRICFIELD");
  G4double ExtractionField = DSStorage::Get()->GetExtractionField();
  
  // in Ar at 3kV/cm all extracted //http://arxiv.org/pdf/astro-ph/0701286.pdf
  //otherwise a probability has to be set
  
  G4int numberofsecondary;//from http://arxiv.org/pdf/1207.2292.pdf
  G4double phgasA = 0.0813; 
  G4double phgasB = 139;
  G4double phgasG = 30.6;
  G4double GASGAP = 0.25*cm;

  G4double phpercm  = phgasA * ExtractionField/(volt/cm) - phgasB*Pressure/bar - phgasG;
  //numberofsecondary = G4int(G4int(floor(YieldFactor))* G4RandGauss::shoot(phpercm*GASGAP/cm,sqrt(phpercm*GASGAP/cm)) + 0.5);
  numberofsecondary = G4int( G4RandGauss::shoot(phpercm*GASGAP/cm,sqrt(phpercm*GASGAP/cm)) + 0.5 );  // floor(YieldFactor) removed because always = 1
                                                                                             // to check 
  numberofsecondary = G4int(numberofsecondary/DSStorage::Get()->GetScaleS2() + 0.5) ;

  DSLog(debugging)<< "\n Exiting from DS2Light::DoIt -- "<< "NumberOfSecondaries = "<<" "<<numberofsecondary<<endlog;

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
    G4double aSecondaryTime = t1,
      SingTripRatio=.1; //guess: revisit!!!
    if(G4UniformRand()<SingTripRatio/(1+SingTripRatio))
      aSecondaryTime -= tauFast*log(G4UniformRand());
    else
      aSecondaryTime -= tauSlow*log(G4UniformRand());

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
G4double DSLight2::GetLArDriftVelocity(G4double tempinput, G4double efieldinput) {
  //from ATLAS internal note LARG-NO-058

  double p1=-0.016863;
  double p2 = -0.0083412;
  double p3 = 0.18088;
  double p4 = 8.9751;
  double p5 = 1.4614;
  double p6 = 0.32891;
  double T= tempinput;
  double T0 = 92.91;
  double E = efieldinput/1000.; //in kV/cm

  double vdrift = (p1*(T-T0)+1)*(p3*E*log(1+p4/E)+p5*pow(E,p6))+p2*(T-T0);

  //G4cout<<" test drift "<<vdrift<<G4endl;
  vdrift*=1e5;//cm/s
  return vdrift;

}




/*
 * $Log: DSLight2.cc,v $
 * Revision 1.1  2014/05/07 12:21:03  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.13  2014/03/31 17:38:17  dfranco
 * latest optics config + DSLigh2 scintillation scaling + DSMaterial fixes
 *
 * Revision 1.12  2014/01/29 13:13:41  perassos
 * Update of the electric field handling and of the Nuclear Recoils generator
 *
 * Revision 1.10  2013/08/16 15:33:44  perassos
 * Added S1 and S2 to DSLight2
 *
 * Revision 1.9  2013/08/07 14:13:20  dfranco
 * Fixed a bug
 *
 * Revision 1.8  2013/08/07 13:16:38  dfranco
 * Modified R0
 *
 * Revision 1.7  2013/08/07 13:03:12  dfranco
 * Added the Lindhard factor
 *
 * Revision 1.6  2013/08/06 13:58:20  dfranco
 * Added 3 variables to store the deposited energy in LAr, scintillator, and water. The last two are not yet implemented. g4rooter has been updated with 3 new variables: tpcene, vetoene, and muene
 *
 * Revision 1.5  2013/08/02 15:45:58  dfranco
 * Further development on DSLight2
 *
 * Revision 1.4  2013/08/02 12:25:23  dfranco
 * Development of new DSLight class
 *
 * Revision 1.3  2013/08/01 14:48:32  dfranco
 * *** empty log message ***
 *
 * Revision 1.2  2013/08/01 14:28:08  dfranco
 * added energy loss data from SRIM/TRIM
 *
 * Revision 1.1  2013/07/25 09:55:54  dfranco
 * Added a second version of NEST, still in development, and actually not working. The default version is DSLight
 *
 *
 */
