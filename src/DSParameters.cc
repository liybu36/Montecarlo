// --------------------------------------------------------------------------//
/** 
 * AUTHOR: Davide Franco
 * 
 */
// --------------------------------------------------------------------------//


//---------------------------------------------------------------------------//

#include "DSParameters.hh"      //Present DS Class Headers 
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//---------------------------------------------------------------------------//


#include <fstream>
#include <sstream>
#include "DSLogger.hh"
#include "DSIO.hh"
#include "DSStorage.hh"


using namespace std ;

DSParameters* DSParameters::me = 0;

// singleton
DSParameters::DSParameters(){
    //fLumirrorSkipNEntries = 10 ; 
    ReadDetectorGeometry();
    ReadDetectorProperties();
    ReadVetoPMTGeometry();
    ReadOpticalProperties();
    ReadOpticsTuning();
}

DSParameters* DSParameters::Get() {
  if (!me)
    me = new  DSParameters();
  return me;
}


void DSParameters::ReadDetectorGeometry() {
  G4String mys;
  //std::ifstream  inDetectorGeometryFile("../data/DSGeometry.dat",std::ios::in);

  DSIO::Get()->GetStreamLogFile() << endl ;
  DSIO::Get()->GetStreamLogFile() << "########## Detector Geometry ###########" << endl ;


  while (getline (DSIO::Get()->GetStreamDSGeometry(), mys ) ) {
  
    DSIO::Get()->GetStreamLogFile() << mys << endl ;

    if (!mys.find("WorldSizeX")) {
      fWorldSizeX = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("WorldSizeY")) {
      fWorldSizeY = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("WorldSizeZ")) {
      fWorldSizeZ = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("DSTankHeight")) {
      fTankHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("DSTankRMax")) {
      fTankRmax = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("SSSRadius")) {
      fSSSRadius = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("CryostatShiftZ")) {
      fCryostatShiftZ = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TrunkShiftZ")) {
      fTrunkShiftZ = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TrunkDiameter")) {
      fTrunkDiameter = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TrunkThickness")) {
      fTrunkThickness = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TrunkTopBottomOffsetZ")) {
      fTrunkTopBottomOffsetZ = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TrunkTopBottomOffsetX")) {
      fTrunkTopBottomOffsetX = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TrunkBottomHeight")) {
      fTrunkBottomHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TrunkMiddleHeight")) {
      fTrunkMiddleHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TrunkDistToCryostatAxis")) {
      fTrunkDistToCryostatAxis = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TPCShiftZ")) {
      fTPCShiftZ = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TeflonSupportDiameter")) {
      fTeflonSupportDiameter = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TeflonSupportThickness")) {
      fTeflonSupportThickness = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TeflonSupportHeight")) {
      fTeflonSupportHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TeflonSupportHeight")) {
      fTeflonSupportHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("PMTAssemblyHeight")) {
      fPMTAssemblyHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TeflonCapDiameter")) {
      fTeflonCapDiameter = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TeflonCapHeight")) {
      fTeflonCapHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("LArBottomLayerThickness")) {
      fLArBottomLayerThickness = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("CathodeWindowDiameter")) {
      fCathodeWindowDiameter = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("CathodeWindowHeight")) {
      fCathodeWindowHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("ITOThickness")) {
      fITOThickness = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("ReflectorInnerDiameter")) {
      fReflectorInnerDiameter = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("ReflectorOuterDiameter")) {
      fReflectorOuterDiameter = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("ReflectorHeight")) {
      fReflectorHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("AboveGridInnerDiameter")) {
      fAboveGridInnerDiameter = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("AboveGridOuterDiameter")) {
      fAboveGridOuterDiameter = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("AboveGridHeight")) {
      fAboveGridHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("GasPocketThickness")) {
      fGasPocketThickness = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("TPBThickness")) {
      fTPBThickness = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("DivingBellHeight")) {
      fDivingBellHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("DivingBellTopHeight")) {
      fDivingBellTopHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("DivingBellOuterDiameter")) {
      fDivingBellOuterDiameter = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("DivingBellInnerDiameter")) {
      fDivingBellInnerDiameter = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("FieldRingsHeight")) {
      fFieldRingsHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("FieldRingsThickness")) {
      fFieldRingsThickness = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("PMTBodyDiameter")) {
      fPMTBodyDiameter = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("PMTBodyHeight")) {
      fPMTBodyHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("PMTHeadDiameter")) {
      fPMTHeadDiameter = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("PMTHeadHeight")) {
      fPMTHeadHeight = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("PMTWallThickness")) {
      fPMTWallThickness = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("PMTWindowThickness")) {
      fPMTWindowThickness = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("PMTSpacing")) {
      fPMTSpacing = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }
    if (!mys.find("PMTOffset")) {
      fPMTOffset = atof(mys.substr(mys.find("=")+1).c_str()) * mm;
    }


  }

  DSLog(trace) << "Detector Geometry loaded" << endlog;
}

void DSParameters::ReadVetoPMTGeometry(){
  G4double pmtnum, theta, phi; // theta = polar angle, phi = azimuthal angle
  DSIO::Get()->GetStreamLogFile() << endl;
  DSIO::Get()->GetStreamLogFile() << "########## Veto PMT Positions ###########" << endl ;
  fVPMTNum.clear() ; 
  fVPMTTheta.clear() ; 
  fVPMTPhi.clear() ; 
  
  int index=0;
  while(DSIO::Get()->GetStreamDSVPMTGeometry() >> pmtnum >> theta >> phi) {
    if(index >= 110) break;
    //DSLog(trace) << "ReadVetoPMTGeometry() " << pmtnum << endlog;
    DSIO::Get()->GetStreamLogFile() << index<< ", pmt " << pmtnum << "\ttheta = "<< theta <<"\tphi = "<< phi << endl ;
    fVPMTNum.push_back(pmtnum);
    fVPMTTheta.push_back(theta);
    fVPMTPhi.push_back(phi);
    index++;
  }
  DSIO::Get()->CloseStreamDSVPMTGeometry();
}


void DSParameters::ReadOpticsTuning() {
  G4String mys;
  DSIO::Get()->GetStreamLogFile() << endl ;
  DSIO::Get()->GetStreamLogFile() << "########## Detector Optics Tuning ###########" << endl ;

  while (getline (DSIO::Get()->GetStreamDSOptics(), mys ) ) {
    DSIO::Get()->GetStreamLogFile() << mys << endl  ;

    if (!mys.find("GridSteelRindScale")) {
      fGridSteelRindScale = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "Scale factor Grid Steel RINDEX: " << fGridSteelRindScale << endlog;
    }
    if (!mys.find("PhotocathodeUVRind")) {
      fPhotocathodeUVRind = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "Photocathode RINDEX UV light: " << fPhotocathodeUVRind << endlog;
    }
    if (!mys.find("PhotocathodeVisRind")) {
      fPhotocathodeVisRind = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "Photocathode RINDEX Visible light: " << fPhotocathodeVisRind << endlog;
    }
    if (!mys.find("FusedSilicaUVRind")) {
      fFusedSilicaUVRind = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "Fused Silica RINDEX UV light: " << fFusedSilicaUVRind << endlog;
    }
    if (!mys.find("FusedSilicaVisRind")) {
      fFusedSilicaVisRind = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "Fused Silica RINDEX visible light: " << fFusedSilicaVisRind << endlog;
    }
    if (!mys.find("TPBUVRind")) {
      fTPBUVRind = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "TPB RINDEX UV light: " << fTPBUVRind << endlog;
    }
    if (!mys.find("TPBVisRind")) {
      fTPBVisRind = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "TPB RINDEX Visible light: " << fTPBVisRind << endlog;
    }
    if (!mys.find("GaseousArgonUVAbs")) {
      fGaseousArgonUVAbs = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "Gaseous Ar Absorption Length UV light: " << fGaseousArgonUVAbs << endlog;
    }
    if (!mys.find("GaseousArgonVisAbs")) {
      fGaseousArgonVisAbs = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "Gaseous Ar Absorption Length Vis light: " << fGaseousArgonVisAbs << endlog;
    }
    if (!mys.find("LiquidArgonUVAbs")) {
      fLiquidArgonUVAbs = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "Liquid Ar Absorption Length UV light: " << fLiquidArgonUVAbs << endlog;
    }
    if (!mys.find("LiquidArgonVisAbs")) {
      fLiquidArgonVisAbs = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "Liquid Ar Absorption Length Visible light: " << fLiquidArgonVisAbs << endlog;
    }
    if (!mys.find("GridSteelUVAbs")) {
      fGridSteelUVAbs = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "GridSteel Absorption Length UV light: " << fGridSteelUVAbs << endlog;
    }
    if (!mys.find("GridSteelVisAbs")) {
      fGridSteelVisAbs = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "GridSteel Absorption Length Visible light: " << fGridSteelVisAbs << endlog;
    }
    if (!mys.find("FusedSilicaUVAbs")) {
      fFusedSilicaUVAbs = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "FusedSilica Absorption Length UV light: " << fFusedSilicaUVAbs << endlog;
    }
    if (!mys.find("FusedSilicaVisAbs")) {
      fFusedSilicaVisAbs = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "FusedSilica Absorption Length Visible light: " << fFusedSilicaVisAbs << endlog;
    }
    if (!mys.find("TPBUVAbs")) {
      fTPBUVAbs = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "TPB Absorption Length UV light: " << fTPBUVAbs << endlog;
    }
    if (!mys.find("TPBVisAbs")) {
      fTPBVisAbs = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "TPB Absorption Length Visible light: " << fTPBVisAbs << endlog;
    }
    if (!mys.find("WLSAbsorptionFactor")) {
      fWLSAbsorptionFactor = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "TPB W.L.Shifting-Length scaling factor: " << fWLSAbsorptionFactor << endlog;
    }
    if (!mys.find("WLSMeanNumberPhotons")) {
      fWLSMeanNumberPhotons = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "TPB mean number of photons: " << fWLSMeanNumberPhotons << endlog;
    }
    if (!mys.find("WLSTimeConstant_ns")) {
      fWLSTimeConstant_ns = atof(mys.substr(mys.find("=")+1).c_str()) ;  
       //DSLog(trace) << "TPB time constant (ns): " << fWLSTimeConstant_ns << endlog;
    }
    if (!mys.find("WLSEfficiency")) {
      fWLSEfficiency = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "TPB Efficiency: " << fWLSEfficiency << endlog;
    }
    if (!mys.find("GArRindexScale")) {
      fGArRindexScale = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "GAr Rayleigh scattering scaling: " << fGArRayleighScale << endlog;
    }
    if (!mys.find("LArRayleighScale")) {
      fLArRayleighScale = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "LAr Rayleigh scattering scaling: " << fLArRayleighScale << endlog;
    }
    if (!mys.find("FSilicaRaylVisLength")) {
      fFSilicaRaylVisLength = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "LAr Rayleigh scattering scaling: " << fLArRayleighScale << endlog;
    }
    if (!mys.find("FSilicaRaylUVLength")) {
      fFSilicaRaylUVLength = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "LAr Rayleigh scattering scaling: " << fLArRayleighScale << endlog;
    }
   if (!mys.find("WithITO")) {
      fWithITO = atoi(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "ITO placement (bool): " << fWithITO << endlog;
    }
    if (!mys.find("WithGasPocket")) {
      fWithGasPocket = atoi(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "GasPocket placement (bool): " << fWithGasPocket << endlog;
     }
    if (!mys.find("WithNewGridModel")) {
      fWithNewGridModel = atoi(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "GasPocket placement (bool): " << fWithNewGridModel << endlog;
    }
    if (!mys.find("LArGridUVRef")) {
      fLArGridUVRef = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "LAr - Grid Surf Reflectivity (1 - abs), UV light: " << fLArGridUVRef << endlog;
    }
    if (!mys.find("LArGridVisRef")) {
      fLArGridVisRef = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "LAr - Grid Surf Reflectivity (1 - abs), vis light: " << fLArGridVisRef << endlog;
    }
    // P. Meyers: eliminated GridReflection, added GridNormalTransparency
    if (!mys.find("GridNormalTransparency")) {
      fGridNormalTransparency = atof(mys.substr(mys.find("=")+1).c_str()) ;  
    }
    if (!mys.find("TeflonTPBUVRef")) {
      fTeflonTPBUVRef = atof(mys.substr(mys.find("=")+1).c_str()) ;  
    }
    if (!mys.find("TeflonTPBVisRef")) {
      fTeflonTPBVisRef = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "Teflon-TPB Surf Reflectivity (1 - abs), vis light: " << fTeflonTPBVisRef << endlog;
    }
    if (!mys.find("TeflonLArUVRef")) {
      fTeflonLArUVRef = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "Teflon-LAr Surf Reflectivity (1 - abs), UV light: " << fTeflonLArUVRef << endlog;
    }
    if (!mys.find("TeflonLArVisRef")) {
      fTeflonLArVisRef = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "Teflon-LAr Surf Reflectivity (1 - abs), vis light: " << fTeflonLArVisRef << endlog;
    }
    if (!mys.find("PMTLArUVRef")) {
      fPMTLArUVRef = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "PMT-LAr Surf Reflectivity (1 - abs), UV light: " << fPMTLArUVRef << endlog;
    }
    if (!mys.find("PMTLArVisRef")) {
      fPMTLArVisRef = atof(mys.substr(mys.find("=")+1).c_str()) ;  
      //DSLog(trace) << "PMT-LAr Surf Reflectivity (1 - abs), vis light: " << fPMTLArVisRef << endlog;
    }
    if (!mys.find("ArTPBVisTran")) {
      fArTPBVisTran = atof(mys.substr(mys.find("=")+1).c_str()) ;
      //DSLog(trace) << "PMT-LAr Surf Reflectivity (1 - abs), vis light: " << fPMTLArVisRef << endlog;
    }
    if (!mys.find("PArRind")) {
      fPArRind = atof(mys.substr(mys.find("=")+1).c_str()) ;
      //DSLog(trace) << "Index of refraction of the pseudo argon: " << PArRind << endlog;
    }
  }

  DSLog(trace) << "#######    Detector Optics loaded      #######" << endlog;
}







void DSParameters::ReadDetectorProperties() {
  
  // Load Detector Properties
  ifstream myPropertiesFile( "../data/detector/DSProperties.dat" );
  string myString;


  if( myPropertiesFile.is_open() ){
    while( getline( myPropertiesFile, myString ) ){

      if(!myString.find("PmtMaxQe")) {
        fPmtMaxQe = atof( myString.substr( myString.find("=")+1 ).c_str() );
	fRealPmtMaxQe = fPmtMaxQe ; 
      }
      //PA QE scaling
      if (DSStorage::Get()->GetFastSimulation() ) fPmtMaxQe = 1. ; 
      if(!myString.find("VPmtMaxQe")) {
	fVPmtMaxQe = atof( myString.substr( myString.find("=")+1 ).c_str() );
      }

    }
    myPropertiesFile.close();
    DSLog(trace) << "Detector Properties loaded" << endlog;  
  }
  else {
    DSLog(trace) << "Couldn't load Detector Properties" << endlog;  
  }


  
  // Load TPC PMT Q.E.
  ifstream myQeFile( "../data/detector/DSPmtQeDS50.dat" );

  if( myQeFile.is_open() ){
    while( getline( myQeFile, myString ) ){

      fPmtQeWl.push_back( atof( myString.substr( 0, myString.find(",")).c_str() ) );
      fPmtQe.push_back  ( atof( myString.substr( myString.find(",")+1 ).c_str() ) );

    }
    myQeFile.close();
    DSLog(trace) << "TPC PMT Q.E. loaded." << endlog;  
  }
  else {
    DSLog(trace) << "Couldn't load the TPC PMT Q.E." << endlog;  
  }
  
  // Get the dimensions of the Lumirror sheath around the cryostat
  G4int indexMaxR = 0;
  fCryoSheathR[0] = 0;
  fCryoSheathR[1] = 0;
  fCryoSheathZ[0] = 0;
  fCryoSheathZ[1] = 0;
  //  DSIO::Get()->CloseStreamDSCryostatProfile(); 
  G4int myNumPointsOuterCryo;
  G4double myOuterCryostatZ[30],myOuterCryostatR[30];
  DSIO::Get()->GetStreamDSCryostatProfile() >> myNumPointsOuterCryo;
  for(G4int i = 0; i < myNumPointsOuterCryo; i++ ) { 
    DSIO::Get()->GetStreamDSCryostatProfile() >> myOuterCryostatZ[i]  >> myOuterCryostatR[i];
    myOuterCryostatZ[i] += fCryostatShiftZ;
    if(i > 0 && 
       myOuterCryostatZ[i] < myOuterCryostatZ[i-1] && 
       myOuterCryostatR[i] < myOuterCryostatR[i-1] &&
       myOuterCryostatR[i] > fCryoSheathR[0])
      {
	indexMaxR = i;
	fCryoSheathR[0] = myOuterCryostatR[i];
      }
  }
  DSIO::Get()->CloseStreamDSCryostatProfile();
  fCryoSheathR[0] = myOuterCryostatR[indexMaxR];
  fCryoSheathR[1] = myOuterCryostatR[indexMaxR+1];
  fCryoSheathZ[0] = myOuterCryostatZ[indexMaxR];
  fCryoSheathZ[1] = myOuterCryostatZ[indexMaxR+1];
  DSLog(debugging) << "Cryostat Sheath Dimensions : " << endlog;
  DSLog(debugging) << "\tfCryoSheathR[0] = " << fCryoSheathR[0] << endlog;
  DSLog(debugging) << "\tfCryoSheathR[1] = " << fCryoSheathR[1] << endlog;
  DSLog(debugging) << "\tfCryoSheathZ[0] = " << fCryoSheathZ[0] << endlog;
  DSLog(debugging) << "\tfCryoSheathZ[1] = " << fCryoSheathZ[1] << endlog;
}

void DSParameters::ReadOpticalProperties() {
  G4double hc = 1239.84172; //eV*nm
  
  // Load Detector Properties
  char lumirrorFileName[128] = "../data/detector/lumirrorReflectance.dat"; 
  ifstream myPropertiesFile( lumirrorFileName );
  G4double wavelength, reflectance;

  // Define constant properties
  G4double epss_const_refl = 0.58;
  G4double untreatedss_const_refl = 0.48;
  G4double alfoil_const_refl = 0.8;
  G4double pmtback_const_refl = 0.78;
  G4double lumirror_const_eff = 0.;
  G4double epss_const_eff = 0.;
  G4double untreatedss_const_eff = 0.;
  G4double alfoil_const_eff = 0.;
  G4double pmtback_const_eff = 0.;
  G4double lumirror_const_speclobe = 0;
  G4double epss_const_speclobe = 0.65;
  G4double untreatedss_const_speclobe = 0.65;
  G4double alfoil_const_speclobe = 0.65;
  G4double pmtback_const_speclobe = 1.;
  G4double lumirror_const_specspike = 0;
  G4double epss_const_specspike = 0.0;
  G4double untreatedss_const_specspike = 0.0;
  G4double alfoil_const_specspike = 0.0;
  G4double pmtback_const_specspike = 0.0;
  G4double lumirror_const_backscat = 0;
  G4double epss_const_backscat = 0.;
  G4double untreatedss_const_backscat = 0.;
  G4double alfoil_const_backscat = 0.;
  G4double pmtback_const_backscat = 0.;
  G4double lumirror_const_rindex = 1.;
  G4double epss_const_rindex = 1.;
  G4double untreatedss_const_rindex = 1.;
  G4double alfoil_const_rindex = 1.;
  G4double pmtback_const_rindex = 1.;
  
  // Load properties into arrays
  G4int i = 0, counter = 0;
  
  if( myPropertiesFile.is_open() ){
    while( myPropertiesFile >> wavelength >> reflectance ){
      counter++;
      //if((counter % fLumirrorSkipNEntries) != 0)
      if((counter % 10) != 0)
	continue;
      fLumirror_Energy.push_back( (hc/wavelength)*eV);
      fElectropolishedStainlessSteel_Energy.push_back( (hc/wavelength)*eV);
      fUntreatedStainlessSteel_Energy.push_back( (hc/wavelength)*eV);
      fAluminumFoil_Energy.push_back( (hc/wavelength)*eV);
      fPMTBack_Energy.push_back( (hc/wavelength)*eV);

      DSLog(debugging) << i << " WL = " << wavelength << " nm\tLUM : " << reflectance/100 << endlog;

      fLumirror_Reflectance.push_back( reflectance/100.);
      fElectropolishedStainlessSteel_Reflectance.push_back( epss_const_refl);
      fUntreatedStainlessSteel_Reflectance.push_back( untreatedss_const_refl);
      fAluminumFoil_Reflectance.push_back( alfoil_const_refl);
      fPMTBack_Reflectance.push_back( pmtback_const_refl);

      fLumirror_Efficiency.push_back( lumirror_const_eff);
      fElectropolishedStainlessSteel_Efficiency.push_back( epss_const_eff);
      fUntreatedStainlessSteel_Efficiency.push_back( untreatedss_const_eff);
      fAluminumFoil_Efficiency.push_back( alfoil_const_eff);
      fPMTBack_Efficiency.push_back( pmtback_const_eff);
      
      fLumirror_SpecularLobe.push_back( lumirror_const_speclobe);
      fElectropolishedStainlessSteel_SpecularLobe.push_back( epss_const_speclobe);
      fUntreatedStainlessSteel_SpecularLobe.push_back( untreatedss_const_speclobe);
      fAluminumFoil_SpecularLobe.push_back( alfoil_const_speclobe);
      fPMTBack_SpecularLobe.push_back( pmtback_const_speclobe);
      
      fLumirror_SpecularSpike.push_back( lumirror_const_specspike);
      fElectropolishedStainlessSteel_SpecularSpike.push_back( epss_const_specspike);
      fUntreatedStainlessSteel_SpecularSpike.push_back( untreatedss_const_specspike);
      fAluminumFoil_SpecularSpike.push_back( alfoil_const_specspike);
      fPMTBack_SpecularSpike.push_back( pmtback_const_specspike);
      
      fLumirror_Backscatter.push_back( lumirror_const_backscat);
      fElectropolishedStainlessSteel_Backscatter.push_back( epss_const_backscat);
      fUntreatedStainlessSteel_Backscatter.push_back( untreatedss_const_backscat);
      fAluminumFoil_Backscatter.push_back( alfoil_const_backscat);
      fPMTBack_Backscatter.push_back( pmtback_const_backscat);
      
      fLumirror_Rindex.push_back( lumirror_const_rindex);
      fElectropolishedStainlessSteel_Rindex.push_back( epss_const_rindex);
      fUntreatedStainlessSteel_Rindex.push_back( untreatedss_const_rindex);
      fAluminumFoil_Rindex.push_back( alfoil_const_rindex);
      fPMTBack_Rindex.push_back( pmtback_const_rindex);
      
      i++;
    }
    myPropertiesFile.close();

    // Veto PMT QEs
    
    G4bool useTPCPMT = false;
    char vqeFileName[128];
    if(useTPCPMT)
      sprintf(vqeFileName,"../data/detector/PMTTPC.dat");
    else
      sprintf(vqeFileName, "../data/detector/vpmtQE.dat");

    myPropertiesFile.open(vqeFileName);
    G4double qe, qe_adjusted, refl;

    // Define constant properties
    G4double vpc_const_refl = 0.1;//0.2;
    G4double vpc_const_eff = 0.;
    G4double vpc_const_speclobe = 1.;
    G4double vpc_const_specspike = 0.;
    G4double vpc_const_backscat = 0.;
    G4double vpc_const_rindex = 1.;

    if(useTPCPMT)
      fVPmtMaxQe = fPmtMaxQe;
    fVPmtMaxQe_adjusted = fVPmtMaxQe/(1-vpc_const_refl*(1-fVPmtMaxQe));

    // Load properties into arrays
    i = 0;
    if( myPropertiesFile.is_open() ){
      while( myPropertiesFile >> wavelength >> qe )
	{
	  if(!useTPCPMT)
	    qe /= 100.;
	  qe_adjusted = qe/(1-vpc_const_refl*(1-qe));
	  if(!useTPCPMT)
	    {
	      qe /= fVPmtMaxQe;
	      qe_adjusted /= fVPmtMaxQe_adjusted;
	    }
	  refl = (1-qe)*vpc_const_refl;
	  refl *= fVPmtMaxQe_adjusted;

	  fVPmtQe.push_back(qe_adjusted);
	  fVPmtQeWl.push_back(wavelength*nm);
	  fVPhotocathode_Energy.push_back((hc/wavelength)*eV);
	  fVPhotocathode_Reflectance.push_back(refl);
	  fVPhotocathode_Efficiency.push_back(vpc_const_eff);
	  fVPhotocathode_SpecularLobe.push_back(vpc_const_speclobe);
	  fVPhotocathode_SpecularSpike.push_back(vpc_const_specspike);
	  fVPhotocathode_Backscatter.push_back(vpc_const_backscat);
	  fVPhotocathode_Rindex.push_back(vpc_const_rindex);

	  DSLog(debugging) << "VETO PMT QE : WL = " << wavelength << " nm\tQE = " << qe << "\tQE' = " << qe_adjusted << "\tr = " << refl << endlog;

	  i++;
	}
    }
    else 
      DSLog(warning) << "WARNING : COULD NOT OPEN VETO PMT QE DATA FILE " << vqeFileName << endlog;
    myPropertiesFile.close();

    // TPC PMT QEs
    char tqeFileName[128] = "../data/detector/PMTTPC.dat"; // "../data/detector/DSPmtQeDS50.dat";
    myPropertiesFile.open(tqeFileName);
    // Define constant properties
    G4double tpc_const_refl = 0.2;
    G4double tpc_const_eff = 0.;
    G4double tpc_const_speclobe = 1.;
    G4double tpc_const_specspike = 0.;
    G4double tpc_const_backscat = 0.;
    G4double tpc_const_rindex = 1.;
    fPmtMaxQe_adjusted = fPmtMaxQe/(1-tpc_const_refl*(1-fPmtMaxQe));
    // Load properties into arrays
    i = 0;
    if( myPropertiesFile.is_open() ){
      while( myPropertiesFile >> wavelength >> qe )
	{
	  qe_adjusted = fPmtMaxQe*qe/(1-tpc_const_refl*(1-fPmtMaxQe*qe));
	  qe_adjusted /= fPmtMaxQe_adjusted;
	  refl = (1-qe)*tpc_const_refl;
	  refl *= fPmtMaxQe_adjusted;

	  //fPmtQe.push_back(qe_adjusted);
	  //fPmtQeWl.push_back(wavelength*nm);
	  fPhotocathode_Energy.push_back((hc/wavelength)*eV);
	  fPhotocathode_Reflectance.push_back(refl);
	  fPhotocathode_Efficiency.push_back(tpc_const_eff);
	  fPhotocathode_SpecularLobe.push_back(tpc_const_speclobe);
	  fPhotocathode_SpecularSpike.push_back(tpc_const_specspike);
	  fPhotocathode_Backscatter.push_back(tpc_const_backscat);
	  fPhotocathode_Rindex.push_back(tpc_const_rindex);

	  DSLog(debugging) << "TPC PMT QE : WL = " << wavelength << " nm\tQE = " << qe << "\tQE' = " << qe_adjusted << "\tr = " << refl << endlog;

	  i++;
	}
    }
    else
      DSLog(warning) << "WARNING : COULD NOT OPEN TPC PMT DATA FILE " << tqeFileName << endlog;
    myPropertiesFile.close();
    
    /*
    // Load properties into material properties table
    fLumirrorMPT = new G4MaterialPropertiesTable();
    fLumirrorMPT->AddProperty("SPECULARLOBECONSTANT",
			      fLumirror_Energy,
			      fLumirror_SpecularLobe,
			      fLumirrorNumEntries);
    fLumirrorMPT->AddProperty("SPECULARSPIKECONSTANT",
			      fLumirror_Energy,
			      fLumirror_SpecularSpike,
			      fLumirrorNumEntries);
    fLumirrorMPT->AddProperty("BACKSCATTERCONSTANT",
			      fLumirror_Energy,
			      fLumirror_Backscatter,
			      fLumirrorNumEntries);
    fLumirrorMPT->AddProperty("REFLECTIVITY",
			      fLumirror_Energy,
			      fLumirror_Reflectance,
			      fLumirrorNumEntries)->SetSpline(true);
    fLumirrorMPT->AddProperty("EFFICIENCY",
			      fLumirror_Energy,
			      fLumirror_Efficiency,
			      fLumirrorNumEntries);
    fLumirrorMPT->AddProperty("RINDEX",
			      fLumirror_Energy,
			      fLumirror_Rindex,
			      fLumirrorNumEntries);
    
    fElectropolishedStainlessSteelMPT = new G4MaterialPropertiesTable();
    fElectropolishedStainlessSteelMPT->AddProperty("SPECULARLOBECONSTANT",
						   fElectropolishedStainlessSteel_Energy,
						   fElectropolishedStainlessSteel_SpecularLobe,
						   fLumirrorNumEntries);
    fElectropolishedStainlessSteelMPT->AddProperty("SPECULARSPIKECONSTANT",
						   fElectropolishedStainlessSteel_Energy,
						   fElectropolishedStainlessSteel_SpecularSpike,
						   fLumirrorNumEntries);
    fElectropolishedStainlessSteelMPT->AddProperty("BACKSCATTERCONSTANT",
						   fElectropolishedStainlessSteel_Energy,
						   fElectropolishedStainlessSteel_Backscatter,
						   fLumirrorNumEntries);
    fElectropolishedStainlessSteelMPT->AddProperty("REFLECTIVITY",
						   fElectropolishedStainlessSteel_Energy,
						   fElectropolishedStainlessSteel_Reflectance,
						   fLumirrorNumEntries)->SetSpline(true);
    fElectropolishedStainlessSteelMPT->AddProperty("EFFICIENCY",
						   fElectropolishedStainlessSteel_Energy,
						   fElectropolishedStainlessSteel_Efficiency,
						   fLumirrorNumEntries);
    fElectropolishedStainlessSteelMPT->AddProperty("RINDEX",
						   fElectropolishedStainlessSteel_Energy,
						   fElectropolishedStainlessSteel_Rindex,
						   fLumirrorNumEntries);

    fUntreatedStainlessSteelMPT = new G4MaterialPropertiesTable();
    fUntreatedStainlessSteelMPT->AddProperty("SPECULARLOBECONSTANT",
					      fUntreatedStainlessSteel_Energy,
					      fUntreatedStainlessSteel_SpecularLobe,
					      fLumirrorNumEntries);
    fUntreatedStainlessSteelMPT->AddProperty("SPECULARSPIKECONSTANT",
					      fUntreatedStainlessSteel_Energy,
					      fUntreatedStainlessSteel_SpecularSpike,
					      fLumirrorNumEntries);
    fUntreatedStainlessSteelMPT->AddProperty("BACKSCATTERCONSTANT",
					      fUntreatedStainlessSteel_Energy,
					      fUntreatedStainlessSteel_Backscatter,
					      fLumirrorNumEntries);
    fUntreatedStainlessSteelMPT->AddProperty("REFLECTIVITY",
					      fUntreatedStainlessSteel_Energy,
					      fUntreatedStainlessSteel_Reflectance,
					      fLumirrorNumEntries)->SetSpline(true);
    fUntreatedStainlessSteelMPT->AddProperty("EFFICIENCY",
					      fUntreatedStainlessSteel_Energy,
					      fUntreatedStainlessSteel_Efficiency,
					      fLumirrorNumEntries);
    fUntreatedStainlessSteelMPT->AddProperty("RINDEX",
					      fUntreatedStainlessSteel_Energy,
					      fUntreatedStainlessSteel_Rindex,
					      fLumirrorNumEntries);

    fAluminumFoilMPT = new G4MaterialPropertiesTable();
    fAluminumFoilMPT->AddProperty("SPECULARLOBECONSTANT",
				  fAluminumFoil_Energy,
				  fAluminumFoil_SpecularLobe,
				  fLumirrorNumEntries);
    fAluminumFoilMPT->AddProperty("SPECULARSPIKECONSTANT",
				  fAluminumFoil_Energy,
				  fAluminumFoil_SpecularSpike,
				  fLumirrorNumEntries);
    fAluminumFoilMPT->AddProperty("BACKSCATTERCONSTANT",
				  fAluminumFoil_Energy,
				  fAluminumFoil_Backscatter,
				  fLumirrorNumEntries);
    fAluminumFoilMPT->AddProperty("REFLECTIVITY",
				  fAluminumFoil_Energy,
				  fAluminumFoil_Reflectance,
				  fLumirrorNumEntries);
    fAluminumFoilMPT->AddProperty("EFFICIENCY",
				  fAluminumFoil_Energy,
				  fAluminumFoil_Efficiency,
				  fLumirrorNumEntries);
    fAluminumFoilMPT->AddProperty("RINDEX",
				  fAluminumFoil_Energy,
				  fAluminumFoil_Rindex,
				  fLumirrorNumEntries);

    fVPhotocathodeMPT = new G4MaterialPropertiesTable();
    fVPhotocathodeMPT->AddProperty("SPECULARLOBECONSTANT",
				   fVPhotocathode_Energy,
				   fVPhotocathode_SpecularLobe,
				   fVQENumEntries);
    fVPhotocathodeMPT->AddProperty("SPECULARSPIKECONSTANT",
				   fVPhotocathode_Energy,
				   fVPhotocathode_SpecularSpike,
				   fVQENumEntries);
    fVPhotocathodeMPT->AddProperty("BACKSCATTERCONSTANT",
				   fVPhotocathode_Energy,
				   fVPhotocathode_Backscatter,
				   fVQENumEntries);
    fVPhotocathodeMPT->AddProperty("REFLECTIVITY",
				   fVPhotocathode_Energy,
				   fVPhotocathode_Reflectance,
				   fVQENumEntries);
    fVPhotocathodeMPT->AddProperty("EFFICIENCY",
				   fVPhotocathode_Energy,
				   fVPhotocathode_Efficiency,
				   fVQENumEntries);
    fVPhotocathodeMPT->AddProperty("RINDEX",
				   fVPhotocathode_Energy,
				   fVPhotocathode_Rindex,
				   fVQENumEntries);

    fPMTBackMPT = new G4MaterialPropertiesTable();
    fPMTBackMPT->AddProperty("SPECULARLOBECONSTANT",
			     fPMTBack_Energy,
			     fPMTBack_SpecularLobe,
			     fLumirrorNumEntries);
    fPMTBackMPT->AddProperty("SPECULARSPIKECONSTANT",
			     fPMTBack_Energy,
			     fPMTBack_SpecularSpike,
			     fLumirrorNumEntries);
    fPMTBackMPT->AddProperty("BACKSCATTERCONSTANT",
			     fPMTBack_Energy,
			     fPMTBack_Backscatter,
			     fLumirrorNumEntries);
    fPMTBackMPT->AddProperty("REFLECTIVITY",
			     fPMTBack_Energy,
			     fPMTBack_Reflectance,
			     fLumirrorNumEntries);
    fPMTBackMPT->AddProperty("EFFICIENCY",
			     fPMTBack_Energy,
			     fPMTBack_Efficiency,
			     fLumirrorNumEntries);
    fPMTBackMPT->AddProperty("RINDEX",
			     fPMTBack_Energy,
			     fPMTBack_Rindex,
			     fLumirrorNumEntries);
    */
    DSLog(trace) << "Detector Properties loaded" << endlog;  
  }
  else {
    DSLog(trace) << "Couldn't load Lumirror Properties file " << lumirrorFileName << endlog;  
  }

}


/////////////////////////////////////////////////////////
//// Method returning the TPC QE , given the WL in nm
/////////////////////////////////////////////////////////
G4double   DSParameters::GetTPCQE(G4double myPhotonWL) {
  
  //  if( myPhotonWL < fPmtQeWl[1] )       return fPmtQe[1]; 
  if( myPhotonWL < fPmtQeWl[1] )       return 0; //return 0 for the QE if WL<200 nm 
  if( myPhotonWL > fPmtQeWl.back()  )  return fPmtQe.back();  
  
  G4int    myIndex = 0;
  G4double _WL0 = 0;
  G4double _WL1 = 0;
  
  for(std::vector<G4double>::iterator it = fPmtQeWl.begin(); it != fPmtQeWl.end(); ++it) {
    _WL1 = *it;
    if(*it > myPhotonWL) break;
    _WL0 = *it;
    myIndex++;
  }
  
  G4double _QE0 = fPmtQe[myIndex - 1]*fPmtMaxQe;
  G4double _QE1 = fPmtQe[myIndex]*fPmtMaxQe;
  // interpolation  
  G4double _m  =  (_QE0 - _QE1)/(_WL0 - _WL1);
  G4double _q  =  _QE1 - _m*_WL1;

  return _m*myPhotonWL + _q;
    
}

std::vector<G4double> DSParameters::ConvertTo4DoubleVector(const char* st) {
  G4double v1;
  G4double v2;
  G4double v3;
  G4double v4;
  std::istringstream is(st);
  is >> v1 >> v2 >> v3 >> v4;
  std::vector<G4double> tmp ;
  tmp.push_back(v1);
  tmp.push_back(v2);
  tmp.push_back(v3);
  tmp.push_back(v4);
  return tmp;
}

/*
 * $Log: DSParameters.cc,v $
 * Revision 1.10  2015/01/17 11:31:53  pagnes
 * PAr model added form optical tuning
 *
 * Revision 1.9  2015/01/07 16:45:53  pagnes
 * changed veto optical properties format from arrays to vectors
 *
 * Revision 1.8  2014/12/22 22:42:04  jbrodsky
 * Fixed bug in filling veto PMT locations from data file. The loop was going one over the limit of the relevant arrays (trying to fill entry 110 in an array of size 110) which caused bad codebehavior and segfaults.
 *
 * Revision 1.7  2014/11/20 17:08:55  dfranco
 * QE of PMTs = 0 if wavelength<200 nm
 *
 * Revision 1.6  2014/11/13 16:47:04  dfranco
 * removed variables which were creating conflicts with the previous version of g4ds10
 *
 * Revision 1.5  2014/11/05 15:47:13  pagnes
 * temporary optics tuning
 *
 * Revision 1.4  2014/10/13 18:43:57  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.3  2014/07/16 08:23:04  pagnes
 * QE scaling to 1.0 added (/ds/manager/fast_simulation xxx)
 *
 * Revision 1.2  2014/06/03 13:31:35  meyers
 * Migrate TPC grid and ITO optics updates to g4ds10
 *
 * Revision 1.1  2014/05/07 12:21:04  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.21  2014/03/19 16:37:26  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.20  2014/03/11 16:49:56  dfranco
 * read optical parameters from external file (DSOptics.dat)
 *
 * Revision 1.19  2013/08/05 12:25:13  swesterd
 * Added G4OpBoundaryProcess.hh so that light can pass through some optical boundaries
 *
 * Revision 1.18  2013/08/05 03:13:52  swesterd
 * some fine tuning of bscint and veto parameters
 *
 * Revision 1.17  2013/06/24 13:43:29  dfranco
 * cout removed
 *
 * Revision 1.16  2013/06/24 13:05:55  dfranco
 * TPC QE values were filled twice: once in the standard way, the second deconvoluting the reflections. The second has been commented
 *
 * Revision 1.15  2013/06/22 07:21:21  dfranco
 * Fixed a bug in the photoelectron absorption in DS50. Back to the previous QE method
 *
 * Revision 1.14  2013/06/19 18:35:28  swesterd
 * added DSScintCelll and made tpc PMTs' QE and reflections work like veto PMTs
 *
 * Revision 1.13  2013/06/05 23:03:32  swesterd
 * moved optical boundary MPTs to DSMaterial and gave the trunks optical boundary properties consistent with untreated stainless steel
 *
 * Revision 1.12  2013/06/04 23:38:49  swesterd
 * Changed the length of the Lumirror reflectance array to be an adjustable size, set it to 60 elements
 *
 * Revision 1.11  2013/06/04 14:11:36  dfranco
 * Added a function returning the TPC QE in DSParameters, applied in the tracking action
 *
 * Revision 1.10  2013/05/27 23:59:02  swesterd
 * added a (currently commented out) Lumirror sheath to the cryostat and introduced DSOpBoundaryProcess to try to figure out why the boundaries are being screwy, with some edits so that it can handle constant and vector properties with freaking out
 *
 * Revision 1.9  2013/05/25 07:58:23  swesterd
 * Got the veto PMT optical boundaries all working along with photocathode optical properties, added PMT quantum efficiency to DSTrackingAction, and added a function to DSTrackingAction that locates and quadratically interpolates points in data, for getting useful QEs
 *
 * Revision 1.8  2013/05/07 23:06:30  swesterd
 * added optical boundaries and Lumirror in the veto
 *
 * Revision 1.7  2013/05/06 14:59:52  perassos
 * Updates on the TPC surface properties and geometry
 *
 * Revision 1.6  2013/05/01 08:20:24  swesterd
 * added boron-loaded scintillator optical properties
 *
 * Revision 1.5  2013/04/26 16:17:19  perassos
 * TPC PMT QE added
 *
 * Revision 1.4  2013/04/19 16:10:09  perassos
 * Added ITO and TPB layers
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
