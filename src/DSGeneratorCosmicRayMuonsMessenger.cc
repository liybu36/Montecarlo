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
#include "DSGeneratorCosmicRayMuonsMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSGeneratorCosmicRayMuons.hh"
#include "DSLogger.hh"

#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>               
#include "G4Tokenizer.hh"

DSGeneratorCosmicRayMuonsMessenger::DSGeneratorCosmicRayMuonsMessenger
(DSGeneratorCosmicRayMuons* fPtclGun)  {
  generator = fPtclGun;
  fDirectory = new G4UIdirectory("/ds/generator/cosmicray/");
  fDirectory->SetGuidance("Control of cosmic rays event generator");

  // height of source
  heightCmd = new 
    G4UIcmdWithADoubleAndUnit("/ds/generator/cosmicray/height",this);
  heightCmd->SetGuidance("Set the z position of the muon shower");
  heightCmd->SetGuidance("Default value: 8.1 m");
  heightCmd->SetParameterName("Height",true,true);
  heightCmd->SetDefaultUnit("cm");
  heightCmd->SetUnitCandidates("mm cm m km");

  // radius of source
  radiusCmd = new 
    G4UIcmdWithADoubleAndUnit("/ds/generator/cosmicray/radius",this);
  radiusCmd->SetGuidance("Set the radius of the muon shower");
  radiusCmd->SetGuidance("Default value: 8.0 m");
  radiusCmd->SetParameterName("Radius",true,true);
  radiusCmd->SetDefaultUnit("cm");
  radiusCmd->SetUnitCandidates("mm cm m km");

  // depth of the Laboratory
  depthCmd = new 
    G4UIcmdWithADoubleAndUnit("/ds/generator/cosmicray/depth",this);
  depthCmd->SetGuidance("Set the depth of the Laboratory");
  depthCmd->SetGuidance("expressed in water equivalent");
  depthCmd->SetGuidance("Default value: 3700 m (Gran Sasso depth)");
  depthCmd->SetParameterName("Depth",true,true);
  depthCmd->SetDefaultUnit("cm");
  depthCmd->SetUnitCandidates("cm m km");
  //depthCmd->SetRange("Depth>=1*km && Depth <=10*km");

  //index of the muon power spectrum
  indexCmd = new 
    G4UIcmdWithADouble("/ds/generator/cosmicray/index",this);
  indexCmd->SetGuidance("Set the spectral index of muons");
  indexCmd->SetGuidance("2.0 --> Exotic sources ");
  indexCmd->SetGuidance("2.7 --> Prompt sources (e.g. charm decay)");
  indexCmd->SetGuidance("3.7 --> Standard spectrum ");
  indexCmd->SetGuidance("Default value: 3.7");
  indexCmd->SetParameterName("Index",true,true);
  //indexCmd->SetRange("Index == 2.0 || Index == 2.7 || Index == 3.7"); 
  
  // lower edge of the energy spectrum
  energyLowCmd = new 
    G4UIcmdWithADoubleAndUnit("/ds/generator/cosmicray/energyLow",this);
  energyLowCmd->SetGuidance("Set the lower edge of the energy spectrum");
  energyLowCmd->SetGuidance("(minimum 100 MeV, maximum 10 GeV)");
  energyLowCmd->SetGuidance("Default value: 1 GeV");
  energyLowCmd->SetParameterName("energyLow",true,true);
  energyLowCmd->SetDefaultUnit("GeV");
  energyLowCmd->SetUnitCandidates("MeV GeV TeV");
  //energyLowCmd->SetRange("energyLow>=100*MeV && energyLow <=10*GeV");

  // upper edge of the energy spectrum
  energyUpCmd = new 
    G4UIcmdWithADoubleAndUnit("/ds/generator/cosmicray/energyUp",this);
  energyUpCmd->SetGuidance("Set the upper edge of the energy spectrum");
  energyUpCmd->SetGuidance("(minimum 100 GeV, maximum 50 TeV)");
  energyUpCmd->SetGuidance("Default value: 10 TeV");
  energyUpCmd->SetParameterName("energyUp",true,true);
  energyUpCmd->SetDefaultUnit("GeV");
  energyUpCmd->SetUnitCandidates("MeV GeV TeV");
  //energyUpCmd->SetRange("energyLow>=100*GeV && energyLow <=50*TeV");

  // name of the file containing the angular spectrum
  fileNameCmd = new 
    G4UIcmdWithAString("/ds/generator/cosmicray/filename",this);
  fileNameCmd->SetGuidance("Name of the file containing the angular spectrum");
  fileNameCmd->SetGuidance("Format: costheta phi pdf");
  fileNameCmd->SetGuidance("(evenly-spaced in costheta)");
  fileNameCmd->SetGuidance("Default: zenith_azimuth.dat");
}


DSGeneratorCosmicRayMuonsMessenger::~DSGeneratorCosmicRayMuonsMessenger() 
{
  delete fDirectory;
  delete heightCmd;
  delete radiusCmd;
  delete depthCmd;
  delete indexCmd;
  delete energyLowCmd;
  delete energyUpCmd;
  delete fileNameCmd;
}

void DSGeneratorCosmicRayMuonsMessenger::SetNewValue(G4UIcommand *command, 
						       G4String newValues) 
{
  if (command == heightCmd)
    {
      generator->SetHalfZ(heightCmd->ConvertToDimensionedDouble(newValues));
    }
  else if (command == radiusCmd) 
    {
      generator->SetRadius(radiusCmd->ConvertToDimensionedDouble(newValues));
    }
  else if (command == depthCmd) 
    {
      double depth(depthCmd->ConvertToDimensionedDouble(newValues));
      if(depth<1*km || depth>10*km)
	{
	  DSLog(error) << "Depth "<<depth<<" out of range! Exit now" << endlog;
	  DSLog(fatal) <<""<< endlog;
	}     
      generator->SetRockDepth(depth);
    }
  else if (command == indexCmd)
    {
      double index(indexCmd->ConvertToDouble(newValues));
      if(index!=2.0 && index!=2.7 && index!=3.7)
	{
	  DSLog(error) << " Spectral index "<<index<<" not permitted! Exit now" << endlog;
	  DSLog(fatal) <<""<< endlog;
	}     
      generator->SetSpectralIndex(index);
    }
   else if (command == energyLowCmd)
    {
      double enelow(energyLowCmd->ConvertToDimensionedDouble(newValues));
      if(enelow<100*MeV || enelow>10*GeV)
	{
	  DSLog(error) << "Energy lower limit "<<enelow<<" out of range! Exit now" << endlog;
	  DSLog(fatal) <<""<< endlog;
	}     
      generator->SetEnergyInf(enelow);
    }
  else if (command == energyUpCmd)
    {
      double eneup(energyUpCmd->ConvertToDimensionedDouble(newValues));
      if(eneup<100*GeV || eneup>50*TeV)
	{
	  DSLog(error) << "Energy upper limit "<<eneup<<" out of range! Exit now" << endlog;
	  DSLog(fatal) <<""<< endlog;
	}     
      generator->SetEnergySup(eneup);
    }
  else if (command == fileNameCmd)
    {
      generator->SetFileName(newValues);
      DSLog(trace) << "File of angular spectrum changed to " << 
	generator->GetFileName() << endlog;
      DSLog(trace) << "Be sure that the format is correct!" << endlog;
    }
}




/*
 * Revision 1.2  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 */
