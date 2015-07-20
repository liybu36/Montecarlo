// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// *********--***********************************************************
//

#include <iostream>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "DSVGeneratorMessenger.hh"
#include "G4UIcmdWithADouble.hh"
#include "DSVGenerator.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "DSLogger.hh"
using namespace std;

DSVGeneratorMessenger::DSVGeneratorMessenger(DSVGenerator* gen){
  generator = gen;
  fDirectory = new G4UIdirectory("/ds/generator/");
  fDirectory->SetGuidance("Control of DSG4Gun event generator");

  fPositionCmd = new G4UIcmdWith3VectorAndUnit("/ds/generator/position",this);
  fPositionCmd->SetGuidance("Set the gun position");
  fPositionCmd->SetUnitCategory("Length");
  fPositionCmd->SetDefaultUnit("cm");
  fPositionCmd->SetUnitCandidates("mm cm m");
  
  fDirectionCmd = new G4UIcmdWith3Vector("/ds/generator/direction",this);
  fDirectionCmd->SetGuidance("Set the gun direction");
  
  fEnergyCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/energy",this);
  fEnergyCmd->SetGuidance("Set the gun energy");
  fEnergyCmd->SetUnitCategory("Energy");
  fEnergyCmd->SetDefaultUnit("MeV");
  fEnergyCmd->SetUnitCandidates("eV keV MeV GeV");
  
  fSphereBulkCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/sphere_radius",this);
  fSphereBulkCmd->SetGuidance("Bulk radius");
  fSphereBulkCmd->SetUnitCategory("Length");
  fSphereBulkCmd->SetDefaultUnit("cm");
  fSphereBulkCmd->SetUnitCandidates("mm cm m");
  
  fSphereBulkRadMinCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/sphere_radius_min",this);
  fSphereBulkRadMinCmd->SetGuidance("Bulk minimum radius");
  fSphereBulkRadMinCmd->SetGuidance("Default: 0");
  fSphereBulkRadMinCmd->SetUnitCategory("Length");
  fSphereBulkRadMinCmd->SetDefaultUnit("cm");
  fSphereBulkRadMinCmd->SetUnitCandidates("mm cm m");

  fSphereSurfCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/surface_radius",this);
  fSphereSurfCmd->SetGuidance("Surface radius");
  fSphereSurfCmd->SetUnitCategory("Length");
  fSphereSurfCmd->SetDefaultUnit("cm");
  fSphereSurfCmd->SetUnitCandidates("mm cm m");

  fSphereCentreCmd = new G4UIcmdWith3VectorAndUnit("/ds/generator/set_center",this);
  fSphereCentreCmd->SetGuidance("Set the sphere position");
  fSphereCentreCmd->SetGuidance("Default: 0. 0. 0. cm");
  fSphereCentreCmd->SetUnitCategory("Length");
  fSphereCentreCmd->SetDefaultUnit("cm");
  fSphereCentreCmd->SetUnitCandidates("mm cm m");

 

  fEnergyDisTypeCmd = new G4UIcmdWithAString("/ds/generator/dist_energy",this);
  G4String candidates = "Mono Lin Pow Exp Gauss Brem BBody Cdg";  
  fEnergyDisTypeCmd->SetCandidates(candidates);
 
  fSetEminCmd      = new G4UIcmdWithADoubleAndUnit("/ds/generator/emin",this);
  fSetEminCmd->SetDefaultUnit("MeV");
  fSetEminCmd->SetUnitCandidates("eV keV MeV GeV");
  
  fSetEmaxCmd      = new G4UIcmdWithADoubleAndUnit("/ds/generator/emax",this);
  fSetEmaxCmd->SetDefaultUnit("MeV");
  fSetEmaxCmd->SetUnitCandidates("eV keV MeV GeV");
  
  fSetAlphaCmd     = new G4UIcmdWithADouble("/ds/generator/alpha",this);
  
  fSetTempCmd      = new G4UIcmdWithADoubleAndUnit("/ds/generator/temp",this);
  fSetTempCmd->SetDefaultUnit("kelvin");
  fSetTempCmd->SetUnitCandidates("kelvin");
  
  fSetEzeroCmd     = new G4UIcmdWithADoubleAndUnit("/ds/generator/ezero",this);
  fSetEzeroCmd->SetDefaultUnit("MeV");
  fSetEzeroCmd->SetUnitCandidates("eV keV MeV GeV");
  
  fSetGradientCmd  = new G4UIcmdWithADouble("/ds/generator/gradient",this);
  
  fSetInterCeptCmd = new G4UIcmdWithADouble("/ds/generator/intercept",this);
  
 
  fParticleCmd = new G4UIcmdWithAString("/ds/generator/particle",this);
  fParticleCmd->SetGuidance("Set the gun type");

  // confine to volume
  fConfineCmd = new G4UIcmdWithAString("/ds/generator/confine",this);
  fConfineCmd->SetGuidance("Confine source to volume (NULL to unset).");
  fConfineCmd->SetGuidance("usage: confine VolName");
  fConfineCmd->SetParameterName("VolName",true,true);
  fConfineCmd->SetDefaultValue("NULL");

  fDisPosTypeCmd = new G4UIcmdWithAString("/ds/generator/postype",this);
  candidates = "Point Plane Surface Volume";
  fDisPosTypeCmd->SetCandidates(candidates);
  
  fDisPosShapeCmd = new G4UIcmdWithAString("/ds/generator/posshape",this);
  candidates = "Square Circle Ellipse Rectangle Sphere Ellipsoid Cylinder Parallelepiped";
  fDisPosShapeCmd->SetCandidates(candidates);
  
  fSetHalfXCmd     = new G4UIcmdWithADoubleAndUnit("/ds/generator/sethalfX",this);
  fSetHalfXCmd->SetUnitCandidates("mm cm m");
  
  fSetHalfYCmd     = new G4UIcmdWithADoubleAndUnit("/ds/generator/sethalfY",this);
  fSetHalfYCmd->SetUnitCandidates("mm cm m");
 
  fSetHalfZCmd     = new G4UIcmdWithADoubleAndUnit("/ds/generator/sethalfZ",this);
  fSetHalfZCmd->SetUnitCandidates("mm cm m");
  
  fRadiusCmd	 = new G4UIcmdWithADoubleAndUnit("/ds/generator/setradius",this);
  fRadiusCmd->SetUnitCandidates("mm cm m");
  
  fRadius0Cmd	 = new G4UIcmdWithADoubleAndUnit("/ds/generator/setradius0",this);
  fRadius0Cmd->SetUnitCandidates("mm cm m");
  
  //define Angular distribution of the source
  fSetAngDistTypeCmd = new G4UIcmdWithAString("/ds/generator/ang_type",this);
  candidates = "iso cos user planar beam1d beam2d focused";
  fSetAngDistTypeCmd->SetCandidates(candidates);

  fSetMinThetaCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/ang_mintheta",this);
  fSetMinThetaCmd->SetUnitCandidates("deg rad");

  fSetMaxThetaCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/ang_maxtheta",this);
  fSetMaxThetaCmd->SetUnitCandidates("deg rad");

  fSetMinPhiCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/ang_minphi",this);
  fSetMinPhiCmd->SetUnitCandidates("deg rad");

  fSetMaxPhiCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/ang_maxphi",this);
  fSetMaxPhiCmd->SetUnitCandidates("deg rad");

  fSetBeamSigmaInAngRCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/ang_sigmar",this);
  fSetBeamSigmaInAngRCmd->SetUnitCandidates("deg rad");

  fSetBeamSigmaInAngXCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/ang_sigmax",this);
  fSetBeamSigmaInAngXCmd->SetUnitCandidates("deg rad");

  fSetBeamSigmaInAngYCmd = new G4UIcmdWithADoubleAndUnit("/ds/generator/ang_sigmay",this);
  fSetBeamSigmaInAngYCmd->SetUnitCandidates("deg rad");

  //  fUserDefAngThetaCmd = new G4UIcmdWith3Vector("/ds/generator/ang_rot1",this);
  //  fUserDefAngPhiCmd = new G4UIcmdWith3Vector("/ds/generator/ang_rot2",this);
  fSetAngRot1Cmd = new G4UIcmdWith3Vector("/ds/generator/ang_rot1",this);
  fSetAngRot2Cmd = new G4UIcmdWith3Vector("/ds/generator/ang_rot2",this);

  fSetFocusPointCmd = new G4UIcmdWith3VectorAndUnit("/ds/generator/ang_focuspoint",this);
  
  
  fEnergyFileCmd = new G4UIcmdWithAString("/ds/generator/energyfile",this);
 
  fUniformTPCCmd = new G4UIcmdWithABool("/ds/generator/tpcdistribution",this);

  fUniformGasPocketCmd = new G4UIcmdWithABool("/ds/generator/gaspocketdistribution",this);
  
  fTPCCenterCmd = new G4UIcmdWithABool("/ds/generator/tpccenter",this);

  fGenInCryostatsCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_cryostats",this);

  fGenInLiquidArgonCmd = new G4UIcmdWithAnInteger("/ds/generator/liquidargon",this);


  fGenInTeflonCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_teflon",this);

  fGenInFusedSilicaCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_fused_silica",this);

  fGenInPMTPhotocathodeCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_pmt_photocathode",this);

  fGenInPMTStemCmd = new G4UIcmdWithAnInteger("/ds/generator/bgd_pmt_stem",this);
  
  fGenInG2Cryostats = new G4UIcmdWithABool("/ds/generator/is_G2_cryostat",this);
  
  fGenInHolderSourceCmd = new G4UIcmdWithAnInteger("/ds/generator/holderSource_on",this);

  fNumberOfParticlesCmd = new G4UIcmdWithAnInteger("/ds/generator/numberofparticles",this);


}


DSVGeneratorMessenger::~DSVGeneratorMessenger() {

  delete fDirectory;
  delete fPositionCmd; 
  delete fDirectionCmd;
  delete fEnergyCmd;
  delete fParticleCmd; 
  delete fSphereBulkCmd; 
  delete fSphereSurfCmd;
  delete fSphereCentreCmd;
  delete fConfineCmd;
  delete fEnergyDisTypeCmd;
  delete fSetEminCmd;
  delete fSetEmaxCmd;
  delete fSetAlphaCmd;
  delete fSetTempCmd;
  delete fSetEzeroCmd;
  delete fSetGradientCmd;
  delete fSetInterCeptCmd; 
  delete fDisPosTypeCmd ;
  delete fDisPosShapeCmd ;
  delete fSetHalfXCmd;  
  delete fSetHalfYCmd;  
  delete fSetHalfZCmd;  
  delete fRadiusCmd;  
  delete fRadius0Cmd;  
  delete fEnergyFileCmd;
  delete fGenInCryostatsCmd;
  delete fGenInLiquidArgonCmd;
  delete fGenInTeflonCmd;
  delete fGenInPMTStemCmd;
  delete fGenInFusedSilicaCmd;
  delete fGenInPMTPhotocathodeCmd;
  delete fGenInG2Cryostats;
  delete fGenInHolderSourceCmd;
  delete fNumberOfParticlesCmd;
  delete fUniformGasPocketCmd;

}


void DSVGeneratorMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue) { 
   if (cmd == fPositionCmd) {         
         generator->SetVParticlePosition(fPositionCmd->ConvertToDimensioned3Vector(newValue));
   } else if (cmd == fDirectionCmd){
         generator->SetVParticleDirection(fDirectionCmd->ConvertTo3Vector(newValue));
   } else if (cmd == fEnergyCmd){
         generator->SetVParticleEnergy(fEnergyCmd->ConvertToDimensionedDouble(newValue));
   } else if (cmd == fParticleCmd){
         if( newValue == "Ar40" )
           generator->SetVParticleDefinition(G4IonTable::GetIonTable()->FindIon( 18, 40, 0, 0 ));
         else
           generator->SetVParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(newValue));
   } else  if (cmd == fSphereBulkCmd){
         generator->SetPositionFlag(2);
         generator->SetIsVolumeDistribution(true);
         generator->SetPosDisType("Volume");
         generator->SetPosDisShape("Sphere");
         generator->SetRadius(fSphereBulkCmd->ConvertToDimensionedDouble(newValue));	 
   } else if (cmd == fSphereSurfCmd){
         generator->SetPositionFlag(2);
         generator->SetIsVolumeDistribution(true);
         generator->SetPosDisType("Surface");
         generator->SetPosDisShape("Sphere");
         generator->SetRadius(fSphereSurfCmd->ConvertToDimensionedDouble(newValue));	 
   } else if (cmd == fSphereCentreCmd){
   	 generator->SetCentreCoords(fSphereCentreCmd->ConvertToDimensioned3Vector(newValue));
   } else if (cmd == fSphereBulkRadMinCmd){
   	 generator->SetRadius0(fSphereBulkRadMinCmd->ConvertToDimensionedDouble(newValue));
   } else if (cmd == fConfineCmd){
         generator->ConfineSourceToVolume(newValue);
   } else if (cmd == fEnergyDisTypeCmd){ 
         DSLog(routine) << "G4Gun Energy distribution " <<newValue  << endlog ;
         generator->SetIsEnergyDistribution(true);
	 generator->SetEnergyDisType(newValue);
   } else if (cmd == fSetEminCmd){ 
         DSLog(routine) << "G4Gun Emin = " <<newValue  << endlog ;
         generator->SetEmin(fSetEminCmd->ConvertToDimensionedDouble(newValue));
   } else if (cmd == fSetEmaxCmd){ 
         DSLog(routine) << "G4Gun Emax = " <<newValue  << endlog ;
         generator->SetEmax(fSetEmaxCmd->ConvertToDimensionedDouble(newValue));
   } else if (cmd == fSetAlphaCmd){ 
         DSLog(routine) << "G4Gun Alpha = " <<newValue  << endlog ;
         generator->SetAlpha(fSetAlphaCmd->ConvertToDouble(newValue));
   } else if (cmd == fSetTempCmd){ 
         DSLog(routine) << "G4Gun Temperature = " <<newValue  << endlog ;
         generator->SetTemp(fSetTempCmd->ConvertToDimensionedDouble(newValue));
   } else if (cmd == fSetEzeroCmd){ 
         DSLog(routine) << "G4Gun Ezero = " <<newValue  << endlog ;
         generator->SetEzero(fSetEzeroCmd->ConvertToDimensionedDouble(newValue));
   } else if (cmd == fSetGradientCmd){ 
         DSLog(routine) << "G4Gun Gradient = " <<newValue  << endlog ;
         generator->SetGradient(fSetGradientCmd->ConvertToDouble(newValue));
   } else if (cmd == fSetInterCeptCmd){ 
         DSLog(routine) << "G4Gun Intercept = " <<newValue  << endlog ;
         generator->SetInterCept(fSetInterCeptCmd->ConvertToDouble(newValue));
   } else if (cmd == fDisPosTypeCmd){ 
         DSLog(routine) << "Position Distribution Type = " <<newValue  << endlog ;
         generator->SetPositionFlag(2);
         generator->SetIsVolumeDistribution(true);
         generator->SetPosDisType(newValue);
   } else if (cmd == fDisPosShapeCmd){ 
         DSLog(routine) << "Position Distribution Shape = " <<newValue  << endlog ;
         generator->SetPositionFlag(2);
         generator->SetIsVolumeDistribution(true);
         generator->SetPosDisShape(newValue);
   } else if (cmd == fSetHalfXCmd){ 
         DSLog(routine) << "HalfX = " <<newValue  << endlog ;
         generator->SetHalfX(fSetHalfXCmd->ConvertToDimensionedDouble(newValue));
   } else if (cmd == fSetHalfYCmd){ 
         DSLog(routine) << "HalfY = " <<newValue  << endlog ;
         generator->SetHalfY(fSetHalfYCmd->ConvertToDimensionedDouble(newValue));
   } else if (cmd == fSetHalfZCmd){ 
         DSLog(routine) << "Half Z = " <<newValue  << endlog ;
         generator->SetHalfZ(fSetHalfZCmd->ConvertToDimensionedDouble(newValue));
   } else if (cmd == fRadiusCmd){ 
         DSLog(routine) << "Radius = " <<newValue  << endlog ;
         generator->SetRadius(fRadiusCmd->ConvertToDimensionedDouble(newValue));
   } else if (cmd == fRadius0Cmd){ 
         DSLog(routine) << "Radius0 = " <<newValue  << endlog ;
         generator->SetRadius0(fRadius0Cmd->ConvertToDimensionedDouble(newValue));
   } else if (cmd == fEnergyFileCmd){ 
         DSLog(routine) << "Random Energy File Name  = " <<newValue  << endlog ;
         generator->SetEnergyFileName(newValue);
   } else if (cmd == fUniformTPCCmd){ 
       if(fUniformTPCCmd->ConvertToBool(newValue)) {
         generator->SetIsVolumeDistribution(true);
         DSLog(routine) << "Random Position Distribution in the TPC  = " <<newValue  << endlog ;
         generator->SetUniformTPC();
       }
   } else if (cmd == fUniformGasPocketCmd){ 
       if(fUniformGasPocketCmd->ConvertToBool(newValue)) {
         generator->SetIsVolumeDistribution(true);
         DSLog(routine) << "Random Position Distribution in the Gas Pocket  = " <<newValue  << endlog ;
         generator->SetUniformGasPocket();
       }
   } else if (cmd == fTPCCenterCmd){ 
       if(fTPCCenterCmd->ConvertToBool(newValue)) {
         DSLog(routine) << "Set position or distribution in the TPC center  = " <<newValue  << endlog ;
         generator->SetTPCCenter();
       }
   } else if (cmd == fGenInCryostatsCmd){
         generator->SetPositionFlag(3);
   } else if (cmd == fGenInTeflonCmd){
         generator->SetPositionFlag(4);
   } else if (cmd == fGenInFusedSilicaCmd){
         generator->SetPositionFlag(5);
   } else if (cmd == fGenInPMTPhotocathodeCmd){
         generator->SetPositionFlag(6);
   } else if (cmd == fGenInPMTStemCmd){
         generator->SetPositionFlag(7);
   } else if (cmd == fGenInLiquidArgonCmd){
         generator->SetIsVolumeDistribution(true);
         generator->SetPositionFlag(8);
   } else if (cmd == fGenInG2Cryostats){
         generator->SetIsG2( fGenInG2Cryostats->ConvertToBool(newValue) );
   } else if (cmd == fGenInHolderSourceCmd){
         generator->SetPositionFlag(9);
	 
   } else if (cmd == fNumberOfParticlesCmd){
         generator->SetVNumberOfParticles( fNumberOfParticlesCmd->ConvertToInt(newValue) );
   }
}
/*
 * $Log: DSVGeneratorMessenger.cc,v $
 * Revision 1.9  2015/01/21 10:21:06  pagnes
 * DSLight3 set as default optics, cout removed
 *
 * Revision 1.8  2015/01/14 16:58:37  dfranco
 * added source vial and correspondent commands by Laura and Bianca; manual updated
 *
 * Revision 1.7  2014/11/13 16:22:22  dfranco
 * fixed a problem with liquid argon space distribution
 *
 * Revision 1.6  2014/10/22 14:03:47  dfranco
 * add method to generate events in Liqud Argon
 *
 * Revision 1.5  2014/07/18 13:54:49  dfranco
 * Added a new position generator in the Gas Pocket; calibrated the numnber of photons generated per electron in S2; added a new command to generate multiple particles in the same position
 *
 * Revision 1.4  2014/05/21 13:47:48  dfranco
 * Fixed bug in DSGeneratorRDMDecayGun.cc and warning in DSVGeneratorMessenger.cc
 *
 * Revision 1.3  2014/05/08 11:00:39  pagnes
 * Added generator in PMTs stem
 *
 * Revision 1.2  2014/05/07 14:27:31  dfranco
 * fixed some bugs and added GdScintillator
 *
 * Revision 1.1  2014/05/07 12:21:05  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.13  2014/04/22 09:52:28  perassos
 * Added the option to simulate bgd from the G2 cryostats
 *
 * Revision 1.11  2014/04/11 10:20:39  perassos
 * Added generation in materials
 *
 * Revision 1.10  2014/02/13 12:03:26  dfranco
 * added new commands for spatial distributions and fixed the manual
 *
 * Revision 1.9  2014/01/29 13:13:41  perassos
 * Update of the electric field handling and of the Nuclear Recoils generator
 *
 * Revision 1.8  2014/01/07 14:11:27  perassos
 * Added the command to simulate Ar recoils
 *
 * Revision 1.7  2013/10/15 16:39:39  dfranco
 * added the possibility to define the energy spectrum reading an ascii file
 *
 * Revision 1.6  2013/08/20 03:25:52  swesterd
 * added G2 TPC geoemtry (not complete) and added monoenergetic energy distribution to generator
 *
 * Revision 1.5  2013/06/22 06:55:45  dfranco
 * change command name from sphere_origin to set_center for spatial distributions
 *
 * Revision 1.4  2013/05/02 16:52:36  dfranco
 * Added Radius and Radius0 methods to the position distributions
 *
 * Revision 1.3  2013/05/02 16:36:13  dfranco
 * Added several position distributions
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
