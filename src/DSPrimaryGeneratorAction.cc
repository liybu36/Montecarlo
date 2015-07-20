#include "DSPrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSPrimaryGeneratorActionMessenger.hh"



DSPrimaryGeneratorAction::DSPrimaryGeneratorAction() {
   fMessenger = new DSPrimaryGeneratorActionMessenger(this);
}

DSPrimaryGeneratorAction::~DSPrimaryGeneratorAction(){
    delete fMessenger;
    delete generator;
}

void DSPrimaryGeneratorAction::GeneratePrimaries(G4Event *event) {
  
  generator->DSGeneratePrimaries(event);  
  

}

/*
 * $Log: DSPrimaryGeneratorAction.cc,v $
 * Revision 1.1  2014/05/07 12:21:04  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
