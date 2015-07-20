#ifndef DSGeneratorRDMUIcmdWithNucleusAndUnit_h
#define DSGeneratorRDMUIcmdWithNucleusAndUnit_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              MGGeneratorRDMUIcmdWithNucleusAndUnit.hh
//
// Version:             0.b.3
// Date:                29/02/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class Description
//
// The G4MGGeneratorRDMUIcmdWithNucleusAndUnit permits input of the nucleus definition
// in terms of its (atomic weight, atomic number, excitation energy).
// Input is expected in the form:
//
//			A, Z, E (energy unit)
//
// where A, Z, E are respectively the atomic weight, atomic number and
// excitation energy.  The energy unit can be any of the geant4 defined energy
// units.  The default is "keV"
//
// class description - end
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//
//
#include "G4UIcommand.hh"
#include "globals.hh"

#include "DSGeneratorRDMNucleus.hh"
////////////////////////////////////////////////////////////////////////////////
//
class DSGeneratorRDMUIcmdWithNucleusAndUnit : public G4UIcommand
{
public:  // with description
  DSGeneratorRDMUIcmdWithNucleusAndUnit
  (const char * theCommandPath,G4UImessenger * theMessenger);
  //    Constructor identifying the command path in the User Interface and the
  //    associated G4UImessenger which will use this G4UIcommand object.
  //
  ~DSGeneratorRDMUIcmdWithNucleusAndUnit();
  //    Desctructor
  //
  DSGeneratorRDMNucleus GetNewNucleusValue(G4String paramString);
  //    Extracts the values A, Z, E and unit from paramString.
  //
  G4double GetNewUnitValue(G4String paramString);
  //    Returns the value of the unit (paramString) as defined in geant4
  //
  G4String ConvertToString(DSGeneratorRDMNucleus nuc, const char * unit);
  //    Converts the Nucleus defined by nuc and the associated units of
  //    energy *unit into a G4String.
  void SetParameterName(const char * theNameA,const char * theNameZ,
                        const char * theNameE,
                        G4bool omittable, G4bool currentAsDefault=true);
  //    Identifies the parameter names associated with A, Z, E
  //
  void SetDefaultValue(DSGeneratorRDMNucleus defVal);
  //    Sets the default Nucleus if the command is invoked without any
  //    parameters.
  void SetUnitCandidates(const char * candidateList);
  //    Sets the list of unit candidates
  //
  void SetDefaultUnit(const char * defUnit);
  //    Sets the default unit if the command is invoked without any
  //    parameters.
};
////////////////////////////////////////////////////////////////////////////////
#endif



/*
 * Revision 1.1  2013/03/22 13:50:09  dfranco
 * Added radioactive decay and radioactive decay chain generators. They require the correspondent stacking actions. Two mac files are included as examples
 *
 */
