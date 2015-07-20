#include "G4UIcommand.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "DSGeneratorRDMUIcmdWithNucleusAndUnit.hh"

#include <sstream>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
DSGeneratorRDMUIcmdWithNucleusAndUnit::DSGeneratorRDMUIcmdWithNucleusAndUnit
(const char * theCommandPath,G4UImessenger * theMessenger)
:G4UIcommand(theCommandPath,theMessenger)
{
  G4UIparameter * intParamA = new G4UIparameter('i');
  SetParameter(intParamA);
  G4UIparameter * intParamZ = new G4UIparameter('i');
  SetParameter(intParamZ);
  G4UIparameter * dblParamE = new G4UIparameter('d');
  SetParameter(dblParamE);
  G4UIparameter * untParam = new G4UIparameter('s');
  SetParameter(untParam);
  untParam->SetParameterName("Unit");
}
////////////////////////////////////////////////////////////////////////////////
//
DSGeneratorRDMUIcmdWithNucleusAndUnit::~DSGeneratorRDMUIcmdWithNucleusAndUnit()
{
  ;
}
////////////////////////////////////////////////////////////////////////////////
//
DSGeneratorRDMNucleus DSGeneratorRDMUIcmdWithNucleusAndUnit::GetNewNucleusValue(G4String paramString)
{
  G4int a;
  G4int z;
  G4double e;
  char unts[30];

  istringstream is(paramString);
  is >> a >> z >> e >>unts;
  G4String unt = unts;

  return DSGeneratorRDMNucleus(a,z,e*ValueOf(unt));
}

G4double DSGeneratorRDMUIcmdWithNucleusAndUnit::GetNewUnitValue(G4String paramString)
{
  G4int a;
  G4int z;
  G4double e;  

  char unts[30];
  
  istringstream is(paramString);
  is >> a >> z >> e  >> unts;

  G4String unt = unts;
  
  return ValueOf(unt);
}

////////////////////////////////////////////////////////////////////////////////
//
G4String DSGeneratorRDMUIcmdWithNucleusAndUnit::ConvertToString(DSGeneratorRDMNucleus def, 
						    const char *unitName)
{
  G4double uv = ValueOf(unitName);

  ostringstream os;
  os << def.GetA() << " " << def.GetZ()
     << " "<< def.GetE()/uv<<" "<< unitName <<  '\0';
  return G4String(os.str());
}                         
////////////////////////////////////////////////////////////////////////////////
//
void DSGeneratorRDMUIcmdWithNucleusAndUnit::SetParameterName
(const char * theNameA, const char * theNameZ,
const char * theNameE,G4bool omittable,G4bool currentAsDefault)
{
  G4UIparameter * theParamA = GetParameter(0);
  theParamA->SetParameterName(theNameA);
  theParamA->SetOmittable(omittable);
  theParamA->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamZ = GetParameter(1);
  theParamZ->SetParameterName(theNameZ);
  theParamZ->SetOmittable(omittable);
  theParamZ->SetCurrentAsDefault(currentAsDefault);
  G4UIparameter * theParamE = GetParameter(2);
  theParamE->SetParameterName(theNameE);
  theParamE->SetOmittable(omittable);
  theParamE->SetCurrentAsDefault(currentAsDefault);
}
////////////////////////////////////////////////////////////////////////////////
//
void DSGeneratorRDMUIcmdWithNucleusAndUnit::SetDefaultValue(DSGeneratorRDMNucleus def)
{
  G4UIparameter * theParamA = GetParameter(0);
  theParamA->SetDefaultValue(def.GetA());
  G4UIparameter * theParamZ = GetParameter(1);
  theParamZ->SetDefaultValue(def.GetZ());
  G4UIparameter * theParamE = GetParameter(2);
  theParamE->SetDefaultValue(def.GetE());
}


void DSGeneratorRDMUIcmdWithNucleusAndUnit::SetUnitCandidates(const char * candidateList)
{
  G4UIparameter * untParam = GetParameter(3);
  G4String canList = candidateList;
  untParam->SetParameterCandidates(canList);
}

void DSGeneratorRDMUIcmdWithNucleusAndUnit::SetDefaultUnit(const char * defUnit)
{
  G4UIparameter * untParam = GetParameter(3);
  untParam->SetOmittable(true);
  untParam->SetDefaultValue(defUnit);
}





/*
 * Revision 1.1  2013/03/22 13:50:09  dfranco
 * Added radioactive decay and radioactive decay chain generators. They require the correspondent stacking actions. Two mac files are included as examples
 *
 */









