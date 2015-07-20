#include "DSGeneratorRDMNucleus.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
////////////////////////////////////////////////////////////////////////////////
//
DSGeneratorRDMNucleus::DSGeneratorRDMNucleus ()
  : a(24), z(11), e(0.0)   //Default to 24Na radioactivity
{;}
///////////////////////////////////////////////////////////////////////////////
//
DSGeneratorRDMNucleus::DSGeneratorRDMNucleus (G4int a1, G4int z1, G4double e1)
{
  //
  //
  a = a1;
  z = z1;
  e = e1;
}
///////////////////////////////////////////////////////////////////////////////
//
DSGeneratorRDMNucleus::~DSGeneratorRDMNucleus ()
{;}
///////////////////////////////////////////////////////////////////////////////
//
std::ostream &operator << (std::ostream &mys, const DSGeneratorRDMNucleus &q)
//
//
// Definition of the insertion operator << to provide the nucleus limits to
// ostream.
//
{
  mys <<"Atomic weight: " <<q.GetA()
    <<"Atomic number: " <<q.GetZ()
    <<"Excitation energy: "<<q.GetE();
  return mys;
}
///////////////////////////////////////////////////////////////////////////////






/*
 * Revision 1.1  2013/03/22 13:50:09  dfranco
 * Added radioactive decay and radioactive decay chain generators. They require the correspondent stacking actions. Two mac files are included as examples
 *
 */


