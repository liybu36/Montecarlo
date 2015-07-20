// --------------------------------------------------------------------------//
/** 
 * AUTHOR: Davide Franco
 * 
 * Taken from  MaGe code which was originally taken from the Babar code 
 * 27-7-2005
 */
// --------------------------------------------------------------------------//


//---------------------------------------------------------------------------//

#include "DSLogger.hh"      //Present DS Class Headers 
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSIO.hh"      //Present DS Class Headers 

//---------------------------------------------------------------------------//


#include <sstream>
#include <cstdlib>

using namespace std;

static ostringstream devnull;

DSLogger::Severity DSLogger::_minSeverity = DSLogger::debugging;
//ostream * DSLogger::_myOstream (&cout);
ostream * DSLogger::_myOstream  = 0 ;
ostream * DSLogger::_myErrstream (&cerr);
ostream * DSLogger::_myNullstream (&devnull);
bool DSLogger::_doPrint = true;

DSLogger::DSLogger(){}

DSLogger::~DSLogger(){}

DSLogger::DSLogger(const DSLogger & other){}

ostream&  DSLogger::msg(DSLogger::Severity severity, 
		       const char* facility )
{
  _doPrint = true;
  if(!_myOstream) _myOstream = new ostringstream ;
  else ((ostringstream*) _myOstream)->str("");
  
  if(severity >= _minSeverity){
    *_myOstream << toString(severity) << ":" << facility << ":";
  }else{
    _doPrint =false;
    return *_myNullstream ;
  }

  if ( severity == fatal ){
    *_myOstream << ::endlog;
    ::abort();
  }
  return *_myOstream;
}

void  DSLogger::endlog(ostream& os){
  if(_doPrint) {
    cout << ((ostringstream*) _myOstream)->str() << endl;
    DSIO::Get()->GetStreamLogFile() << ((ostringstream*) _myOstream)->str() << endl;
  
  }
//  if(_doPrint)*_myOstream << endl;
}

const char* DSLogger::toString(DSLogger::Severity sever){
  switch (sever) {
  case -2:
    return "Debug";
    break;
  case -1:
    return "Develop";
    break;
  case 0:
    return "Trace";
    break;
  case 1:
    return "Routine";
    break;
  case 2:
    return "Warning";
    break;
  case 3:
    return "Error";
    break;
  case 4:
    return "Fatal";
    break;
  }
  return "Fatal";
}

DSLogger::Severity DSLogger::toEnum(const  std::string& level){
  if(level == "development") return development ;
  if(level == "debugging")   return debugging ;
  if(level == "trace")       return trace;
  if(level == "routine")     return routine;
  if(level == "fatal")       return fatal;
  if(level == "warning")     return warning;
  if(level == "error")       return error;
  return fatal;
}

ostream& endlog(ostream& os){
  DSLogger::endlog(os);
  return os;
}
/*
 * $Log: DSLogger.cc,v $
 * Revision 1.1  2014/05/07 12:21:04  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.3  2013/04/12 14:36:44  cocco
 * Modified script in order to setup correct environment variables on LNGS cluster
 * Added include cstdlib to handle ::abort function in DSLogger.cc
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
