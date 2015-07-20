// --------------------------------------------------------------------------//
/** 
 * AUTHOR: Davide Franco
 * 
 * Taken from the MaGe code which was originally taken from the Babar code 
 * 27-7-2005
 */
// --------------------------------------------------------------------------//

#ifndef _BXLOGGER_HH
#define _BXLOGGER_HH

#include <iostream>
#include <string>

#define ERRLINE_HACK_1(line)   #line
#define ERRLINE_HACK_2(line)   ERRLINE_HACK_1(line)

#ifdef DSLog
#undef DSLog
#endif
#define DSLog(sev) DSLogger::msg( DSLogger::sev, __FILE__ "(" ERRLINE_HACK_2(__LINE__) ")" )

//#define DSLog(sev) DSLogger::msg( DSLogger::sev, __FILE__ "(" ERRLINE_HACK_2(__LINE__) ")", 0 )

//---------------------------------------------------------------------------//

class DSLogger 
{
public:

  // a severity enum is defined; use only these values
  enum Severity {debugging=-2, development=-1, trace=0, routine, warning, error, fatal};
   //	fatal:		The message is related to a condition preventing
  //			further execution of the program.  ErrLogger will
  //			terminate the program.  Programmers should not call
  //			abort or exit themselves.

  //	error:		A condition exists such that requested result
  //			or action can not be produced.  This is a serious

  //	warning:	The result is produced, but may not be
  //			what's desired due to an unexpected condition

  //	routine:	Nothing known to be wrong with the result;
  //			messages that are always produced in normal
  //			operation

  //	trace:		Messages about the flow of program control
  //			and which optional operations took place.
  //			(This is the default if nothing is defined)
  
  //	debugging:	Information in addition to the above

 // members
  static Severity GetSeverity() { return _minSeverity; }
  static std::ostream& msg( DSLogger::Severity severity, 
		       const char* facility);//, const char* code );

  static void SetSeverity(Severity sever){ _minSeverity = sever;}
  static void endlog(std::ostream& s);
//  static void endlog();
  static Severity toEnum(const std::string& );

protected:
 //default constructor
  DSLogger();

  //copy constructor
  DSLogger(const DSLogger &);

  //destructor
  ~DSLogger();

  //public interface

  //protected members
protected:


  //private  members
private:
  
  static const char* toString(Severity);

  static std::ostream* _myOstream;
  static std::ostream* _myErrstream;
  static std::ostream* _myNullstream;

  static Severity _minSeverity;

  static bool _doPrint;
};
std::ostream& endlog(std::ostream& s);

#endif
/*
 * $Log: DSLogger.hh,v $
 * Revision 1.1  2014/05/07 12:20:53  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
