// cmdline.cpp
//
#include <cstdlib>
#include "cmdline.hpp"


//////
// *** STATIC MEMBERS:  CmdOpt Class
//

CmdOpt::OpMap CmdOpt::opMap_ = init_map();


//////
// *** STATIC METHODS:  CmdOpt Class
//


//////
// METHOD:  CmdOpt :: PackageArgCArgV()         (static)
//
//   Static method that returns an OptList (a vector of CmdOpt
//   objects) representing the parsing of argv[] into token,value
//   pairs.
//
CmdOpt::OptList CmdOpt::PackageArgCArgV(int argc, char * argv[]) {

  std::vector<Text> arglist;
  for (int i=1; i<argc; i++) {          // Reformat argv[]
    arglist.push_back(argv[i]);         // into vector<string> 
  }                                     //
  OptList ret;
  while (!arglist.empty()) {            // Gobble the vector<string>
    ret.push_back(CmdOpt(arglist));     // into new CmdOpt objects
  }                                     //
  if (false) {
    for (Index i = 0; i < ret.size(); i++) {
      ret[i].Display();                 // Debug stuff
    }                                   //
  }
  return ret;

}

//////
// METHOD:  CmdOpt :: OutputAllRecognizedTokens()   (static)
//
//   Outputs to a char ostream a list of all defined tokens in the
//   OpMap map.  Used to list valid options in the help routine in
//   main.cpp.  (A very sparse online help...)
//
void CmdOpt::OutputAllRecognizedTokens(std::ostream * out) {

  int k = OPT_UNK+1;
  for (; k<OPT_NUMOPTS; k++) {
    OpMap::const_reverse_iterator it = opMap_.rbegin();
    Count found = 0;
    for (; it!=opMap_.rend(); it++) {
      if (it->second == k) {
        if (found++ > 0) {*out << ", ";} 
        else             {*out << "  ";}
        *out << it->first;
      }
    }
    if (found > 0) {
      *out << "\n";
    }
  }
}


/////
// CONSTRUCTOR:  CmdOpt()       ( NOTICE: ** MODIFIES INPUT ** )
//
//   Parses command line arguments into a parameter literal (mToken)
//   and its corresponding value (mValue) if applicable, and then
//   matches mToken to an OP_ID code useful in switch-case clauses.
//
//   Performs some VERY first-level validation of mValue
//   contents. (Namely, strips surrounding whitespace and trailing
//   commas (but not leading commas).)
//
//   MODIFIES incoming 'arglist' argument by gobbling the elements
//   used to construct the CmdOpt object.  This allows the calling
//   routine to contine creating yet more CmdOpt objects from the
//   remaining elements of 'arglist'
//
CmdOpt::CmdOpt(std::vector<Text> &arglist)
{
  Text FirstArg = arglist.at(0);
  Text SecondArg = "";
  if (arglist.size() > 1) {
    SecondArg = arglist.at(1);
  }

  if (!ValidTokenFormat(FirstArg)) {    // If NOT a dash-token:
    mToken = "";                        // Then retain arg as a VALUE of
    mValue = FirstArg;                  // a NOOP token.
    mOpID = OPT_NOOP;
    arglist.erase (arglist.begin());
  }
  else {    // Else we are a properly formatted dash token:
            // (could be single or double dash)

    if (FirstArg.at(1) == '-')  { // *** Check for "long" (double dash) option
                                  // *** (ie, "--frequency")

      if (FirstArg.find('=') != Text::npos) {     // Long options with equals,
        std::size_t eqpos = FirstArg.find('=');   // ie --frequency=2.6
        mToken = FirstArg.substr(0,eqpos);
        mValue = FirstArg.substr(eqpos+1,FirstArg.size()-1);
        arglist.erase (arglist.begin());
      }
      else if ((arglist.size() > 1) && !ValidTokenFormat(SecondArg)) {
        mToken = FirstArg;                        // Long options w/o equals,
        mValue = SecondArg;                       // ie --frequency 2.6
        arglist.erase (arglist.begin(),arglist.begin()+2);
      }
      else {                    // Long options which do not require a value, 
        mToken = FirstArg;      // i.e., --something --something-else
        mValue = "";
        arglist.erase (arglist.begin());
      }

    }
    else {    // *** Else we are "short" (single dash) option(s)
              // *** Process these:
      
      if (FirstArg.size() == 2) {       // Checks for short singletons, ie -F
        mToken = FirstArg;
        if ((arglist.size() > 1) && !ValidTokenFormat(SecondArg)) {
          mValue = SecondArg;           // Short options with value, ie -F 2.6
          arglist.erase (arglist.begin(),arglist.begin()+2);
        } else {
          mValue = "";                  // Short options without value
          arglist.erase (arglist.begin());
        }
      } else { 
        mToken = FirstArg.substr(0,2);  // Combos of short options, ie -NFG
        mValue = "";
        FirstArg.erase (1,1);           // Kill first character after dash
        arglist.at(0) = FirstArg;       // but leave remainder for next
      }                                 // iteration to catch.
    }
  } //
  ///   We now have an mToken and an mValue
  //
  // OK, Now map command line argument (mToken) to a corresponding
  // OP_ID code enumerated in cmdline.hpp:
  //
  OpMap::iterator FoundOp;
  FoundOp = opMap_.find(mToken);        // Initial search  
  if (FoundOp != opMap_.end()) {        // Test if found match
    mOpID = opMap_[mToken];             // Use match if found,
  } else {                              // else use special OP_ID for
    mOpID = OPT_UNK;                    // unknown token.
  }
  //
  // Now a VERY basic validation of mValue.  The following code
  // guarantees that mValue has (1) no surrounding whitespace, and (2)
  // no trailing commas.  NOTE: Only 0x20 spaces are counted as whitespace.
  //
  mOrigValue = mValue;
  int prvlen = 0;
  int curlen = mValue.size();
  while ((curlen > 0)&&(curlen != prvlen)) {  // Whittle till stable or empty
    if(mValue.at(0)==' ') {
      mValue.erase(0,1);
    } else if (mValue.at(curlen-1)==' ') {
      mValue.erase(curlen-1,1);
    } else if (mValue.at(curlen-1)==',') {    // Note: trailing ,'s will be
      mValue.erase(curlen-1,1);               // trimmed, but leading ,'s
    }                                         // will not.
    prvlen=curlen;
    curlen = mValue.size();
  }

} // END CmdOpt::CmdOpt()
///


//////
// METHOD:  CmdOpt :: GetValueCount()
//
//   Assuming mValue to be a comma separated list, returns number of
//   values in said list.  Basically, just counts commas.  Perhaps we
//   want to get more sophisticated in the future, but this is how
//   it's done now.
//
Count CmdOpt::GetValueCount() const {
  if (mValue.length() == 0) return 0;
  Count count = 0;  // else count commas
  for (std::size_t offset = mValue.find(',');
       offset != Text::npos;
       offset = mValue.find(',', offset+1)) {
    ++count;
  }
  return (count+1);
}


//////
// METHOD:  CmdOpt :: ValidTokenFormat()        (static)
//
//   True if string looks like a token (ie "-N" or "--token")
//
bool CmdOpt::ValidTokenFormat(Text toktxt) {

  if (toktxt.length() < 2) {    // Token MUST contain at least two
    return false;               // characters (eg, can't be shorter
  }                             // than "-N", for example).

  if (toktxt.at(0) != '-') {    // First character of a token MUST
    return false;               // be a "-".
  }                             //

  if (toktxt.at(1) != '-') {    // (single-dash option)

    if (isdigit(toktxt.at(1))) {        // Token name cannot be or start
      return false;                     // with a digit (creates ambiguity
    }                                   // with negative numbers)

    // fall-through

  } 
  else {                        // (double-dash option)

    if ((toktxt.length() > 2) && (isdigit(toktxt.at(2)))) {
      return false;                     // Again, token name cannot be
    }                                   // or begin with a digit.

    // fall-through; Note that the blank long-option "--" DOES pass
    // this test, which is by design.  (Some programs use "--" to
    // indicate end-of-options, eg to follow with filenames which
    // might *look* like options, eg a file named "-N".  We don't use
    // that feature here, but I am hereby reserving the option of
    // using it.

  }

  return true;  // If we get to here, there were no disqualifiers, and
                // the text LOOKS like it's an option token.  Return
                // true.

}//
//


//////
// METHOD:  CmdOpt :: ValidFixedPointFloat()    (static)
//
//   Verifies string represents a valid fixed-point float.
//
bool CmdOpt::ValidFixedPointFloat(Text numtxt) {
  if (numtxt.length()==0) {
    return false;
  }
  int dotcount = 0;
  if (numtxt.at(0) == '.') {
    dotcount++;
  }
  if (isdigit(numtxt.at(0)) || numtxt.at(0)=='+'
                            || numtxt.at(0)=='-' || numtxt.at(0)=='.') {
    for (Index i=1; i < numtxt.size(); i++) {
      if (isdigit(numtxt.at(i))||numtxt.at(i)==('.')) {
        if (numtxt.at(i) == '.') {dotcount++;}
      }
      else {return false;}
    }
    if (dotcount > 1) {return false;}
    else {return true;}
  }
  else {
    return false;
  }
}//
//


//////
// METHOD:  CmdOpt :: ValidInteger()            (static)
//
bool CmdOpt::ValidInteger(Text numtxt) {
  if (ValidFixedPointFloat(numtxt) == true) {
    int dotcount = 0;
    for (Index i=0; i < numtxt.size(); i++) {
      if (numtxt.at(i) == '.') {dotcount++;}
    }
    if (dotcount == 0) {return true;}
    else {return false;}
  }
  else {
    return false;
  }
}//
//


//////
// METHOD:   CmdOpt :: ValidScientific()                    (static)
//
//   Returns true if text denotes a number in scientific notation,
//   with or without an exponent. (I.e., "3.14" and "3.14e0" both
//   return true.)
//
bool CmdOpt::ValidScientific(Text numtxt) {
  int ecount = 0;
  int epos = 0;
  for (Index i=0; i < numtxt.size(); i++) {
    if (numtxt.at(i) == 'e' || numtxt.at(i) == 'E') {
      epos = i;
      ecount++;
    }
  }
  if (ecount == 1) {    // One (and only one) 'E': 
    Text substr1 = numtxt.substr(0,epos);
    Text substr2 = numtxt.substr(epos+1);
    return (ValidFixedPointFloat(substr1)
                 && ValidInteger(substr2));
  }
  else if (ecount == 0) { // No 'E's - check for fixed point:
    return ValidFixedPointFloat(numtxt);
  }
  else {  // There are too many 'E's 
    return false;
  }
}


//////
// METHOD:  CmdOpt :: PeekValue()
//
bool CmdOpt::PeekValue() const {
  return (!mValue.empty());
}


//////
// METHOD:  CmdOpt :: PopValue_Integer()
// METHOD:  CmdOpt :: PeekValue_IsInteger()
// METHOD:  CmdOpt :: PeekValue_Integer()
//
int CmdOpt::PopValue_Integer(bool required, int retval) {
  retval = PeekValue_Integer(required, retval);
  KillValueThroughSep();
  return retval;
}//
//
bool CmdOpt::PeekValue_IsInteger() const {
  try {PeekValue_Integer(true,0);}
  catch(...) {return false;}
  return true;
}//
//
int CmdOpt::PeekValue_Integer(bool required,  // Should empty mValue throw?
                              int  defval     // Default return val if empty
                             ) const {        //
  Text peektxt = GetValueBeforeSep();
  if (peektxt.empty()) {                // Check for empty mValue,
    if (required) {                     //   ...and throw if appropriate,
      throw(Runtime("Required value not provided;"
                    + Text(" expected integer value.")));}
  }                                     //   ...else GOTO endif
  else { // (peektext not empty)        // Otherwise try decoding as int:
    char suffix = *peektxt.rbegin();    // Get last char, possible multiplier
    peektxt.resize(peektxt.size()-1);   // Trim last char from peektxt
    switch (suffix) {
    case 'K':
      peektxt.append("000");
      break;
    case 'M':
      peektxt.append("000000");
      break;
    case 'B':
      peektxt.append("000000000");
      break;
    default:
      peektxt.append(1,suffix);         // Restore last char if not a suffix
      break;
    }
    if (ValidInteger(peektxt)) {        // Check validity of peektxt as int
      defval = atoi(peektxt.c_str());   // Use as return val if valid
    } else {                            // Else throw exception.
      throw(Runtime("Invalid value: cannot interpret '"
                    + peektxt + "' as type Integer."));
    }
  } // endif
  return defval;
}//
//


//////
// METHOD:  CmdOpt :: PopValue_Real()
// METHOD:  CmdOpt :: PeekValue_Real()
// METHOD:  CmdOpt :: PeekValue_IsReal()
//
Real CmdOpt::PopValue_Real(bool required, Real retval) {
  retval = PeekValue_Real(required, retval);
  KillValueThroughSep();
  return retval;
}//
//
bool CmdOpt::PeekValue_IsReal() const {
  try {PeekValue_Real(true,0);}
  catch(...) {return false;}
  return true;    
}//
//
Real CmdOpt::PeekValue_Real(bool required, Real defval) const {
  Text peektxt = GetValueBeforeSep();
  if (peektxt.empty()) {                // Check for empty mValue,
    if (required) {                     //   ...and throw if appropriate,
      throw(Runtime(                    //
      "Required value not provided; expected numeric value."));}
  }                                     //   ...else GOTO endif
  else { // (peektext not empty)        // Otherwise try decoding as Real:

    if (ValidScientific(peektxt)) {
      defval = strtod(peektxt.c_str(),NULL);
    } else if (peektxt=="INF" || peektxt=="inf") {
      defval = 1.0/0.0;
    } else if (peektxt=="-INF" || peektxt=="-inf") {
      defval = -1.0/0.0;
    }
    else {
      throw(Runtime("Invalid value: cannot interpret '"
                    +peektxt+"' as type Real."));
    }
  } // endif
  return defval;
}//
//


//////
// METHOD:  CmdOpt :: PopValue_XYZ()
//
R3::XYZ CmdOpt::PopValue_XYZ() {
  Real x = PopValue_Real();
  Real y = PopValue_Real();
  Real z = PopValue_Real();
  return R3::XYZ(x,y,z);
}//
//


//////
// METHOD: CmdOpt :: PopValue_Text()
// METHOD: CmdOpt :: PeekValue_Text()
//
Text CmdOpt::PopValue_Text(bool required, Text retval) {
  retval = PeekValue_Text(required, retval);
  KillValueThroughSep();
  return retval;
}//
//
Text CmdOpt::PeekValue_Text(bool required, Text defval) const {
  Text peektxt = GetValueBeforeSep();
  if (peektxt.empty()) {                // Check for empty mValue,
    if (required)                       //   ...and throw if appropriate,
      {throw(Runtime(                   //
      "Required value not provided; expected text value."));}
  }                                     //   ...else GOTO endif & return defval
  else { // (peektext not empty)        // Otherwise,
    defval = peektxt;                   // use peektxt for return value
  }
  return defval;
}//
//


//////
// METHOD:  CmdOpt :: GetValueBeforeSep()
// METHOD:  CmdOpt :: KillValueThroughSep()
//
Text CmdOpt::GetValueBeforeSep() const {
  std::size_t commapos = mValue.find_first_of(',');
  Text txt = mValue.substr(0,commapos);
  return txt;
}//
//////
void CmdOpt::KillValueThroughSep() {
  std::size_t commapos = mValue.find_first_of(',');
  if (commapos==Text::npos) {
    mValue = "";
  } else {
    mValue = mValue.substr(commapos+1,Text::npos);
  }
}//
//
