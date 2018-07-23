// cmdline.hpp
//
// This file develops the CmdOpt class.  This is used to interpret 
// commandline arguments, and to pass them user changable parameters.
//
// ~S. Walsh
//
//
#ifndef CMDLINE_H_
#define CMDLINE_H_
//
#include <iostream>
#include <vector>
#include <map>
#include "geom.hpp"

//////
// CLASS:  CmdOpt
//
//   A class to encapsulate command-line options in an easy-to-process
//   way.
//
class CmdOpt {
public:

  // ::::::::::::::::::::::::::::::::::::::
  // ::: Embedded Types  (CmdOpt Class) :::
  // ::::::::::::::::::::::::::::::::::::::

  enum cmdopts_e {      // Provides ID codes for each option type

    OPT_NOOP,           // For missing tokens
    OPT_UNK,            // Unknown tokens
    OPT_FREQ,           // Frequency
    OPT_NUMBER,         // Number of Phonons
    OPT_TTLIVE,         // Phonon Time-to-Live
    OPT_TOA,            // Take-off angle degree

    OPT_MODARGS,        // Model-compiled args
    OPT_COMPSELECT,     // Model-compiled selector index
    OPT_CYLRANGE,       // Cylinder range
    OPT_FLATTEN,        // Apply Earth-flattening xform to depth coords
    OPT_EARTHRAD,       // Set Earth radius for EFT and curved ECS's

    OPT_EVENT_MT,       // Event-source type / moment-tensor
    OPT_EVENT_LOC,      // Event-source location

    OPT_OVR_MFP,        // Override MFP's for Scatterers
    OPT_NODEFLECT,      // Suppress deflection on scattering

    OPT_REPORTS,        // Micro-reports on or off
    OPT_OUTDIR,         // Directory in which to write output files
    OPT_REPORT_FILE,    // Specify output file for reports
    OPT_OCSRAW,         // Use raw (untransformed) ICS coords for output
    OPT_SEISBINS,       // Number of time bins/cycle for seismometers
    OPT_SEISBINSIZE,    // Seismometer binwidth in seconds
    OPT_SEIS,           // Single seismometer
    OPT_SEISARRAY,      // Seismometer array
    OPT_SEIS_P2P,       // Seismometer Point-to-Point array
    OPT_SEIS_P2PW,      // Seismometer Point-to-Point array, wavelength units

                        // Symbols beginning with OPTM_ are "mission"
                        // options, and declare what actions we want
                        // the software to take.

    OPTM_HELP,          // Print help message and immediately exit
    OPTM_DUMPGRID,      // Print a plaintext dump of the grid
    OPTM_PARAMOUTFN,    // Output model params to filename
    OPTM_RTTEST,        // Perform R/T Coef test
    OPTM_EVENTTEST,     // Print analytics of event source
    OPTM_RUNSIM,        // Run simulation (on by default)

    OPT_NUMOPTS         // Dummy symbol for iterations. (Kindof a hack.)
                        // (Used in OutputAllRecognizedTokens)
  };

  typedef std::vector<CmdOpt> OptList;  // OptList: An array of CmdOpt
                                        // objects


protected:

  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Static Member Data  (CmdOpt Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::

  typedef std::map                // OpMap: A map type linking strings
          <Text,cmdopts_e>        // strings that might be passed on a
          OpMap;                  // command line to the OP_codes that
                                  // represent them.

  static OpMap opMap_;            // An OpMap. Used to allow CmdOpt
                                  // objects to self-identify their
                                  // OP_code IDs.


  // ::::::::::::::::::::::::::::::::::::::
  // ::: Static Methods  (CmdOpt Class) :::
  // ::::::::::::::::::::::::::::::::::::::

  inline static OpMap init_map();   
            // (Defined in this file after class
            // block.) This function populates the
            // OP_code map.  Called automatically by
            // initialization of member opMap_.

public:

  static OptList PackageArgCArgV(int, char **);
            // When given an argc and argv representation of command
            // line args, this function return a vector of CmdOpt
            // objects representing those same args.

  static void OutputAllRecognizedTokens(std::ostream * out);
            // Prints all recognized tokens as a comma separated
            // list. Used in --help output.


protected:
  ;
  // :::::::::::::::::::::::::::::::::::::::::::
  // ::: Private Member Data  (CmdOpt Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::

  cmdopts_e mOpID;          // Input argument mapped

  Text mToken;              // Input Argument
  Text mValue;              // Input Value (if applicable)
  Text mOrigValue;          // Unaltered Input Value (does not change
                            // with "pop" operations)

public:

  // ::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (CmdOpt Class) :::
  // ::::::::::::::::::::::::::::::::::::

  CmdOpt(std::vector<Text> &arglist);


public:

  // ::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (CmdOpt Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::

  int GetID() const {return mOpID;}
  Text GetTokenText() const {return mToken;}
  Text GetValueText() const {return mValue;}
  Text GetOptionText() const {return mToken + " " + mOrigValue;}
  Count GetValueCount() const;  // Assuming mValue to contain a comma
                                // separated list, returns the number
                                // of values in said list.


  // :::::::::::::::::::::::::::::::
  // ::: Public Member Functions :::
  // :::::::::::::::::::::::::::::::

  void Display() const
  {
    std::cout << mToken << "       " 
              << mValue << "       " 
              << mOpID <<  std::endl;
  }

  // :::::::::::::::::::::::::::::::::::
  // ::: Number Validation Functions :::
  // :::::::::::::::::::::::::::::::::::

  static bool ValidFixedPointFloat(Text numtxt);
                // Verifies number is a valid float
  
  static bool ValidInteger(Text numtxt);
                // Verifies number is an integer

  static bool ValidScientific(Text numtxt);
                // Verifies number is in proper scientific notation

  static bool ValidTokenFormat(Text toktxt);
                // True if string value looks like a token (as opposed
                // to a (possibly negative) number, or a filename, or
                // something not a token).


  // :::::::::::::::::::::::::::
  // ::: Pop Value Functions :::
  // :::::::::::::::::::::::::::
  //
  //    These functions remove and process a subset of the string data
  //    contained in mValue.  They assume that mValue is a comma-
  //    separated list of stringified data.  They remove up to the
  //    first comma (if any), and return the value interpreted as the
  //    desired numeric or other type. mValue is then left containing
  //    everything to the right of the comma.
  //

  bool PeekValue() const;       // True if mValue not empty. Used to test
                                // presence of value(s) that can be Peeked/
                                // Popped, without regard to value type.

  int PeekValue_Integer(bool required=true, int defval=0) const;
  Real PeekValue_Real(bool required=true, Real defval=0) const;
  Text PeekValue_Text(bool required=true, Text defval="") const;
  //        PeekValue_ functions: Return a "value" from the leftmost chunk
  //        of mValue (up to a comma separator) interpreted as the requested
  //        type. If the mValue chunk cannot be interpretted as the
  //        particular type, then an exception is thrown.  If the chunk is
  //        empty, then return the provided default value (if
  //        required==false) or else throw an exception.


  int PopValue_Integer(int defval) {return PopValue_Integer(false,defval);}
  Real PopValue_Real(Real defval)  {return PopValue_Real(false,defval);}
  Text PopValue_Text(Text defval)  {return PopValue_Text(false,defval);}
  int PopValue_Integer(bool required=true, int defval=0);
  Real PopValue_Real(bool required=true, Real defval=0);
  Text PopValue_Text(bool required=true, Text defval="");
  R3::XYZ PopValue_XYZ();
  //        The PopValue_ functions return values from the leftmost chunk of
  //        mValue, and then gobble mValue up through the chunk seperator
  //        (comma).  If called with no args, then the chunk must contain a
  //        value, or else an exception occurs.  If called with single arg,
  //        then it is the default value and user-provided value is
  //        optional.  In either case, an exception is thrown if the chunk
  //        is non-empty but not valid for the requested type.


  bool PeekValue_IsInteger() const;
  bool PeekValue_IsReal() const;
  //        These return true if the leftmost mValue chunk contains text
  //        that can be validly interpretted as the requested type. Returns
  //        false if the chunk is either empty or not valid for the
  //        requested type.  These function do NOT throw exceptions.


  Text GetValueBeforeSep() const;         // Returns text before separator
  void KillValueThroughSep();             // Shortens mValue by
                                          // killing up through first
                                          // separator


};// class CmdOpt
///


//////
// METHOD:   CmdOpt :: init_map()   (static, inline)
//
//   Populates an OpMap object and returns it to the caller. Used to
//   initialize the static opMap_ member.
//
CmdOpt::OpMap CmdOpt::init_map() {

  OpMap map;

  map[""] = OPT_NOOP;                               

  // Universe:
  map["-F"] =                      OPT_FREQ;  // Frequency
  map["--frequency"] =             OPT_FREQ;  //  ''
  map["-N"] =                    OPT_NUMBER;  // Number of Phonons
  map["--num-phonons"] =         OPT_NUMBER;  //  ''
  map["-T"] =                    OPT_TTLIVE;  // Phonon Time to Live
  map["--timetolive"] =          OPT_TTLIVE;  //  ''
  map["-A"] =                       OPT_TOA;  // Take-off Angle degree
  map["--toa-degree"] =             OPT_TOA;  //  ''

  // Model:
  map["--grid-compiled"] =   OPT_COMPSELECT;  // Compiled model select
  map["--model-args"] =         OPT_MODARGS;  // Model args value list
  map["--model-compiled-args"]= OPT_MODARGS;  //  '' (deprecated)
  map["--range"] =             OPT_CYLRANGE;  // Cylinder Range
  map["--cylinder-range"] =    OPT_CYLRANGE;  //  '' (deprecated)
  map["--flatten"] =            OPT_FLATTEN;  // Earth-flattening transform
  map["--earthrad"] =          OPT_EARTHRAD;  // Earth Radius
  map["--earthradius"] =       OPT_EARTHRAD;  //  ''

  // Event Source:
  map["-E"] =                  OPT_EVENT_MT;  // Src type (Moment Tens)
  map["--source"] =            OPT_EVENT_MT;  //  ''
  map["-L"] =                 OPT_EVENT_LOC;  // Source location
  map["--source-loc"] =       OPT_EVENT_LOC;  //  ''

  // Sim Parameters:
  map["--mfpoverride"] =        OPT_OVR_MFP;  // Override MFPs
  map["--overridemfp"] =        OPT_OVR_MFP;  //  ''
  map["--nodeflect"] =        OPT_NODEFLECT;  // No scattering deflection
  map["--no-deflect"] =       OPT_NODEFLECT;  //  ''

  // Reporting:
  map["--reports"] =            OPT_REPORTS;  // Reports
  map["--report-file"] =    OPT_REPORT_FILE;  // Filename for reports
  map["--output-dir"] =          OPT_OUTDIR;  // Directory for output
  map["--ocsnotransform"] =      OPT_OCSRAW;  // Output Coords no-transform
  map["--ocsraw"] =              OPT_OCSRAW;  //  ''
  map["--binspercycle"] =      OPT_SEISBINS;  // Bins per Cycle
  map["--bins"] =              OPT_SEISBINS;  //  ''  [DEPRECATED FORM]
  map["--binsize"] =        OPT_SEISBINSIZE;  // Bin size in seconds
  map["--seismometer"] =           OPT_SEIS;  // Single seismometer
  map["--seis-array"] =       OPT_SEISARRAY;  // Seis-array
  map["--seis-p2p"] =          OPT_SEIS_P2P;  // Point-to-Point Array
  map["--seis-p2pw"] =        OPT_SEIS_P2PW;  // Point-to-Point Array

  // Mission:
  map["--help"] =                 OPTM_HELP;  // Print help message
  map["--dump-grid"] =        OPTM_DUMPGRID;  // Plaintext grid dump
  map["--mparams-outfile"] =OPTM_PARAMOUTFN;  // Params output filename
  map["--rtcoef-test"] =        OPTM_RTTEST;  // R/T coef test
  map["--event-test"] =      OPTM_EVENTTEST;  // Analyze event source
  map["--run-simulation"] =     OPTM_RUNSIM;  // Run simulation (default)
  map["--run-sim"] =            OPTM_RUNSIM;  //  ''
  
 
  return map;

}




///
#endif //#ifndef CMDLINE_H_
//
