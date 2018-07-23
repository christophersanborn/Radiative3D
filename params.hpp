// params.hpp
//
// This file defines a number of "parameter" classes used to pass
// behavioral arguments in a simplified way.  Useful when a large
// number of parameters needs to be passed to a function or class
// constructor, or otherwise handed around.
//
// Currently just defines the MissionParams class, but in future I may
// (or may not) move the ModelParams or other similar classes over to
// here.
//
// As a general rule, it's the parameter classes that I want the user
// to be able to control via the command line that I develop here, so
// that only this file needs to be included by the module where
// command-line processing occurs (most likely main.cpp). If a
// parameter class is used purely internally, it should be defined
// alongside the class that it parameterizes, rather than here.
//
//
#ifndef PARAMS_H_
#define PARAMS_H_
//

//////
// CLASSES: Definitions
//
// INCLUDING:
//
//   o  class MissionParams
//

//////
// CLASS:   ::::  MissionParams  ::::
//
//   Basically a set of flags encoding the actions desired by the user
//   when they run the program.  Usually, they just want to run the
//   simulation, but sometimes they may want to run tests or diags.
//
class MissionParams {
public:

  bool bRunSim;                 // Run simulation
  bool bHelpMsg;                // Print help message and exit
  bool bDumpGrid;               // Output grid as plaintext
  bool bOutputModParamsOctv;    // Output Model params in GNU/Octave format
  bool bRTCoefTest;             // Perform R/T coefficient test
  bool bSourcePatternTest;      // Inspect source patterns

  std::string FNModParamsOctv;  // Filename for octave model params

  MissionParams() :
    bRunSim     (true),
    bHelpMsg    (false),
    bDumpGrid   (false),
    bOutputModParamsOctv (false),
    bRTCoefTest (false),
    bSourcePatternTest (false),
    FNModParamsOctv ("paramsout.octv")
  {}

  bool ModelBuildRequired() const {
    // Return true if any flags indicate the need for a constructed
    // Model object.
    return (bRunSim || bDumpGrid);
  }

};

///
#endif //#ifndef PARAMS_H_
//
