#include <iostream>
#include <cstdlib>      /* for exit() */
#include <cmath>
#include <sstream>
#include <string>
#include <stdexcept>
#include "params.hpp"
#include "model.hpp"
#include "dataout.hpp"  /* for SuppressAllReports() */
#include "rtcoef.hpp"
#include "cmdline.hpp"
//
#ifndef REVISION_NUM
#define REVISION_NUM "N/A"
#endif

void print_banner();                    // Prints startup banner
void print_help();                      // Prints help message
void process_option(CmdOpt &, 
                    ModelParams &,      // Process a command-line
                    MissionParams &);   // option


//////
// FUNCTION:  main()
//
int main(int argc, char *argv[]) {

  print_banner();

  MissionParams mission;
  ModelParams MParams;

  dataout.SuppressAllReports();   // Don't give the play-by-play; just
                                  // dump the seismic traces at the
                                  // end.

  // ***
  // ** Process Command-Line Options:
  // *

  CmdOpt::OptList opt_list                  // Get "options" from
    = CmdOpt::PackageArgCArgV(argc, argv);  // "arguments"

  for (Index i = 0;                         //
       i < opt_list.size(); i++) {          // Process all the Options
                                            //
    try {process_option(opt_list[i], MParams, mission);}
    catch (std::exception &e) {
      std::cout << "** Error processing command-line option: "
                << opt_list[i].GetOptionText() << "\n** Message: "
                << e.what() << "\n** Exiting...\n";
      exit(1);
    }
  } // Done processing command-line args

  // ***
  // ** Now DO the things the user asked for:
  // *

  if (mission.bHelpMsg) {     // Print help message and exit
    print_help();             //
    exit(0);
  }

  if (true) {                 // Dump Model Parameters
    MParams.Output();         //
  }

  if (mission.bOutputModParamsOctv) {   // Output model params for Octave
    std::stringstream outfilename;
    outfilename << dataout.GetOutputDirectory();
    if (outfilename.str().size() > 0) {
      outfilename << "/";   // directory separator, if needed
    }
    outfilename << mission.FNModParamsOctv;
    std::ofstream outfile(outfilename.str().c_str());
    MParams.OutputOctaveText(&outfile);
    outfile.close();
  }

  if (mission.bRTCoefTest) {  // Reflection/Transmission Test
    RTCoef::RunRTCoefTest(100, 10, 8, 4, 8, 4, 2);
  }

  if (mission.bSourcePatternTest) { // Examine source radiation pattern

    std::cout << "@@ __EVENT_SOURCE_TEST__" << std::endl;

    S2::S2Set * pTOAA   // Degree-3 Tesselsphere for pattern-testing
      = new S2::TesselSphere(S2::TESS_ICO, 3);
    PhononSource::Set_TOA_Array(pTOAA);

    PhononSource S
      = ShearDislocation(MParams.EventSourceMT,
                         ECS.Convert(MParams.EventSourceLoc));

    S.output_differential_probabilities(RAY_NA, RAY_P, "c");
    S.output_differential_probabilities(RAY_NA, RAY_SH, "-");
    S.output_differential_probabilities(RAY_NA, RAY_SV, "y");
                                // Characters tell GMT circle, dash,
                                // or vertical bar, respectively

  }

  if (mission.ModelBuildRequired()) {   // (true if any flags that
    Text IfErrorMsg;                    // require a constructed model
    try{                                // are true)

      IfErrorMsg = "while constructing Earth model:";
      Model Mod(MParams);               // Construct Model

      IfErrorMsg = "during model retrospective output:";
      if (mission.bDumpGrid) {          // Output grid as plaintext
        Mod.GetGridRef().DumpGridToAscii();
      }
      
      if (true) {                       // Output all scattering objects
        Scatterer::PrintAllScatteringStats();
      }

      IfErrorMsg = "during simulation execution:";
      if (mission.bRunSim) {            // RUN SIMULATION
        Mod.RunSimulation();            //
      }

    } catch (std::exception &e) {
      std::cout << "**\n** Error " << IfErrorMsg << "\n"
                << "** What: " << e.what() << "\n** Exiting...\n";
      exit(1);
    }

  }
 
////
}// end main()
//


//////
// FUNCTION:  print_banner()
//
//   Print a startup banner at program start.
//
void print_banner() {

  std::ostream * out = &std::cout;

  *out << "**" << std::endl
       << "**  Radiative3D - "
       << "A code for radiative transport in 3D Earth models" << std::endl
       << "**" << std::endl
       << "**  (c) 2020 Christopher J. Sanborn and the" << std::endl
       << "**           Solid Earth Geophysics Team" << std::endl
       << "**           at the University of Connecticut, Storrs, CT."
                        << std::endl
       << "**           https://github.com/christophersanborn/Radiative3D"
                        << std::endl
       << "**" << std::endl;
  *out << "**  BUILD STATS:  "
       << "Revision number: " << REVISION_NUM << std::endl
       << "**                "
       << "Floating-point representation: " << (8*sizeof(Real))
       << "-bit" << std::endl
       << "**" << std::endl
       << "**" << std::endl;


}


//////
// FUNCTION:  print_help()
//
//   Print a help message
//
void print_help() {

  std::ostream * out = &std::cout;

  *out << "\nRadiative3D Manual Page is available at:\n"
       << "https://github.com/christophersanborn/Radiative3D/blob/master/doc/MANUAL.md\n\n";

  *out << "Recognized command-line options:\n\n";

  CmdOpt::OutputAllRecognizedTokens(out);

  *out << std::endl;

}


//////
// FUNCTION:  process_option()
//
//  Takes a CmdOpt object and sets the relevent model or operational
//  parameters accordingly.  When called in a loop that iterates over
//  a set of CmdOpt objects that have been initialized from argc/argv,
//  this has the effect of "processing" all the command-line
//  arguments.
//
void process_option(CmdOpt & opt, ModelParams & params, 
                    MissionParams & mission) {

  switch (opt.GetID()) {

  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_FREQ:        // *** Frequency:
                                // ***
    //
    params.Frequency = opt.PopValue_Real();
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_NUMBER:      // *** Number of Phonons:
                                // ***
    //
    params.NumPhonons = opt.PopValue_Integer();
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_REPORTS:     // *** Reports (On or Off):
                                // ***
    //
    if (!opt.PeekValue()) {               // If token-only, 
      dataout.SuppressAllReports(false);  // then same as ALL_ON
    }
    else {                                // Else start from ALL_OFF and
      dataout.SuppressAllReports(true);   // add in what's requested
      while (opt.PeekValue()) {
          //
          Text keywd = opt.PopValue_Text();
          if (keywd == "ALL_ON") {dataout.SuppressAllReports(false);}
          else if (keywd == "ALL_OFF") {dataout.SuppressAllReports(true);}
          else if (keywd == "GEN") {dataout.SuppressGenerate(false);}
          else if (keywd == "SCT") {dataout.SuppressScatter(false);}
          else if (keywd == "REF") {dataout.SuppressReflect(false);}
          else if (keywd == "COL") {dataout.SuppressCollect(false);}
          else if (keywd == "CEL") {dataout.SuppressTransfer(false);}
          else if (keywd == "LST") {dataout.SuppressLost(false);}
          else if (keywd == "TMO") {dataout.SuppressTimeout(false);}
          else if (keywd == "INV") {dataout.SuppressInvalid(false);}
          else if (keywd == "SCATTERS") {
            dataout.SuppressGenerate(false);
            dataout.SuppressScatter(false);
            dataout.SuppressReflect(false);
          }
          else {
            throw(Runtime("Valid report keywords are: ALL_ON, ALL_OFF, GEN, "
                   + Text("SCT, REF, COL, CEL, LST, TMO, INV, or SCATTERS.")));
          }
        //
      }// end while
    } // end if
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_REPORT_FILE: // *** Reports Filename
                                // ***
    dataout.SetReportsFile(opt.PopValue_Text());
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_OUTDIR:      // *** Output-file Directory
                                // ***
    dataout.SetOutputDirectory(opt.PopValue_Text());
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_OCSRAW:      // *** Output Coords (OCS)
                                // *** no-transform
    ECS.SetOCSMapping(EarthCoords::OUT_NOTRANSFORM);
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_TTLIVE:      // *** Phonon Time-to-Live:
                                // ***
    //
    params.PhononTTL = opt.PopValue_Real();
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_TOA:         // *** Take-Off Angle Degree:
                                // ***
    //
    params.TOA_Degree = opt.PopValue_Integer();
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_COMPSELECT:  // *** Grid-Compiled Model Selector:
                                // ***
    //
    params.GridSource = ModelParams::GRID_COMPILED;
    params.CompiledSelector = opt.PopValue_Integer(0);
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_MODARGS:     // *** Grid-Compiled Model Args:
                                // ***
    //
    while (opt.PeekValue()) {
      params.CompiledArgs.push_back(opt.PopValue_Real());
    }
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_SEISARRAY:   // *** Seismometer Array:
                                // ***
    {
      //double AZ = opt.PopValue_Real();
      //double R0 = opt.PopValue_Real();
      //double gap =  opt.PopValue_Real();
      //double gather =  opt.PopValue_Real();
      //double n_seis =  opt.PopValue_Real();

      //double angle = (Geometry::Pi/180)*AZ;
      //R3::XYZ uaz (sin(angle), cos(angle), 0);
      //R3::OrthoAxes ori (0,0,0);
      //R3::XYZ R0V (sin(angle)*R0, cos(angle)*R0, 0);
      //for(int j=0; j < n_seis; j++)
      //  {
      //    R3::XYZ loc = uaz.ScaledBy(j*gap) + R0V;
      //    params.AddSeismometerByWavelength(loc,ori,gather);
      //  }
      throw (Runtime(
              "Arg: --seisarray: Disabled; use --seis-p2p instead."
            ));
            // Disabled until I figure out best way to handle arg list
            // given that it is not known at this point what
            // coordinate system was chosen by user.
    }
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_SEIS_P2P:    // *** Point-to-Point Seis Array:
    {                           // ***
      Count ValCount = opt.GetValueCount();
      bool twogathers = (ValCount==10);

      R3::XYZ Origin = opt.PopValue_XYZ();
      R3::XYZ Dest = opt.PopValue_XYZ();
      Real Offset = opt.PopValue_Real();
      Real Gather1 = opt.PopValue_Real();
      Real Gather2 = Gather1;
      if (twogathers) Gather2 = opt.PopValue_Real();
      int nSeis = opt.PopValue_Integer();

      const R3::XYZ Span = Origin.VectorTo(Dest);
      const R3::XYZ Dir = Span.Unit();
      const R3::XYZ Begin = Origin + Dir.ScaledBy(Offset);
      const Real Dist = Begin.VectorTo(Dest).Mag();
      const Real Gap = (nSeis>1) ? (Dist / (nSeis-1))
                                 : 0.0;

      for (int i=0; i<nSeis; i++) {
        R3::XYZ Loc = Begin + Dir.ScaledBy(i*Gap);
        Real rangefrac = (Origin.VectorTo(Loc).Mag())/Span.Mag();
        Real Gather = Gather1 + rangefrac*(Gather2-Gather1);

        EarthCoords::Generic ELoc(Loc.x(),Loc.y(),Loc.z()); // Hack...
        static bool warning_issued = false; //
        if (!warning_issued) {              // Only warn once.
          std::cout << "Warning: Seis array interpolation ignores coordinate "
                    << "system and could produce distorted results.\n";
          warning_issued = true;
        } //TODO: Seismoter allocation needs an overhaul...

        if (twogathers) {
          params.AddSeismometerFixedRadius(ELoc,ModelParams::AX_RTZ,Gather);
        } else {
          // TODO: Remove conditional and retain ONLY fixed-radius call.
          // DEPRECATION NOTE: Use of this option (SEIS_P2P) to
          // specify by wavelength when only one gather radius is
          // given is retained temporarily to not break existing
          // scripts, but this usage is deprecated in favor of using
          // SEIS_P2PW when units should be in wavelengths.
          params.AddSeismometerByWavelength(ELoc,ModelParams::AX_RTZ,Gather);
        }          
      }
    }
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_SEIS_P2PW:   // *** Point-to-Point Seis Array:
    {                           // ***    (units are wavelengths)
      Count ValCount = opt.GetValueCount();
      bool twogathers = (ValCount==10);

      R3::XYZ Origin = opt.PopValue_XYZ();
      R3::XYZ Dest = opt.PopValue_XYZ();
      Real Offset = opt.PopValue_Real();
      Real Gather1 = opt.PopValue_Real();
      Real Gather2 = Gather1;
      if (twogathers) Gather2 = opt.PopValue_Real();
      int nSeis = opt.PopValue_Integer();

      const R3::XYZ Span = Origin.VectorTo(Dest);
      const R3::XYZ Dir = Span.Unit();
      const R3::XYZ Begin = Origin + Dir.ScaledBy(Offset);
      const Real Dist = Begin.VectorTo(Dest).Mag();
      const Real Gap = (nSeis>1) ? (Dist / (nSeis-1))
                                 : 0.0;

      for (int i=0; i<nSeis; i++) {
        R3::XYZ Loc = Begin + Dir.ScaledBy(i*Gap);
        Real rangefrac = (Origin.VectorTo(Loc).Mag())/Span.Mag();
        Real Gather = Gather1 + rangefrac*(Gather2-Gather1);

        EarthCoords::Generic ELoc(Loc.x(),Loc.y(),Loc.z()); // Hack...
        std::cout << "Warning: Seis array interpolation ignores coordinate "
                  << "system and could produce distorted results.\n";
        //TODO: Seismoter allocation needs an overhaul...

        params.AddSeismometerByWavelength(ELoc,ModelParams::AX_RTZ,Gather);
      }
    }
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_SEISBINS:    // *** Time Bins per Cycle:
                                // ***
    //
    params.TimeBinsPerCycle = opt.PopValue_Real();
    params.TimeBinSize = 0;
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_SEISBINSIZE: // *** Time bin width in seconds
                                // ***
    //
    params.TimeBinsPerCycle = 0;
    params.TimeBinSize = opt.PopValue_Real();
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_CYLRANGE:    // *** Cylinder Range:
                                // ***
    //
    params.CylinderRange = opt.PopValue_Real();
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_FLATTEN:     // *** Earth-flattening transform
                                // ***
    //
    ECS.SetEarthFlattening(true);
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_EARTHRAD:    // *** Earth Radius
                                // ***
    //
    ECS.SetEarthRadius(opt.PopValue_Real());
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_EVENT_LOC:   // *** Event Location:
    {                           // ***
      Real loc_x = opt.PopValue_Real();
      Real loc_y = opt.PopValue_Real();
      Real loc_z = opt.PopValue_Real();
      params.EventSourceLoc = EarthCoords::Generic(loc_x,loc_y,loc_z);
    }
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_EVENT_MT:    // *** Source Moment-Tensor:
                                // ***
    {
      Text SourceType = opt.PopValue_Text();

      if (SourceType == "EQ") {
        params.EventSourceMT = Tensor::USGS(0,-1,1,0,0,0);
      }
      else if (SourceType == "EXPL") {
        params.EventSourceMT = Tensor::USGS(1,1,1,0,0,0);
      }
      // OP_EVENT_MT: USGS Option:
      else if (SourceType == "USGS") {
        try {
          double rr = opt.PopValue_Real();
          double tt = opt.PopValue_Real();
          double pp = opt.PopValue_Real();
          double rt = opt.PopValue_Real();
          double rp = opt.PopValue_Real();
          double tp = opt.PopValue_Real();
          params.EventSourceMT = Tensor::USGS(rr,tt,pp,rt,rp,tp);
        } catch(...) {
          throw(Runtime(
            "USGS keyword expects six numeric moment tensor elements."));
        }
      }
      // OP_EVENT_MT: SDR Option:
      else if (SourceType == "SDR") {
        std::vector<Real> MTargs;
        while (opt.PeekValue()) {
          MTargs.push_back(opt.PopValue_Real());
        }
        if(MTargs.size() == 3){         // SDR only
          params.EventSourceMT
            = Tensor::SDR(MTargs.at(0),MTargs.at(1),MTargs.at(2));
        }
        else if (MTargs.size() == 4) {  // SDR, and Iso
          Real iso = MTargs.at(3);
          if (iso > 1.0 || iso < -1.0) {
            iso = Tensor::SDR::IsoFracFromIsoAngle(iso);
          } // Assume iso is an angle not a fraction if outside [-1,1].
          params.EventSourceMT
            = Tensor::SDR(MTargs.at(0), MTargs.at(1),
                          MTargs.at(2), iso);
        }
        else if (MTargs.size() == 5) {  // SDR, Iso, and Moment
          Real iso = MTargs.at(3);
          Real moment = MTargs.at(4);
          if (moment==0.0) {  // Note: New argument order. This catches
            Real temp = iso;  // old scripts using old arg order (unless they
            iso = moment;     // specified a non-zero isofrac, but there
            moment = temp;    // have not been many of those.)
            std::cout         //
              << "Warning: SDR Event Specification: " 
              << "Swapped iso and moment arguments.\n"
              << "Warning: (Note new argument order: "
              << "--source=SDR,strike,dip,rake,iso,moment)\n";
          } // Catch possible confusion over argument order.
          if (iso > 1.0 || iso < -1.0) {
            iso = Tensor::SDR::IsoFracFromIsoAngle(iso);
          } // Assume iso is an angle not a fraction if outside [-1,1].
          params.EventSourceMT
            = Tensor::SDR(MTargs.at(0), MTargs.at(1),
                          MTargs.at(2), iso, moment);
        }
        else {
          throw(Runtime(
            "SDR keyword expects SDR,strike,dip,rake[,iso[,moment]]."));
        }
      }          
      else {
        throw(Runtime(SourceType 
                      + " is not a valid source mechanism keyword."));
      }
    }
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_OVR_MFP:     // *** Override MFPs
    {                           // ***
    Real mfpP = opt.PopValue_Real();
    Real mfpS = opt.PopValue_Real();
    Scatterer::OverrideMFP(mfpP,mfpS);
    } break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPT_NODEFLECT:   // *** No Scattering Deflection
                                //
    Scatterer::SetNoDeflect();
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPTM_HELP:       // *** Print Help Message
                                // ***
    mission.bHelpMsg = true;
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPTM_RUNSIM:     // *** Run Simulation
                                // ***   (on by default, this
                                // ***   counteracts an override)
    mission.bRunSim = true;
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPTM_DUMPGRID:   // *** Plaintext Grid Dump
                                // ***
    mission.bDumpGrid = true;
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPTM_PARAMOUTFN: // *** Output Model Params to Octave File
                                // ***
    mission.bOutputModParamsOctv = true;
    mission.FNModParamsOctv = opt.PopValue_Text();
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPTM_RTTEST:     // *** R/T Coefficient Test
                                // ***
    mission.bRTCoefTest = true;
    mission.bRunSim = false;
    break;


  ////////////////////////////////////////////////////////////////////
  case CmdOpt::OPTM_EVENTTEST:  // *** Event Source Analytics
                                // ***
    mission.bSourcePatternTest = true;
    mission.bRunSim = false;
    break;


  ////////////////////////////////////////////////////////////////////
  default:                      // *** Unrecognized OP_ID:
                                // ***
    //
    throw(Runtime("Unrecognized option.")); 
    break;
  }

}


