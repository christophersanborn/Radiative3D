// model.hpp
//
// This file develops a class that represents the "whole shebang", as
// it were.  This class encapsulates the physical/computational
// "model" of the elastodynamic universe we are trying to simulate.
// This one brings together all the other pieces, and puts them
// together into a functioning system.  The radiative transport "main
// loop" is contained within this class.
//
// This file participates in the global namespace. Since function
// main() is expected to work directly and closely with this object
// class, namespace encapsulation is deemed unnecessary.
//
//
#ifndef MODEL_H_
#define MODEL_H_
//
#include <iostream>
#include <vector>
#include <atomic>
#include "events.hpp"
#include "scatterers.hpp"
#include "grid.hpp"   /* includes ecs.hpp, elastic.hpp, geom.hpp */
#include "media.hpp"  /* includes raytype.hpp */

// In this file:
// CLASS DEFINITIONS FOR:
//
//   o  Class ModelParams
//   o  Class Model
//
// Search on "CLASS:" to jump between class definitions in this 
// file.
//

// *** FORWARD DECLARATIONS:
//
class ModelParams;
class Model;

// *** TYPEDEFS:
//

// *** CLASS DEFINITIONS:
//

//////
// CLASS:   ::: ModelParams :::
// ENCAPS:  Model Parameters
//
//  Encapsulates the set of parameters that define the Earth model or
//  that direct the initialization of the model. Usually, the members
//  of this class will be set by, say, function main() based on the
//  set of command line args passed by the user.  Then the whole set
//  of parameters can be passed en masse to the Model() constructor.
//
//  A secondary purpose of this class is to assign reasonable defaults
//  to the parameters if they remain unset by the user.
//
class ModelParams {

  friend class Model; // Model initialization code needs access to
                      // model parameters.

public:

  // :::::::::::::::::::::::::::::::::::::::::
  // ::: Member Types  (ModelParams Class) :::
  // :::::::::::::::::::::::::::::::::::::::::

  enum grid_source_e {
        GRID_UNSPEC,     // Unspecified by user.
        GRID_FROMFILE,   // Read grid from file.
        GRID_COMPILED,   // User-written Grid::ConstructGridManual()
                         // function (in user.cpp).
        GRID_UNIT_CYL,   // Built-in Unit test based on cylinder model cells.
        GRID_UNIT_TETRA, // Built-in Unit test based on tetra model cells.
  };                     // 

  enum axes_scheme_e {   // Seismometer Axes Schemes
        AX_ENU,          // Orient to East, North, Up
        AX_RTZ           // Orient to Radial, Transverse, Up(Z) w.r.t. Source
  };                     //

private:
  struct SeisRequest {
    EarthCoords::Generic   Location;  // Where is the seismometer
    axes_scheme_e       Orientation;  // Three axes of seismometer
    Real GatherRadiusInner[RAY_NBT];  // Inner radius of gather zone
    Real GatherRadiusOuter[RAY_NBT];  // Outer radius ''   ''    ''
    bool   RadiiUnitsAreWavelengths;  // TRUE: Radii in wavelengths
                                      // FALSE: Radii in same units
                                      // as model uses (probably km)
  };


public:
  ;
  // :::::::::::::::::::::::::::::::::::::::::::::
  // ::: Member Variables  (ModelParams Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::
  //
  //                            Access to these members is mostly
  //                            'public' because this is a
  //                            facilitator class for passing lots
  //                            of parameters to the Model
  //                            constructor.
  //
  //                            Some members may be 'private'
  //                            because their values require
  //                            validity checking.  In this case,
  //                            Set() methods will be provided.
  //

  Tensor::Tensor EventSourceMT;  // Event Source Moment Tensor
  EarthCoords::Generic EventSourceLoc;  // Event Source Location
  int               TOA_Degree;  // Determines number of Take-off
                                 //  Angles on which to define event
                                 //  source functions
  long              NumPhonons;  // Number of phonons to cast
  Real               PhononTTL;  // Phonon Time-to-Live, this is the
                                 //  max propagation sim-time before a
                                 //  phonon is to be abandoned
  Real               Frequency;  // Phonon frequency in Hertz
  Real        TimeBinsPerCycle;  // Establishes bin-width for the
                                 //  seismometer objects relative to
                                 //  the Phonon frequency.
  Real             TimeBinSize;  // Bin width in seconds. TimeBinSize
                                 //  and TimeBinsPerCycle are mutually
                                 //  exclusive, use GetBinSize()
                                 //  accessor to ensure correct
                                 //  readout.
                                 //
  grid_source_e     GridSource;  // Where do we get the grid from?
                                 //
  Real           CylinderRange;  // The range parameter used by the
                                 // simplified "cylinder"-type quasi-
                                 // 1D earth models (ignored for
                                 // tetra-based models).
                                 //
  int         CompiledSelector;  // Int value to be passed to user-
                                 // compiled grid builder to select
                                 // among multiple models.
  std::vector<Real>              // A set of user-defined arguments to
         CompiledArgs;           // be passed to user-written
                                 // "compiled" model grids.

  //
  //                            The following members are private
  //                            because they need access methods to
  //                            set.
  //

private:

  std::vector<SeisRequest>        // The list of requested
                 mSReqList;       // seismometers


public:

  // :::::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (ModelParams Class) :::
  // :::::::::::::::::::::::::::::::::::::::::

  ModelParams() :
    EventSourceMT    ( Tensor::USGS(0,-1,1,0,0,0) ),
    EventSourceLoc   ( 0,0,-1              ),
    TOA_Degree       ( 7                   ),
    NumPhonons       ( 10                  ), // 
    PhononTTL        ( 60.0                ), // Seconds
    Frequency        ( 4.0                 ), // Hertz
    TimeBinsPerCycle ( 0.0                 ), // Zero signals unset
    TimeBinSize      ( 2.0                 ), //
    GridSource       ( GRID_UNSPEC         ),
    CylinderRange    ( 600.0               ), // 600 km
    CompiledSelector ( 0                   )
  {}

  ~ModelParams() {}


  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Set Methods  (ModelParams Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  void AddSeismometerByWavelength(EarthCoords::Generic location, 
                                  axes_scheme_e orientation,
                                  Real radius_wl);

  //            Creates a request for a Seismometer in which the gather
  //            radius is specified in wavelengths instead of the usual
  //            kilometers. This results in a seismometer whose P and S
  //            gather radii are different.
  //

  void AddSeismometerFixedRadius(EarthCoords::Generic location, 
                                 axes_scheme_e orientation,
                                 Real radius);
  //
  //            Creates a request for a Seismometer with a fixed gather
  //            radius, specified in the usual length units (probably
  //            kilometers).  In this case, both the P and S gather
  //            radii are the same.
  //

  void AddSeismometerRing(EarthCoords::Generic location, 
                          axes_scheme_e  orientation,
                          Real inner_radius,
                          Real outer_radius);
  //
  //            Creates a request for a Seismometer with both an inner
  //            and an outer gather radius.  This results in
  //            "ring"-shaped gather region for phonons.  Can be useful
  //            for binning on range rather than azimuth.  Radii
  //            specified in extensive units (probably kilometers), and
  //            are the same for both P and S phonons.
  //

  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (ModelParams Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  Real GetBinSize() const {     // Returns TimeBinSize unless TimeBinsPerCycle
    if (TimeBinsPerCycle==0) {  // is set (ie is non-zero).  Resolves mutual
      return TimeBinSize;       // exclusivity between TBS and TBPC.
    } else {                    //
      return (1.0/(Frequency*TimeBinsPerCycle));
    }
  }

  Real GetBinsPerCycle() const {            // (Answer will agree with
    return (1.0/(Frequency*GetBinSize()));  //  GetBinSize())
  }                                         //


  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Data Reveal Methods  (ModelParams Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::

  void Output() const;
  void OutputOctaveText(std::ostream *) const;

};


//////
// CLASS:   ::: Model :::
//
// ENCAPS:  
//
//   Encapsulates a model of our simulated physical universe.  This
//   object, in a way, coordinates the actions of all the others.
//   It's primary jobs are:
//
//     o  Initialize global and class-static objects and variables,
//        including the discrete array(s) of take-off angles, etc.
//
//     o  Initializes the random number generator
//
//     o  Build and initialize the collection of MediumCell objects
//        that will represent our physical Earth model
//
//     o  Build the PhononSource object to represent our event source,
//        and link it to the MediumCell in which it resides.
//
//     o  Provide a set of methods to run one or more simulation
//        scenarios, eg. "spray X-number of phonons," etc.
//
//
// INPUTS:
//
//   The constructor takes a ModelParams object as input. This
//   contains a wealth of information to direct the construction and
//   simulation of our Earth model and simulation scenarios.  The
//   contents of the ModelParams class are largely determined by
//   either options provided on the command-line or by the defaults
//   hard-coded in the ModelParams class.  Among the parameters might
//   be the filenames of various Earth-model data files which will
//   provide further input to this class for model generation or other
//   purposes.
//
// LIMITATIONS:
//
//   Because the constructor is responsible for setting "global"
//   (actually class-static) parameters in other classes (eg. the the
//   phonon frequency in the ScatterParams class), there should never
//   be more than one Model object instantiated in any single runtime
//   environment, as they cannot coexist without stepping on each
//   other's toes.  This should not be a major limitation, as at
//   present I cannot think of any reason to need more than one model
//   object in a single running instance of the program.
// 
class Model {
private:;

  // :::::::::::::::::::::::::::::::::::::::
  // ::: One-off Typedefs  (Model Class) :::
  // :::::::::::::::::::::::::::::::::::::::

  typedef std::vector<MediumCell*>    // Array of MediumCell Pointers
                MediumCellPtrArray;   // 

  typedef std::vector<CellFace*>      // Array of CellFace Pointers
                  CellFacePtrArray;   //


  // :::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Private Member Variables  (Model Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::

  S2::S2Set        * mpTOA_Array;   // Take-off Angle Array

  ShearDislocation * mpEventSource; // Event Source (Current
                                    // incarnation of the code
                                    // supports only a single event
                                    // source.  But that should be
                                    // all that's needed.)

  Grid                     mGrid;   // The lattice of GridNodes that
                                    // define the geometry and
                                    // physical properties of the
                                    // model.

  MediumCellPtrArray  mCellArray;   // Cell array (pointers - the actual
                                    // cells will be constructed at
                                    // runtime).  These will comprise
                                    // the actual Earth model.

  CellFacePtrArray mSurfaceFaces;   // We maintain a list CellFaces that
                                    // comprise the Earth's surface
                                    // (used for surface-locating for,
                                    // e.g., seismometer placement).

  long               mNumPhonons;   // Number of phonons to generate
                                    // and then simulate.

  ///
  // Threads and coordination:

  std::atomic<long> mPhononsRemain; // Threads deduct from this counter
                                    // until all phonons simulated.
  void SimulationThread();


public:

  // :::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Model Class) :::
  // :::::::::::::::::::::::::::::::::::

  Model(const ModelParams & par = ModelParams());

  ~Model();


  // :::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (Model Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::
  //
  //            (And other methods that return
  //             information about the model)
  //

  MediumCell * FindCellContainingPoint(const R3::XYZ & loc) const;
  //                Returns a pointer to the MediumCell object
  //                that "best" contains the R3 point 'loc'.

  R3::XYZ FindSurface(R3::XYZ loc) const;
  //                Find point on Earth surface directly above/below the
  //                given location (used to pin seismometers to surface
  //                when their supplied z component might not perfectly
  //                place them at surface level.)


  // :::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Internals Reveal Methods  (Model Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::

  const Grid & GetGridRef() {   // Reveal a read-only ref to grid.
    return mGrid;               // (Used so main() can ask the grid
  }                             // to print itself for diag purposes.)


  // ::::::::::::::::::::::::::::::::::::::::::::
  // ::: Sim Execution Methods  (Model Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::

  void RunSimulation();   // Run the simulation outer loop.

  void RunScatterPatternTest(ScatterParams SPar, 
                             raytype inray, long count) const {
    Scatterer Scat(SPar);
    Scat.test_random_rayset(inray, count);
  }

  void RunEventPatternTest(long count) const {
    mpEventSource->output_random_rayset(count);
  }


private:

  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Internal Helper Functions  (Model Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::

  void BuildCellArray_Cylinder();
                    // One of several model-building helpers.  This
                    // one builds up the cell array of Cylinder
                    // cell-types.

  void BuildCellArray_WCGTetra();
  void BuildCellArray_SphericalShells();

  static Index WCGBaseIndexFromIJK(Index i, Index j, Index k,
                                   Count nI, Count nJ, Count nK);

  void WCGBuildBasicPattern(const std::vector<const GridNode *> & nodes,
                            bool mirror, bool isSurface);

  void WCGLinkBlocksForward(Index Block, Index Adjacent,
                            CellFace::face_id_e faceid, bool mirror);


};

////
///
#endif //#ifndef MODEL_H_
//
