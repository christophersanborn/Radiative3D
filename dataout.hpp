// dataout.hpp
//
// This file develops the DataReporter class, and declares a global
// DataOut object based on that class, to process data output of all
// simulation results in the form desired by the user.  The simulation
// inner loops call upon this object to report all events of POTENTIAL
// interest, and then this object determines which events to actually
// report, and in what format to report them, based on the wishes of
// the user.  An appropriate set of ReportXXX() methods will be
// developed for use by the simulation code, and a set of methods that
// control reporting preferences will also be developed.  Ultimately,
// code in main.cpp will translate command-line args into these
// preferences.
//
// This header file will need to be included in the .cpp
// implementation files of any modules that use the DataOut object,
// (the phonon module being the prime example), but should not need to
// be included by any other header files.
//
#ifndef DATAOUT_H__
#define DATAOUT_H__
//
#include <string>
#include <iostream>
#include <fstream>
#include "geom.hpp"
#include "tensors.hpp"
#include "raytype.hpp"

//////
// CLASSES: -- Forward Declarations --
//
//   FROM OTHER HEADERS:  (Referenced here by pointer only - no
//                         need to include full header.)

class Phonon;      /* Defined in phonons.hpp */
class MediumCell;  /* Defined in media.hpp */


//////
// CLASS: Seismometer
//
//   Encapsulates a seismometer, records a time-series of energy
//   accumulator bins. Responsible for resolving incoming phonons into
//   their components of motion.  Binned data can be used to produce
//   seismic envelopes, or, (with additional bookkeeping not yet
//   implemented), full-waveform amplitude synthetics.
//
//   The dimensional quantity actually accumulated in the bins is
//   Energy per Time^2 per Area.  An arriving phonon is assumed to
//   have an Energy-per-Time proportional to it's amplitude-squared.
//   This energy-per-time is then divided by the gather-area of the
//   seismometer and the temporal width of a single time bin, to
//   normalize on the probability of catching a single phonon in a
//   particular bin by the particular seismometer.
//
//   The signal reported by the seismometer at simulation-end is then
//   suitable to use for calculations such as the amount of energy
//   incident on a given target (say, a building foundation),
//   according to the following procedure: (1) Integrate the signal in
//   time over the source-time function of the source event, (2)
//   Integrate in time over the period that the target absorbs energy,
//   (3) Multiply by the area of the target (or integrate over area of
//   target if probing the signal at differantial locations makes
//   sense... which it probably wouldn't for "small" targets).  The
//   result should then be an energy quantity.
//
class Seismometer {
protected:

  // ::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Supporting Data Types  (Seismometer Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::

  enum axes_e {AXIS_X, AXIS_Y, AXIS_Z, NUM_AXES};

  class BinRecord {  // Encapsulates a single time record in our
  public:            // seismic trace. Energy is binned by time window
                     // and distributed by seismic axis and raytype.
    Real        mEnergyAxes[NUM_AXES];          // Energy by axis
    Real      mEnergyByType[RAY_NUMBASICTYPES]; // Energy by raytype
    unsigned   mCountByType[RAY_NUMBASICTYPES]; // Phonon count by rt

    BinRecord() {               // Make sure bins get initialized
      mEnergyAxes[AXIS_X] = 0;  // to zero on construction
      mEnergyAxes[AXIS_Y] = 0;
      mEnergyAxes[AXIS_Z] = 0;
      mEnergyByType[RAY_P] = 0;
      mEnergyByType[RAY_S] = 0;
      mCountByType[RAY_P] = 0;
      mCountByType[RAY_S] = 0;
    }
  };


protected:

  // :::::::::::::::::::::::::::::::::::::::::::::
  // ::: Member Variables  (Seismometer Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::

  R3::XYZ                mLoc;  // Location of seismometer
  R3::XYZ             mAxesX1;  // X1 Axis (East-West, eg.)
  R3::XYZ             mAxesX2;  // X2 Axis (North-South, eg.)
  R3::XYZ             mAxesX3;  // X3 Axis (Vertical, eg.)
                                //     Constructor responsible for
                                //     ensuring axes are orthonormal.
  std::string       mAxesDesc;  // Alpha code describing axes orientation, to
                                // be included in output files. A freeform
                                // string, but typical usage expects a three
                                // letter code such as "ENZ" for East, North,
                                // Up, or "RTZ" for Radial, Transverse, Z.
                                // Default is "UNK" for Unknown.

  Real             mBeginTime;  // Time when recording "turns on"
  Real      mRadiusO[RAY_NBT];  // Gather radius for phonons, indexed by
                                // raytype.  (Outer radius)
  Real      mRadiusI[RAY_NBT];  // Gather radius for phonons, indexed by
                                // raytype.  (Inner Radius)
  Real         mArea[RAY_NBT];  // Area of the gather region. Used for
                                // dimensional scaling of energy bins.
  bool           mPassthrough;  // True if this is a "passthrough" seismometer,
                                // meaning its gather area might overlap with
                                // another, and therefore a "caught" phonon
                                // should continue searching remainder of the
                                // seismometer list for additional catches.
                                // (Default is 'true')

  BinRecord       * mTimeBins;  // ARRAY of BinRecs. Comprises the
                                // seismic trace. Dynamically
                                // allocated.

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Class-Static Member Variables  (Seismometer Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //
  //        These values need runtime initialization prior to
  //        the construction of any Seismometer objects.  The
  //        responsibility to do this initialization lies with
  //        the Model constructor in model.cpp.
  //

  static Real    cmTimePerBin;  // Amount of time (seconds) in each bin
  static Count      cmNumBins;  // Number of bins in a time series

  static R3::XYZ   cmEventLoc;  // Knowing event location allows
                                // azimuths and ranges to be
                                // calculated for seismometers.
  static Tensor::Tensor cmEventMT;  // Event moment-tensor is not
                                    // needed internally, but we
                                    // include it as metadata when we
                                    // output seismic traces for use
                                    // in plotting and annotating.


public:
  ;  
  // :::::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Seismometer Class) :::
  // :::::::::::::::::::::::::::::::::::::::::

  Seismometer(const R3::XYZ loc,        // Seismometer Location
              const R3::XYZ axes_X1,    // Direction of X1 Axis
              const Real inner_radius[RAY_NBT],
              const Real outer_radius[RAY_NBT],
              std::string axes_desc = "UNK");
              // Note: X3 Axis is assumed to be "Up", thus only neces-
              // sary to provide X1 axis, as X2 can be compute from X1
              // and X3.

  ~Seismometer();


  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Set Methods  (Seismometer Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  static void SetTimeBinParameters(Real bintime, Real recordtime) {
    cmTimePerBin = bintime;
    cmNumBins    = floor(recordtime / bintime);
  }

  static void SetEventParameters(R3::XYZ loc, Tensor::Tensor MT) {
    cmEventLoc = loc;
    cmEventMT = MT;
  }
  

  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Do-Something Methods  (Seismometer Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  bool CatchPhonon(const Phonon & phon);
  //            Checks whether the Phonon is within the gather radii,
  //            and if it is, we ingest it and increment the appropriate
  //            bins.  Returns true if we ingested the phonon, and false
  //            if we either (a) did NOT ingest or (b) we DID ingest,
  //            but this Seismometer is marked as a "passthrough".
  //

  void OutputDescription(std::ostream * out, const std::string & IDstr);
  //            Write out a plaintext self-description of this seismometer
  //            (location and other properties, etc.) to the output stream
  //            pointed to by 'out'.  IDstr is a string label prepended to
  //            each output line for post-op filtering.
  //

  void OutputTrace(std::ostream * out, const std::string & IDstr);
  //            Outputs a seismic trace by writing out the data from the
  //            time bins, one record per line.
  //

  void OutputOctaveText(std::ostream * out);
  //            Outputs the contents of a Seismometer object,
  //            including both meta and trace data, in GNU/Octave's
  //            native "text" format.
  //


};


//////
// CLASS: DataReporter
//
class DataReporter {
protected:

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Private Member Variables  (DataReporter Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::

  // *** Files and Filenames:
  //

  std::ostream *  mposReports;  // Output stream for reports (points
                                // to either stdout or an ofstream).
  std::ofstream   mofsReports;  // This is the file, if stdout not
                                // requested.

  std::string        msOutDir;  // Output Directory Name

  // *** Report Stream:
  //

  std::string     mIDGenerate; // Text ID label for GENERATE reports
  bool       mbReportGenerate; // True if GENERATE data wanted

  std::string     mIDCollect;  // Text ID label for COLLECT reports
  bool       mbReportCollect;  // True if COLLECT data wanted

  std::string     mIDReflect;  // Text ID label for REFLECT reports
  bool       mbReportReflect;  // True if REFLECT data wanted

  std::string     mIDTransfer; // Text ID label for REFLECT reports
  bool       mbReportTransfer; // True if REFLECT data wanted

  std::string     mIDScatter;  // Text ID label for SCATTER reports
  bool       mbReportScatter;  // True if SCATTER data wanted

  std::string     mIDLost;     // Text ID label for LOST reports
  bool       mbReportLost;     // True if LOST data wanted
  unsigned long  mNumLost;     // How many phonons discarded as LOST

  std::string     mIDTimeout;  // Text ID label for TIMEOUT reports
  bool       mbReportTimeout;  // True if TIMEOUT data wanted
  unsigned long  mNumTimeout;  // How many phonons discarded for TIMEOUT

  std::string     mIDInvalid;  // Text ID label for INVALID reports
  bool       mbReportInvalid;  // True if INVALID data wanted
  unsigned long  mNumInvalid;  // How many phonons discarded as INVALID

  std::ostream *mpOSSeisTrace; // Output stream for seismometer traces
  std::string    mIDSeisTrace; // Text ID label for seismometer traces
  std::ofstream mOSFSeisTrace; // Output FILE stream for "        "
  bool      mbReportSeisTrace; // True if seismometer output wanted

  std::vector<Seismometer*> mSeismometers;  
            // The Seismometer collection (A vector of pointers.  We
            // grow the vector with the AddSeismometer() member, which
            // is responsible for allocating Seismometers with "new".
            // We shrink it in the destructor, which deletes the
            // Seismometers.)


public:

  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (DataReporter Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::

  DataReporter() :
    mposReports      ( &std::cout ),  // Default to stdout; can override
    msOutDir         ( "."        ),  // Default to current dir
    mIDGenerate      ( "GEN: "    ),
    mbReportGenerate ( true       ),
    mIDCollect       ( "COL: "    ),
    mbReportCollect  ( true       ),
    mIDReflect       ( "REF: "    ),
    mbReportReflect  ( true       ),
    mIDTransfer      ( "CEL: "    ),
    mbReportTransfer ( true       ),
    mIDScatter       ( "SCT: "    ),
    mbReportScatter  ( true       ),
    mIDLost          ( "LST: "    ),
    mbReportLost     ( true       ),
    mNumLost         ( 0          ),
    mIDTimeout       ( "TMO: "    ),
    mbReportTimeout  ( true       ),
    mNumTimeout      ( 0          ),
    mIDInvalid       ( "INV: "    ),
    mbReportInvalid  ( true       ),
    mNumInvalid      ( 0          ),
    mpOSSeisTrace    ( &std::cout ),  // Default to stdout; can override
    mIDSeisTrace     ( "SEIS: "   ),
    mbReportSeisTrace (true       )
  {
    mOSFSeisTrace.open("seis_traces_asc.dat");
    mpOSSeisTrace = &mOSFSeisTrace;
  }

  ~DataReporter() {
    for (unsigned i = 0; i < mSeismometers.size(); i++) {
      delete mSeismometers[i];  // Delete the Seismometer collection
    }
    mOSFSeisTrace.close();
    mofsReports.close();
  }


  // ::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Set Methods  (DataReporter Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::

  void SuppressAllReports(bool suppress = true);
  void SuppressGenerate(bool suppress = true) {mbReportGenerate = !suppress;}
  void SuppressScatter(bool suppress = true) {mbReportScatter = !suppress;}
  void SuppressReflect(bool suppress = true) {mbReportReflect = !suppress;}
  void SuppressCollect(bool suppress = true) {mbReportCollect = !suppress;}
  void SuppressTransfer(bool suppress = true) {mbReportTransfer = !suppress;}
  void SuppressLost(bool suppress = true) {mbReportLost = !suppress;}
  void SuppressTimeout(bool suppress = true) {mbReportTimeout = !suppress;}
  void SuppressInvalid(bool suppress = true) {mbReportInvalid = !suppress;}


  void SetOutputDirectory(std::string);

  void SetReportsFile(std::string);


  // ::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (DataReporter Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::

  std::string GetOutputDirectory() const {
    return msOutDir;
  }


  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Do-Something Methods (DataReporter Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  void AddSeismometer(R3::XYZ loc, R3::XYZ axisX1,
                      const Real inner_radius[RAY_NBT],
                      const Real outer_radius[RAY_NBT],
                      std::string axes_desc) {
      mSeismometers
        .push_back(new Seismometer(loc, axisX1, 
                                   inner_radius, 
                                   outer_radius, axes_desc));
  }

  void OutputPostSimSummary();
  //            Spits out a summary of all aggregate data collected
  //            during the simulation run.  Such aggregate data
  //            includes, e.g., the seismic traces.
  //


  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Data Report Methods  (DataReporter Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::
  //
  //            These are the methods through which the
  //            model simulation code reports events to the
  //            DataReporter object.
  //

  void ReportNewEventPhonon(const Phonon &);
                                  // Phonon generated at event source

  void ReportPhononTimeout(const Phonon &);
                                  // Phonon exceeded lifetime

  void ReportScatterEvent(const Phonon &);
                                  // Phonon scattered

  void ReportPhononCollected(const Phonon &);
                                  // Phonon arrival at collection face

  void ReportReflection(const Phonon &);
                                  // Phonon reflected at cell boundary

  void ReportCellToCell(const Phonon &);
                                  // Phonon left one cell entered another

  void ReportLostPhonon(const Phonon &);
                                  // Phonon exited model

  void ReportInvalidPhonon(const Phonon &);
                                  // Phonon accrued negative or NaN time or
                                  // length or failed to advance in last N
                                  // iterations.

private:

  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Private Utility Methods  (DataReporter Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::

  void output_phonon_dataline(std::ostream *,
                              const std::string &,
                              const Phonon &);


}; // class DataReporter
////
extern DataReporter dataout;  // Global Object
///                           // (initialized in global.cpp)

///
#endif //#ifndef DATAOUT_H__
//
