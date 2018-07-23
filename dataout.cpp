// dataout.cpp
//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "dataout.hpp"
#include "phonons.hpp"
#include "media.hpp"
#include "ecs.hpp"

// In this file:
// CLASS IMPLEMENTATIONS FOR:
//
//   o  Class Seismometer
//   o  Class DataReporter
//
// Search on "&&&&" to jump between class implementations in this
// file.
//


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  Seismometer                                         ****
// ****                                                              ****
//

//
// Static Member Initialization:  (Seismometer Class)
//
Real   Seismometer::cmTimePerBin = 1.0;   // Responsibility to set
Count  Seismometer::cmNumBins    = 100;   // these at runtime lies
R3::XYZ Seismometer::cmEventLoc;          // with Model constructor.
Tensor::Tensor Seismometer::cmEventMT;    //


//////
// CONSTRUCTOR:   Seismometer()
//
Seismometer::Seismometer(const R3::XYZ loc, 
                         const R3::XYZ axes_X1,
                         const Real inner_radius[RAY_NBT],
                         const Real outer_radius[RAY_NBT],
                         std::string axes_desc) :
  mLoc          ( loc       ),
  mAxesDesc     ( axes_desc ),
  mBeginTime    ( 0.0       ), 
  mPassthrough  ( true      )
{
  mAxesX3 = ECS.GetUp(loc);             // Upwards, determined by ECS
  mAxesX2 = mAxesX3.Cross(axes_X1);     // Mutually perpendicular, unscaled
  mAxesX2 = (mAxesX2.IsSquaredZero()) ? ECS.GetNorth(loc) // Fall back on North
                                      : mAxesX2.Unit();   // Or normalize
  mAxesX1 = mAxesX2.Cross(mAxesX3);     // All three now orthonormal.

  mRadiusI[RAY_P] = inner_radius[RAY_P];
  mRadiusI[RAY_S] = inner_radius[RAY_S];
  mRadiusO[RAY_P] = outer_radius[RAY_P];
  mRadiusO[RAY_S] = outer_radius[RAY_S];

  mArea[RAY_P] = (  mRadiusO[RAY_P]*mRadiusO[RAY_P]
                  - mRadiusI[RAY_P]*mRadiusI[RAY_P] ) * Geometry::Pi;
  mArea[RAY_S] = (  mRadiusO[RAY_S]*mRadiusO[RAY_S]
                  - mRadiusI[RAY_S]*mRadiusI[RAY_S] ) * Geometry::Pi;

  mTimeBins = new BinRecord[cmNumBins];   // Heap allocation of the
                                          // seismic trace.

}


//////
// DESTRUCTOR:   ~Seismometer()
//
Seismometer::~Seismometer() {
  delete [] mTimeBins;
}


//////
// METHOD:  Seismometer :: CatchPhonon()
// 
//   Catch and ingest a phonon, incrementing appropriate bins, if
//   the phonon falls within the gather radii.
//
// CALLED FROM:
//
//   Called from DataReporter object whenever a Phonon interacts
//   with a "Collection" face.
//
// RETURN VALUE:
//
//   Returns true if phonon is "caught," unless the Seismometer is a
//   "pass-through" Seismometer, in which case we return false
//   regardless.  (Generally, the DataReporter stops the Seismometer
//   search on a hit, but the pass-through flag allows for overlapping
//   seismometers to gather within a common region.  It is important
//   that all pass-through seismometers are front loaded in the
//   Seismometer list.)
//
bool Seismometer::CatchPhonon(const Phonon & phon) {

  bool retval;
  bool within_window = true;    // Tracks whether inside time window
  bool within_radius = true;    // Tracks whether inside gather radius
                                //
                                //    (Assume both true; then falsify)
                                //
  raytype ph_type = phon.GetRaytype();

  //
  // ::   Determine bin index:
  //
  //    : First, compute corrected arrival time:
  //
  //      (Colision-time, i.e. when the ray hits the surface, is not
  //       equivalent to the time when the wavefront would interact
  //       with the seismometer. This is because we collect phonons
  //       over a non-zero "gather radius", and the ray is likely not
  //       a "direct-hit" to the seismometer.  We assume there is a
  //       parallel ray that intersects the seismometer exactly, and
  //       compute the extra distance that ray has to travel to make
  //       the intersection.  We are using the assumption that the ray
  //       is an approximate represention of a plane wave, which
  //       should be valid over distance scales of about a wavelength
  //       or two.
  //
  //       NOTE: We want to make sure we DON'T do the distance
  //       correction if this is a "ring" seismometer (with both an
  //       inner and outer gather radius), since this type of virtual
  //       seismometer is NOT based on the assumption that the mLoc
  //       member codes the actual seismometer location.)
  //
  Real arv_time = phon.GetTimeAlive();
  Real correction = 0;

  if (mRadiusI[ph_type] <= 0) {         // Skip for ring seismometer
    R3::XYZ ToSeism = phon.GetLocation().VectorTo(mLoc);
    correction = ToSeism.Dot(phon.GetDirection());  // Distance correction
    correction = correction / phon.Velocity();      // Time correction
  }
  arv_time += correction;

  Real ScaledTime = ((arv_time - mBeginTime) /
                     cmTimePerBin);     // In units of "bin size" and
                                        // relative to the recording
                                        // start-time offset
  if (ScaledTime < 0.0) 
    {within_window = false;}

  Index binindex = std::floor(ScaledTime);      // Convert to integer
                                                // (negatives will wrap
                                                // around)
  if (binindex >= cmNumBins)
    {within_window = false;}

  //
  // ::   Check intersection with gather region:
  //
  Real Distance = phon.GetLocation().DistFrom(mLoc);
  if (Distance > mRadiusO[ph_type])
    {within_radius = false;}
  if (Distance < mRadiusI[ph_type])
    {within_radius = false;}

  //
  // ::   Determine return value:
  //
  retval = within_radius;
  if (mPassthrough)         // mPassthrough overrides default
    {retval = false;}       // behavior

  //
  // :: Bail if outside time-window or outside the radii window
  // :: (nothing to record in that case)
  //
  if (!within_window) {return retval;}
  if (!within_radius) {return false;}

  //
  // :: Record energy by component axis:
  //
  R3::XYZ dopm = phon.DirectionOfMotion();  // (Dir of particle motion)
  Real xfrac = dopm.Dot(mAxesX1);   //
  Real yfrac = dopm.Dot(mAxesX2);   //
  Real zfrac = dopm.Dot(mAxesX3);   //
  xfrac *= xfrac;                   // Energy fractions along axes
  yfrac *= yfrac;                   //
  zfrac *= zfrac;                   //
  Real amp = phon.GetAmplitude();   // Amplitude (units of root energenicity)
  Real energy   = amp*amp;          // Energy
  energy /= cmTimePerBin;           // Energy per time
  energy /= mArea[ph_type];         // Energy per time per area
  Real energy_x = energy * xfrac;   // Component fractions
  Real energy_y = energy * yfrac;   //
  Real energy_z = energy * zfrac;   //

  mTimeBins[binindex].mEnergyAxes[AXIS_X] += energy_x;
  mTimeBins[binindex].mEnergyAxes[AXIS_Y] += energy_y;
  mTimeBins[binindex].mEnergyAxes[AXIS_Z] += energy_z;

  //
  // :: Record energy by ray type:
  //
  mTimeBins[binindex].mEnergyByType[ph_type] += energy;

  //
  // :: Record Phonon-count by ray type:
  //
  mTimeBins[binindex].mCountByType[ph_type] += 1;

  return retval;

}


//////
// METHOD:  OutputDescription()
//
void Seismometer::OutputDescription(std::ostream * out, 
                                    const std::string & IDstr) {
  using std::setw;
  EarthCoords::Generic oloc = ECS.OutConvert(mLoc);

  *out << IDstr
       << "Location:  (X,Y,Z) = ( "
       << oloc.x1() << ", "
       << oloc.x2() << ", "
       << oloc.x3() << " )\n";
  *out << IDstr
       << "Orientation: \n";
  *out << IDstr
       << "Gather Radius (Outer):  (P: "
       << setw(14) << mRadiusO[RAY_P] << " )  (S: "
       << setw(14) << mRadiusO[RAY_S] << " )\n";
  *out << IDstr
       << "Gather Radius (Inner):  (P: "
       << setw(14) << mRadiusI[RAY_P] << " )  (S: "
       << setw(14) << mRadiusI[RAY_S] << " )\n";


}

//////
// METHOD:  OutputTrace()
//
void Seismometer::OutputTrace(std::ostream * out, 
                              const std::string & IDstr) {
  using std::setw;

  *out << IDstr
       << "Bin-size is "
       << cmTimePerBin << " seconds.  "
       << cmNumBins << " bin records follow:  ("
       << cmTimePerBin*cmNumBins << " seconds trace time)\n";
  *out << IDstr
       << "#### BEGIN TRACE ####\n";

  for (unsigned int i = 0; i < cmNumBins; i++) {
    *out << "    "
         << setw(14) << mTimeBins[i].mEnergyAxes[AXIS_X]
         << setw(14) << mTimeBins[i].mEnergyAxes[AXIS_Y]
         << setw(14) << mTimeBins[i].mEnergyAxes[AXIS_Z]
         << setw(14) << mTimeBins[i].mEnergyByType[RAY_P]
         << setw(14) << mTimeBins[i].mEnergyByType[RAY_S]
         << setw(10) << mTimeBins[i].mCountByType[RAY_P]
         << setw(10) << mTimeBins[i].mCountByType[RAY_S]
         << std::endl;
  }

  *out << IDstr
       << "#### END TRACE ####\n";

}

//////
// METHOD:  OutputOctaveText()
//
//   Writes contents of the Seismometer object, including meta and
//   trace data, in a format readable by GNU/Octave.
//
void Seismometer::OutputOctaveText(std::ostream * out) {

  using std::setw;
  EarthCoords::Generic oloc = ECS.OutConvert(mLoc);
  EarthCoords::Generic eloc = ECS.OutConvert(cmEventLoc);
  EarthCoords::Generic oax1 = ECS.OutConvertDirectional(mLoc,mAxesX1);
  EarthCoords::Generic oax2 = ECS.OutConvertDirectional(mLoc,mAxesX2);
  EarthCoords::Generic oax3 = ECS.OutConvertDirectional(mLoc,mAxesX3);

  *out << "# Generated by Radiative3D for input into GNU/Octave" 
       << std::endl;
  *out << "# name: SEIS\n"                      // Struct Header
       << "# type: scalar struct\n"             //
       << "# ndims: 2\n"
       << "1 1\n" 
       << "# length: 14 \n\n";
  *out << "# name: Location\n"                  // Seismometer Location
       << "# type: matrix\n"                    //
       << "# rows: 1\n"
       << "# columns: 3\n"
       << oloc.x1() << " "
       << oloc.x2() << " "
       << oloc.x3() << "\n\n\n";
  *out << "# name: AxesX1\n"                    // Axis X1
       << "# type: matrix\n"                    //
       << "# rows: 1\n"
       << "# columns: 3\n"
       << oax1.x1() << " "
       << oax1.x2() << " "
       << oax1.x3() << "\n\n\n";
  *out << "# name: AxesX2\n"                    // Axis X2
       << "# type: matrix\n"                    //
       << "# rows: 1\n"
       << "# columns: 3\n"
       << oax2.x1() << " "
       << oax2.x2() << " "
       << oax2.x3() << "\n\n\n";
  *out << "# name: AxesX3\n"                    // Axis X3
       << "# type: matrix\n"                    //
       << "# rows: 1\n"
       << "# columns: 3\n"
       << oax3.x1() << " "
       << oax3.x2() << " "
       << oax3.x3() << "\n\n\n";
  *out << "# name: AxesDesc\n"                  // Axes Description
       << "# type: string\n"                    //
       << "# elements: 1\n"
       << "# length: " << mAxesDesc.length()
       << "\n" << mAxesDesc << "\n\n\n";
  *out << "# name: GatherRadius\n"              // Gather Radii
       << "# type: matrix\n"                    //
       << "# rows: 2\n"
       << "# columns: 2\n"
       << mRadiusI[RAY_P] << " " << mRadiusO[RAY_P] << "\n" 
       << mRadiusI[RAY_S] << " " << mRadiusO[RAY_S] << "\n"
       << "\n\n";
  *out << "# name: TimeWindow\n"                // Time Window
       << "# type: matrix\n"                    //
       << "# rows: 1\n"
       << "# columns: 2\n"
       << "0 "
       << cmTimePerBin*cmNumBins << " \n\n\n";
  *out << "# name: NumBins\n"                   // Number of Bins per Series
       << "# type: scalar\n"                    //
       << cmNumBins << "\n\n\n";
  *out << "# name: EventLoc\n"                  // Event Location
       << "# type: matrix\n"                    //
       << "# rows: 1\n"
       << "# columns: 3\n"
       << eloc.x1() << " "
       << eloc.x2() << " "
       << eloc.x3() << "\n\n\n";
  *out << "# name: EventMT\n"                   // Event Moment Tensor
       << "# type: matrix\n"                    //
       << "# rows: 3\n"
       << "# columns: 3\n"
       << cmEventMT.xx() << " " 
          << cmEventMT.xy() << " " 
          << cmEventMT.xz() << " \n" 
       << cmEventMT.yx() << " " 
          << cmEventMT.yy() << " " 
          << cmEventMT.yz() << " \n" 
       << cmEventMT.zx() << " " 
          << cmEventMT.zy() << " " 
          << cmEventMT.zz() << " \n" 
       << "\n" << "\n";
  *out << "# name: Frequency\n"                 // Frequency
       << "# type: scalar\n"                    //
       << MediumCell::GetFrequencyHertz() << " \n\n\n";
  *out << "# name: TraceXYZ\n"                  // XYZ Trace
       << "# type: matrix\n"                    //
       << "# rows:" << cmNumBins << "\n"
       << "# columns: 3 \n";
       for (unsigned int i = 0; i < cmNumBins; i++) {
         *out << "  "
              << setw(14) << mTimeBins[i].mEnergyAxes[AXIS_X] << " "
              << setw(14) << mTimeBins[i].mEnergyAxes[AXIS_Y] << " "
              << setw(14) << mTimeBins[i].mEnergyAxes[AXIS_Z]
              << " \n"; 
       }
  *out << "# name: TracePS\n"                   // PS Trace
       << "# type: matrix\n"                    //
       << "# rows:" << cmNumBins << "\n"
       << "# columns: 2\n";
       for (unsigned int i = 0; i < cmNumBins; i++) {
         *out << " "
              << mTimeBins[i].mEnergyByType[RAY_P] << " "
              << mTimeBins[i].mEnergyByType[RAY_S] 
              << " \n"; 
       }

  *out << "# name: CountPS\n"                   // Phonon-Count Trace
       << "# type: matrix\n"                    //
       << "# rows:" << cmNumBins << "\n"
       << "# columns: 2\n";
       for (unsigned int i = 0; i < cmNumBins; i++) {
         *out << " "
              << mTimeBins[i].mCountByType[RAY_P] << " "
              << mTimeBins[i].mCountByType[RAY_S]
              << " \n";
       }

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  DataReporter                                        ****
// ****                                                              ****
//

//////
// METHOD:   DataReporter :: SuppressAllReports()
//
void DataReporter::SuppressAllReports(bool suppress) {
   
  bool out = false;
  if (suppress == false) {
    out = true;
  }
  
  mbReportGenerate = out;
  mbReportCollect  = out;
  mbReportReflect  = out;
  mbReportTransfer = out;
  mbReportScatter  = out;
  mbReportLost     = out;
  mbReportTimeout  = out;

}


//////
// METHOD:   DataReporter :: SetOutputDirectory()
//
void DataReporter::SetOutputDirectory(std::string outdir) {

  msOutDir = outdir;

  // TODO: verify directory exists first

  // TODO: Currently, there is nothing to enforce sequence when
  // setting output filenames and output directories.  If
  // --report-file precedes --output-dir on the command line, the user
  // may be surprised to find his report file in the cwd instead of
  // the requested directory. Need some infrastructure development to
  // make that work.

}


//////
// METHOD:   DataReporter :: SetReportsFile()
//
//   If called, opens file named in fname with the mofsReports member,
//   and sets mposReports to point to it (instead of to stdout as
//   initialized in constructor).  If this member is never called,
//   then any reports generated go to stdout by default.
//
void DataReporter::SetReportsFile(std::string fname) {

  std::stringstream ss;

  if (msOutDir.size() > 0) {
    ss << msOutDir << "/";
  }
  ss << fname;
  mofsReports.open(ss.str().c_str());
  // TODO:  Check for and handle errors opening file
  mposReports = &mofsReports;
 
}


//////
// METHOD:   output_phonon_dataline()               [Access: Private]
//
//   Writes out the current state of a phonon object out as an ascii
//   dataline.
//
inline void DataReporter
::output_phonon_dataline(std::ostream * out, 
                         const std::string & IDstr, const Phonon & phon) {
  using std::setw;
  const char * RT[2] = {"P", "S"};
  EarthCoords::Generic outloc = ECS.OutConvert(phon.mLoc);

  *out << IDstr
       << setw( 6) << phon.mSID << " "
       << setw( 2) << RT[phon.mType] << " "
       << " ttpl:( "
       << setw(10) << phon.mTimeAlive << " "
       << setw(10) << phon.mPathLength << " "
       << ")   xyz:( " 
       << setw(12) << outloc.x1() << " "
       << setw(12) << outloc.x2() << " "
       << setw(12) << outloc.x3() << " "
       << ")   thph:( " 
       << setw(12) << phon.mDir.Theta() << " "
       << setw(12) << phon.mDir.Phi() << " )  a:( "
       << setw(11) << phon.mAmplitude << " )";

  if (true) {   // Optionally include address of cell on the
                // output line.
                // TODO: Make conditional on user debug-level choice
       *out << "  cell: " << phon.mpCell;
  }

  // Endline
  *out << std::endl;

}


//////
// METHOD:  ReportNewEventPhonon()
//
void DataReporter::ReportNewEventPhonon(const Phonon & phon) {
  if (mbReportGenerate) {
    output_phonon_dataline(mposReports, mIDGenerate, phon);
  }
}


//////
// METHOD:  ReportPhononTimeout()
//
void DataReporter::ReportPhononTimeout(const Phonon & phon) {
  if (mbReportTimeout) {
    output_phonon_dataline(mposReports, mIDTimeout, phon);
  }
}


//////
// METHOD:  ReportScatterEvent()
//
void DataReporter::ReportScatterEvent(const Phonon & phon) {
  if (mbReportScatter) {
    output_phonon_dataline(mposReports, mIDScatter, phon);
  }
}


//////
// METHOD:  ReportPhononCollected()
//
void DataReporter::ReportPhononCollected(const Phonon & phon) {

  // ::::
  // ::     Output immediate report of
  // ::     Phonon hitting surface:
  //
  if (mbReportCollect) {
    output_phonon_dataline(mposReports, 
                           mIDCollect, 
                           phon);
  }//
  // ::::
  // ::     Query each Seismometer for
  // ::     Phonon Detection:
  //
  for (unsigned i = 0; i < mSeismometers.size(); i++) {
    Seismometer & Seis = *mSeismometers[i];
    bool Stop = Seis.CatchPhonon(phon);
    if (Stop) {   // (Search until
      break;      //    we hit one.)
    }             //
  }

}


//////
// METHOD:  ReportReflection()
//
void DataReporter::ReportReflection(const Phonon & phon) {
  if (mbReportReflect) {
    output_phonon_dataline(mposReports, mIDReflect, phon);
  }
}

//////
// METHOD:  ReportCellToCell()
//
void DataReporter::ReportCellToCell(const Phonon & phon) {
  if (mbReportTransfer) {
    output_phonon_dataline(mposReports, mIDTransfer, phon);
  }
}

//////
// METHOD:  ReportLostPhonon()
//
void DataReporter::ReportLostPhonon(const Phonon & phon) {
  if (mbReportLost) {
    output_phonon_dataline(mposReports, mIDLost, phon);
  }
}

//////
// METHOD:  OutputPostSimSummary()
//
void DataReporter::OutputPostSimSummary() {
  std::cerr << "Printing Post-Sim Summary: \n";

  // ::::::
  // :: Print summaries and traces of all Seismometers:
  // :

  std::ostream * out = mpOSSeisTrace;

  for (unsigned int i = 0; i < mSeismometers.size(); i++) {

    Seismometer & Seis = (*mSeismometers[i]);   // Convenience
    char prev_fill = out->fill(' ');  // Save previous fill char

    *out << mIDSeisTrace 
         << "//" << std::setw(54) << "//\n";
    *out << mIDSeisTrace 
         << "------------|   Seismometer Number: "
         << std::setw(3) << std::setfill('0') << i << std::setfill(' ')
         << "   |------------\n";
    *out << mIDSeisTrace 
         << "//" << std::setw(54) << "//\n";

    out->fill(prev_fill);             // Restore fill character

    Seis.OutputDescription(out, mIDSeisTrace);
    Seis.OutputTrace(out,mIDSeisTrace);

    if (false) {  // Also output into individual files
      // TEMP-CODE: Output seismic traces to hard-coded filenames.
      // TODO: Re-write so that file output and the file names are
      //       user-choosable options.
      std::stringstream ss;
      ss << "seis_" 
         << std::setw(3) << std::setfill('0') << i << std::setfill(' ')
         << "_asc.dat";
      std::ofstream seisfile;
      seisfile.open(ss.str().c_str());
      Seis.OutputDescription(&seisfile, mIDSeisTrace);
      Seis.OutputTrace(&seisfile,mIDSeisTrace);
      seisfile.close();
    }

    if (true) {  // Also output into files in GNU/Octave format
      // TEMP-CODE: Hard-coded filenames
      // TODO: Use a user-chooseable filename pattern
      std::stringstream ss;
      if (msOutDir.size() > 0) {
        ss << msOutDir << "/";
      }
      ss << "seis_" 
         << std::setw(3) << std::setfill('0') << i << std::setfill(' ')
         << ".octv";
      std::ofstream seisfile;
      seisfile.open(ss.str().c_str());
      Seis.OutputOctaveText(&seisfile);
      seisfile.close();
    }

  }

}
