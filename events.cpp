// events.cpp
//
#include <cmath>
#include "events.hpp"
#include "phonons.hpp"
#include "dataout.hpp"

//
// CLASS IMPLEMENTATIONS:

//////
// CLASS:  Sources::ShearDislocation
// 

//////
// CONSTRUCTOR:  ShearDislocation()
//
//   Builds out the probability arrays for P, SH, SV emission for an
//   event source described by a rank-2 moment tensor (includes
//   shear-dislocation, CLVD, explosion type sources), as described by
//   Aki and Richards box 9.10 of Quantitative Seismology, 1980.
//
//   Note that A&R use a reference frame where the coordinate axes
//   x,y,z point in the directions North, East, Down (NED).  There is
//   no guarantee that the Earth model at the location of the source
//   respects this convention.  The code below will blindly compute
//   the probability distributions in the simulation-frame XYZ system
//   without any regard for "which way is up".  If there is a
//   discrepency between the local NED directions and the
//   simulation-frame XYZ coordinates, then it is here assumed that
//   the NED moment tensor provided by the user has already been
//   rotated (perhaps in the model-building code) to compensate.  This
//   can easily be achieved with the aid of the ECS (Earth Coordinate
//   System) libraries.  See Model() constructor for actual
//   implementation.
//
//   The 'Loc' (location) argument is retained by the event source
//   object for the purposes of phonon generation only. It is not
//   needed to get the local NED directions, as it is assumed the
//   relevant transformation has already been applied.
//
ShearDislocation::ShearDislocation(Tensor::Tensor MT, R3::XYZ Loc) :
  // Initialize base class, allocate probability matrices: 
  PhononSource(1,RAY_NUMTYPES),
  // Set up event-source specific class vars:
  mLoc(Loc),
  mAmp_P(1.0),
  mAmp_S(1.0),
  mpCell(0)
{

  //
  // Moment Tensor Elements:
  //

  Real mxx = MT.xx();
  Real myy = MT.yy();
  Real mzz = MT.zz();
  Real mxy = MT.xy();  // Assumes symmetric (xy == yx)
  Real mxz = MT.xz();  // ''
  Real myz = MT.yz();  // ''

  //
  // RADIATION PATTERNS:
  //
  for (int ctr = 0; ctr < nTOA; ctr++) {
    // Take-off angle:
    Real theta = (*pTOA)[ctr].Theta();
    Real az    = (*pTOA)[ctr].Phi();
    //
    // P-WAVE PATTERN:
    // Calculate the probability *amplitude* at this theta,phi:
    //
    Real prob;
    prob = pow(sin(theta),2) 
      * (mxx*pow(cos(az),2) + mxy*sin(2*az) + myy*pow(sin(az),2) - mzz)
      + 2*sin(theta)*cos(theta)
      * (mxz * cos(az)  +  myz * sin(az)) 
      + mzz;
    prob = prob * prob;    // Actual probability is amplitude squared.
    mPDists[RAY_P].SetRelativeProb(ctr, prob);  // Store in the array.
    //
    // SH-WAVE PATTERN:
    // Calculate the probability *amplitude* at this theta,phi:
    //
    prob = sin(theta) * (0.5*sin(2*az)*(myy-mxx) + cos(2*az)*mxy)
                + cos(theta) * (cos(az)*myz - sin(az)*mxz);
    prob = prob * prob;    // Actual probability is amplitude squared.
    mPDists[RAY_SH].SetRelativeProb(ctr, prob); // Store in the array.
    //
    // SV-WAVE PATTERN:
    // Calculate the probability *amplitude* at this theta,phi:
    //
    prob = sin(theta)*cos(theta) 
      * (mxx*pow(cos(az),2) + mxy*sin(2*az) + myy*pow(sin(az),2) - mzz)
      + (1.0 - 2*pow(sin(theta),2)) * (mxz*cos(az) + myz*sin(az));
    prob = prob * prob;    // Actual probability is amplitude squared.
    mPDists[RAY_SV].SetRelativeProb(ctr, prob); // Store in the array.
  }
  //
  // And now the (relative) whole space probabilities:
  //
  mWholeProbs[0].SetRelativeProb(RAY_P, mPDists[RAY_P].GetMagnitude());
  mWholeProbs[0].SetRelativeProb(RAY_SH, mPDists[RAY_SH].GetMagnitude());
  mWholeProbs[0].SetRelativeProb(RAY_SV, mPDists[RAY_SV].GetMagnitude());
  //
}//
//

//////
// METHOD:
///@brief
///
///   Generates a new Phonon at the event source location and
///   probabilistically assigns it an inital direction, polarization,
///   etc.
///
Phonon ShearDislocation::GenerateEventPhonon() {

  Phonon P = GenerateRandomPhonon(RAY_NA);

  P.SetLocation(mLoc);    // Inform phonon of its location
  P.InsertInto(mpCell);   // And in which cell it is

  // (TODO: code to set relative amplitudes.)

  dataout.ReportNewEventPhonon(P);  // Let the whole world know

  return P;

}

