// sources.cpp
//
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <cstdlib>      /* rand(), srand(), RAND_MAX, exit() */
#include "sources.hpp"
#include "phonons.hpp"

using namespace std;


// ************************************************
// *** CLASS IMPLEMENTATION:  Sources::PhononSource
// ************************************************
//

//////
// Static Members:
//
int PhononSource::nTOA = 0;
S2::S2Set * PhononSource::pTOA = NULL;


//////
// Static Methods:
//

//////
// Constructors:
//

//////
// CONSTRUCTOR:  PhononSource(int, int)
///@brief
///
///   Initializes a two-tiered probability map controlling initial
///   phonon state and take-off angle, and, (todo) initial phase and
///   polarization angle.
///
///   Sets up the two Probability Distribution (ProbDist)
///   vectors: #mPDists and #mWholeProbs.  The two parameters define the
///   number of input conditions, and the number of output conditions.
///   These determine the size of the vectors and the length of the
///   WholeProbs arrays.
///
///   The mPDists vector will have nraytypes_out number of elements,
///   and each contianed ProbDist object will have length nTOA.
///
///   The mWholeProbs vector will have nraytypes_in number of elements,
///   and each contained ProbDist object will have length
///   @p 'nraytypes_out'.
///
///   Note that this BASE CLASS constructor constructs and allocates
///   (resizes) the ProbDist vectors, but does NOT populate them with
///   actual probability values.  That is the task of derived classes,
///   which will be specialized to the particular type of source
///   (e.g. event-source, scattering, etc.) being modeled.
///
///   NOTE: It is assumed that the static member nTOA has been properly
///   set prior to construction of the first PhononSource object.  If
///   not... there will be problems.
///
PhononSource::PhononSource(int nraytypes_in,  //!< Num input raytypes
                           int nraytypes_out  //!< Num output raytypes
                          ) :
  mPDists     (nraytypes_out, ProbDist(nTOA)          ),
  mWholeProbs (nraytypes_in,  ProbDist(nraytypes_out) )
{}


//////
// METHOD:   PhononSource :: output_differential_probabilities()
//
//   Used for diagnostic purposes
//
void PhononSource::output_differential_probabilities(raytype inray, 
                                                     raytype outray,  
                                                     const char* symbol,
                                                     Real scalefactor) {
  ProbDist & Dist = mPDists[outray];
  Real rel_prob = mWholeProbs[inray]
                   .GetDiffProb(outray); // Relative scaling
  S2::S2Set & takeoffs = (*pTOA);        // Alias for takeoff-angle array

  for (int ctr = 0; ctr < nTOA; ctr++) {
    Real prob = rel_prob * Dist.GetDiffProb(ctr);
    cout << setw(16) << takeoffs[ctr].Lon() 
         << setw(16) << takeoffs[ctr].Lat() 
         << setw(16) << sqrt(prob * nTOA) * scalefactor
         << "  " << symbol << endl; // Symbol holds a GMT flag marker code.
  }
}


void PhononSource::output_random_rayset(int nrays, raytype inray) {
  S2::S2Set & takeoffs = (*pTOA); // An alias for the takeoff angle array
  vector< vector<int> > counts;
  counts.clear();
  counts.resize(3);
  vector<int> & P_count = counts[RAY_P];
  vector<int> & SH_count = counts[RAY_SH];
  vector<int> & SV_count = counts[RAY_SV];
  P_count.clear();
  SH_count.clear();
  SV_count.clear();
  P_count.resize(nTOA,0);
  SH_count.resize(nTOA,0);
  SV_count.resize(nTOA,0);

  // Cast rays
  for (int ctr = 0; ctr < nrays; ctr++) {
    raytype rt = (raytype) mWholeProbs[inray].GetRandomIndex();
    Index toa_index = mPDists[rt].GetRandomIndex();
    counts[rt][toa_index]++;
  } 

  // Output Rays:
  // P_waves:
  for (int ctr = 0; ctr < nTOA; ctr++) {
    int num = P_count[ctr];
    const char * symbol = "c";
    if (num != 0) {
      Real val = (Real) num / nrays;
      cout << setw(16) << takeoffs[ctr].Lon() 
           << setw(16) << takeoffs[ctr].Lat() 
           << setw(16) << sqrt(val * 3.0 * nTOA) * 0.2
           << "  " << symbol << endl; // Symbol holds a GMT flag marker code.
    }
  }
  // SH_waves:
  for (int ctr = 0; ctr < nTOA; ctr++) {
    int num = SH_count[ctr];
    const char * symbol = "-";
    if (num != 0) {
      Real val = (Real) num / nrays;
      cout << setw(16) << takeoffs[ctr].Lon() 
           << setw(16) << takeoffs[ctr].Lat() 
           << setw(16) << sqrt(val * 3.0 * nTOA) * 0.2
           << "  " << symbol << endl; // Symbol holds a GMT flag marker code.
    }
  }
  // SV_waves:
  for (int ctr = 0; ctr < nTOA; ctr++) {
    int num = SV_count[ctr];
    const char * symbol = "y";
    if (num != 0) {
      Real val = (Real) num / nrays;
      cout << setw(16) << takeoffs[ctr].Lon() 
           << setw(16) << takeoffs[ctr].Lat() 
           << setw(16) << sqrt(val * 3.0 * nTOA) * 0.2
           << "  " << symbol << endl; // Symbol holds a GMT flag marker code.
    }
  }

}


//////
// METHOD: PhononSource::GenerateRandomPhonon()
//
Phonon PhononSource::GenerateRandomPhonon(raytype inray) {
  
  // Determine output raytype:
  raytype rt = (raytype) mWholeProbs[inray].GetRandomIndex();

  // Get a Take-off Angle (TOA) index:
  Index toa_index = mPDists[rt].GetRandomIndex();

  // Resolve TOA index to a ThetaPhi direction:
  S2::ThetaPhi dir = (*pTOA)[toa_index];

  // Return a Phonon constructed from direction and raytype:
  return Phonon(dir, rt);

}
