// scatterers.cpp
//
#include <cmath>        /* log(), sqrt() */
#include <iomanip>      /* setw() */
#include <cstdlib>      /* rand(), RAND_MAX */
#include "scatterers.hpp"
#include "phonons.hpp"

                        /* Assume srand() has been called elsewhere */

using namespace std;

// in this file:
// CLASS IMPLEMENTATIONS FOR:
//
//   o  Class Scatterer
//


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  Scatterer                                           ****
// ****                                                              ****
//

//
// Static Member Initialization:  (Scatterer Class)
//

Scatterer * Scatterer::cm_ll_first = 0; // Linked list pointers
Scatterer * Scatterer::cm_ll_last  = 0; //
Count Scatterer::cm_ll_count = 0;  // Linked list size
bool Scatterer::cm_MFPOverride_b = false;
Real Scatterer::cm_MFPOverrides[] = {0,0};
bool Scatterer::cm_NoDeflect_b = false;

//////
// STATIC METHOD:  Scatterer :: GetScattererMatchingParams()
//
//   Return pointer to an existing Scatterer if one can be found that
//   suitable matches the provided ScatterParams object, otherwise
//   create a new Scatterer and return pointer if no such match can be
//   found.
//
Scatterer *
Scatterer::GetScattererMatchingParams(ScatterParams par) {

  if (cm_MFPOverride_b && cm_NoDeflect_b) { // Check for overrides
    par = ScatterParams(Elastic::HSneak(1.0,0.0,1.0,1.0),1.0,1.0);
    // Squash parameters to dummy value, so that one and only one
    // object will ever be built.
  } //

  Scatterer * result = 0;
  Scatterer * fromlist = cm_ll_first;

  while (fromlist != 0) { // First scan linked list of existing
                          // Scatterers for possible match:
                          //
    Real diff = par.CompareRoughly(fromlist->mParams);
                          // A very rough comparison metric...
    if (diff <= 0) {      // A looser cuttoff may be in order...
      result = fromlist;
      break;
    } else {              // else keep walking the list
      fromlist = fromlist->mpllNext;
    }
  }

  if (result != 0) {      // If we found a match, return it
    return result;        //
  }                       // else continue as follows:

  result = new Scatterer(par);      // Allocate new scatterer
  cm_ll_count++;                    //

  if (cm_ll_first == 0) {           // And append to linked list...
    cm_ll_first = result;           // (set pointer to first if needed)
  }

  result->mpllPrev = cm_ll_last;    // link new object to previous oject;
  result->mpllNext = 0;             // mark new object as last in list;
  cm_ll_last = result;              // update pointer to last;

  if (result->mpllPrev != 0) {      // and link previous object (if any)
    result->mpllPrev->mpllNext = result;  // to new object;
  }

  return result;

}


//////
// CONSTRUCTOR:  Scatterer()
//
Scatterer::Scatterer(ScatterParams par) :
  PhononSource(NUM_INTYPES, NUM_OUTTYPES),  
                        // Init base class for two input raytypes and
                        // four output raytypes.
  mParams(par)          // Record params for later use if needed
{

  PopulateProbDists(par);   // Populate probability distributions
  PopulateWholeProbs();     //

  if (!cm_MFPOverride_b) {  // Populate Mean Free Path values
    ComputeMFPs();          //
  } else {
    mMeanFreeP[RAY_P] = cm_MFPOverrides[RAY_P];   // Override values
    mMeanFreeP[RAY_S] = cm_MFPOverrides[RAY_S];   // if requested
  }

  if (!cm_NoDeflect_b) {    // Populate Dipole Moments
    ComputeDipoles();       //
  } else {
    mDipoles[RAY_P] = 1.0;  // Note: these are informational values, over-
    mDipoles[RAY_S] = 1.0;  // riding them does not prevent deflection. A
  }                         // check in GRSPh() does that.

  //
}//
//


//////
// METHOD:  Scatterer :: PopulateProbDists()         [Helper Function]
//
//   Computes and populates the mPDists members.  Called by the
//   constructor, and whenever the (memory-hungry) mPDists member
//   might need to be re-computed (e.g. after a cache-miss if the
//   distribution has been flushed).
//
void Scatterer::PopulateProbDists(ScatterParams par) {

  S2::S2Set & toa = (*pTOA);    // (Convenience alias)

  m_spol.clear();       // spol is an array of polarization angles for
  m_spol.resize(nTOA);  // the S->S conversionevents.  Technically,
                        // there are probabilities of going to either
                        // SH (spol=90) or SV (spol=0), but we only
                        // treat the probability of going to S
                        // generally, and code a polarization angle
                        // that represents an appropriate linear
                        // combination of 0 and 90, (based on the
                        // individual SV/SH probabilities, which are
                        // known to the GSATO function.)
  //
  // RADIATION PATTERNS:
  //   (Populate PF arrays)
  Real gpp, gps, gsp, gss, spolv;
  for (int k = 0; k < nTOA; k++) {
    par.GSATO(toa[k], gpp, gps, gsp, gss, spolv); // Get g-values and
    mPDists[GPP].SetRelativeProb(k, gpp);         // spol value
    mPDists[GPS].SetRelativeProb(k, gps);
    mPDists[GSP].SetRelativeProb(k, gsp);
    mPDists[GSS].SetRelativeProb(k, gss);
    m_spol[k]  = spolv;
  }

}//
//


//////
// METHOD:  Scatterer :: PopulateWholeProbs()        [Helper Function]
//
//   Computes and populates the WholeProbs probability arrays.
//
//   Assumes mPDists[] has already been populated.
//
void Scatterer::PopulateWholeProbs() {

  mWholeProbs[IN_P].SetRelativeProb(GPP, mPDists[GPP].GetMagnitude());
  mWholeProbs[IN_P].SetRelativeProb(GPS, mPDists[GPS].GetMagnitude());
  mWholeProbs[IN_P].SetRelativeProb(GSP, 0); // (S->P precluded for IN_P)
  mWholeProbs[IN_P].SetRelativeProb(GSS, 0); // (S->S precluded for IN_P)

  mWholeProbs[IN_S].SetRelativeProb(GPP, 0); // (P->P precluded for IN_S)
  mWholeProbs[IN_S].SetRelativeProb(GPS, 0); // (P->S precluded for IN_S)
  mWholeProbs[IN_S].SetRelativeProb(GSP, mPDists[GSP].GetMagnitude());
  mWholeProbs[IN_S].SetRelativeProb(GSS, mPDists[GSS].GetMagnitude());

}//
//


//////
// METHOD:  Scatterer :: ComputeMFPs()               [Helper Function]
//
//   Computes and stores Mean Free Paths.
//
//   Assumes mWholeProbs[] has already been populated.
//
void Scatterer::ComputeMFPs() {

  //
  // Get the *inverse* MFPs  (Basically the G0 values)
  //
  Real imfp_p = mWholeProbs[IN_P].GetMagnitude();
  Real imfp_s = mWholeProbs[IN_S].GetMagnitude();
  //
  // Compute mean-free-paths:
  //
  //   (The imfp's currently contain the sums of the G-functions.  The
  //   "probability of scattering per unit length" is the "surface
  //   average" of the G-functions. Note that the individual G-values
  //   have NOT been scaled by the size of their area element (this
  //   remains on the TODO: list - for now we make the approximation
  //   that each ToA represents an equal area), thus instead of
  //   dividing by 4Pi to get the surface average, we actually need to
  //   divide by nTOA.)  //TODO: TesselSphere needs to code area
  //   elements.
  //
  imfp_p /= nTOA;  // TODO: This will change when above TODO fixed
  imfp_s /= nTOA;
  mMeanFreeP[IN_P] = 1.0 / imfp_p;
  mMeanFreeP[IN_S] = 1.0 / imfp_s;

}//
//


//////
// METHOD:   Scatterer :: ComputeDipoles()
//
//   Computes the dipole moments of the scattering shapes. This
//   provides a useful way to characterize the scattering shapes on a
//   spectrum between "forward scattering" and "back scattering."
//
//   Returns a real-valued number between -1.0 and 1.0 where numbers
//   approaching 1.0 mean dominantly forward-scattering (minimal path
//   deviation), -1.0 means dominant back-scattering (~180*
//   reversals), and 0.0 means either equal contribution from both
//   forward and backwards or else the dominant scattering direction
//   is side-deflection.
//
//   Computation is achieved by multiplying each differential
//   probability by cos(theta) where theta is the deflection angle
//   (co-lattitude from forward axis).
//
//   Assumes mPDists[] and mWholeProbs[] have already been computed
//   and are available for analysis.
//
void Scatterer::ComputeDipoles() {

  const S2::S2Set & toa = (*pTOA);          // Alias (for convenience)
  Real moments[NUM_OUTTYPES] = {0,0,0,0};   // Four moments

  for (int idx = 0; idx < nTOA; idx++) {
    Real costh = toa[idx].z();  // (Unit vec, so z is cos theta)
    moments[GPP] += costh * mPDists[GPP].GetDiffProb(idx);
    moments[GPS] += costh * mPDists[GPS].GetDiffProb(idx);
    moments[GSP] += costh * mPDists[GSP].GetDiffProb(idx);
    moments[GSS] += costh * mPDists[GSS].GetDiffProb(idx);
  }

  mDipoles[IN_P] = moments[GPP] * mWholeProbs[IN_P].GetDiffProb(GPP)
                 + moments[GPS] * mWholeProbs[IN_P].GetDiffProb(GPS);
  mDipoles[IN_S] = moments[GSP] * mWholeProbs[IN_S].GetDiffProb(GSP)
                 + moments[GSS] * mWholeProbs[IN_S].GetDiffProb(GSS);

}//
//


//////
// METHOD:  Scatterer :: SetMeanFreePathsPS()
//
//   Though the MFPs are computed automatically by the constructor,
//   this method can be used to change them after the fact.  This is
//   mainly for diagnostic and development/debugging purposes, or for
//   creating very specialized and not-necessarily-physical test
//   models. It should not be part of a normal grid-based model
//   building process.
//
//   NOTE: C.f. this function with class static OverrideMFP().  This function
//   overrides the MFP values for a particular object.  OverrideMFP()
//   overrides all objects at construction time.
//
void Scatterer::SetMeanFreePathsPS(Real mfpP, Real mfpS) {
  mMeanFreeP[RAY_P] = mfpP;
  mMeanFreeP[RAY_S] = mfpS;
}//
//


//////
// METHOD:  Scatterer :: GetRandomPathLength()
//
// PURPOSE: Return a randomly-generated path length to a scattering
//          event.  Lengths are generated in such a way that the
//          average, or "mean free path," of many iterations of this
//          function will tend towards the values stored in
//          MeanFreePath[].
//
Real Scatterer::GetRandomPathLength(raytype intype) {

  Real r = ((Real) rand()) / ((Real) RAND_MAX + 1); 
                // A double in the range [0,1)
  r = 1.0 - r;  // A double in the range (0,1]
                // TODO: This is a little ugly. Maybe compute r like
                // we do in rtcoef.cpp instead. (~CJS 20140514)

  return -log(r) * mMeanFreeP[intype];

}


//////
// METHOD:  Scatterer :: GetRandomScatteredRelativePhonon()
//
// PURPOSE: Generate and return a "relative" phonon at random from the
//          scattering probability distributions.  The "relative"
//          phonon can be used to "bend and rotate" the path of the
//          incoming phonon.
//
Phonon Scatterer::GetRandomScatteredRelativePhonon(raytype intype) {
  S2::S2Set & toa = (*pTOA);             // Alias
  raytype out_types[4] = {RAY_P, RAY_S,  // Maps conversion types
                          RAY_P, RAY_S}; // (GPP, GPS, GSP, GSS) to
                                         // output raytypes (RAY_P or
                                         // RAY_S)

  if (cm_NoDeflect_b) {     // First, check if deflection is overridden:
    Phonon nodeflect = Phonon(S2::ThetaPhi(0,0),intype);
    nodeflect.SetPolarization(0); // And return a non-deflected
    return nodeflect;             // result if so.
  }                               // Otherwise, continue:

  // Get output conversion type and corresponding output ray type:
  out_types_e conv = (out_types_e) mWholeProbs[intype].GetRandomIndex();
  raytype ort = out_types[conv];

  // Get take-off-angle index for output ray:
  Index toa_index = mPDists[conv].GetRandomIndex();

  // Ray polariztion depends on conversion type:
  Real pol = 0;
  switch (conv) {
  case GPP : // Polarization meaningless for P-waves; default to zero.
    pol = 0;
    break;
  case GSP :
    pol = 0;
    break;
  case GPS : // This conversion always makes theta-hat polarization.
    pol = 0; // (0 for theta-hat, Pi/2 for phi-hat)
    break;
  case GSS : // This conversion looks up polarization in spol array.
    pol = m_spol[toa_index];
    break;
  default:   // (Shouldn't end up here - this is just a #pragma
    break;   //  for NUM_OUTTYPES - avoids compiler warning.)
  }

  // Generate the Phonon:
  Phonon phon = Phonon(toa[toa_index], ort);
  phon.SetPolarization(pol);

  return phon;

}


//////
// METHOD:  Scatterer :: test_random_rayset()
//
// PURPOSE: Print on stdout a randomly-generated set of rays,
// formatted as take-off angles and symbol codes suitable for GMT
// plotting.  This is for testing/validating of scatter patterns.
//
void Scatterer::test_random_rayset(raytype intype, int count) {
  S2::S2Set & toa = (*pTOA); 
  vector< vector<int> > bins;
  bins.clear();
  bins.resize(3);               // For raytypes RAY_P, RAY_SH, RAY_SV
  for (int k=0; k<3; k++) {
    bins[k].clear();
    bins[k].resize(nTOA,0);
  }
  //
  // FILL BINS:
  //  
  for (int ctr=0; ctr<count; ctr++) {
    Phonon phon = GetRandomScatteredRelativePhonon(intype);
    raytype rt = phon.GetFullRaytype();
    int toa_index = toa.GetBestIndexFromPoint(phon.GetDirection());
    bins[rt][toa_index]++;
  }
  //
  // OUTPUT:
  //
  const char * symbol[3];
  symbol[RAY_P]  = "c";    // GMT circle code
  symbol[RAY_SH] = "-";    // GMT horiz. bar code
  symbol[RAY_SV] = "y";    // GMT vert. bar code
  for (int irt=0; irt<3; irt++) {
    for (int itoa=0; itoa<nTOA; itoa++) {
      int num = bins[irt][itoa];
      if (num != 0) {
        Real val = (Real) num / count;
        val = sqrt(val * 3.0 * nTOA) * 0.2;
        cout << setw(16) << toa[itoa].Lon()
             << setw(16) << toa[itoa].Lat()
             << setw(16) << val << " "
             << symbol[irt] << endl;
      }
    }
  }

}

//////
// METHOD:  Scatterer :: PrintAllScatteringStats()            [static]
//
//   Calls the PrintStats() method on each Scatterer object in the
//   class-maintained linked list of Scatterer objects.
//
void Scatterer::PrintAllScatteringStats() {     // STATIC METHOD

  Scatterer * next = cm_ll_first;   // First object in linked list

  cout << "#  BEGIN SCATTERER DUMP:\n";
  cout << "#  Overrides: "
       << ((cm_MFPOverride_b) ?
           ((cm_NoDeflect_b) ? "MFPs, Deflections" : "MFPs") :
           ((cm_NoDeflect_b) ? "Deflections" : "None"))
       << std::endl;
  cout << "#  "
       << "    nu    eps        a    kappa         el       gam0    "
       << "     address      MFP (P)    MFP (S)      DM (P)    DM (S)\n";
  cout << "#  "
       << "====== ====== ======== ======== ========== ==========    "
       << "============    =========  =========    ========  ========\n";

  while (next != 0) {
    next->PrintStats();
    next = next->mpllNext;
  }

  cout << "#  "
       << "====== ====== ======== ======== ========== ==========    "
       << "============    =========  =========    ========  ========\n";
  cout << "#  END SCATTERERS\n";

}


//////
// METHOD:  Scatterer :: PrintStats()
//
//   Prints to stdout a breakdown of various stats of the scattering
//   object, including the parameters it was constructed on, and
//   various derived quantities, such as mean-free-paths.
//
void Scatterer::PrintStats() {

  const int width1 = 8;
  const int width2 = 12;

  cout << "   "
       << setprecision(6)
       << setw(width1-2) << mParams.GetNu() << " "
       << setw(width1-2) << mParams.GetEps() << " "
       << setw(width1)   << mParams.GetA() << " "
       << setw(width1)   << mParams.GetKappa() << " "
       << setw(width1+2) << mParams.GetL() << " "
       << setw(width1+2) << mParams.GetGam0() << "    "
       << setw(width2)   << this << "    "
       << setw(width1+1) << mMeanFreeP[RAY_P] << "  "
       << setw(width1+1) << mMeanFreeP[RAY_S] << "    "
       << setprecision(4)
       << setw(width1)   << mDipoles[RAY_P] << "  "
       << setw(width1)   << mDipoles[RAY_S] << "\n";


}
