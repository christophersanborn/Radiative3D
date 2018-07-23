// scatterers.hpp
//
// This file develops a set of classes that represent seismic
// scattering events as "sources" of re-radiated elastodynamic
// phonons.
//
// The classes developed here derive from the base class developed in
// sources.hpp.  Including this file in your code module will
// automatically include sources.hpp as well.
//
//
#ifndef SCATTERERS_H_
#define SCATTERERS_H_
//
#include "sources.hpp"
#include "scatparams.hpp"

//////
// CLASSES: Forward Declarations
//

//////
// TYPEDEFS:
//

//////
// CLASSES: Definitions
//

//////
// CLASS:  Scatterer -FROM- PhononSource
///@brief
///
///   This class represents "Scatterers," as treated theoretically in the
///   manner described by Sato and Fehler.  The Scatterer object derives
///   from, and thus is a kind of, the PhononSource object.  Ie, the
///   Scatterer object is a "source" of scattered phonons.
///
/// ## MEMORY FOOTPRINT
///
///   This class, like any class derived from PhononSource, could be
///   considered somewhat "heavy".  It maintains a large array of
///   probabilities correllating to each possible takeoff-angle, and the
///   Scatterer class in particular also maintains a similar array of
///   polarization values.  Because of this, care and consideration should
///   be given to reducing the total number of Scatterer objects that need
///   to be allocated.  In a model with many volume cells, giving each one
///   its own scatterer might not be a good idea.
///
/// ## ALLOCATION PROCEDURE
///
///   To ameliorate this situation, this class proveds a static method
///   called GetScattererMatchingParams(), the purpose of which is to
///   maintain a list of allocated Scatterers, and to return a pointer to
///   an existing Scatterer if one that sufficiently matches the requested
///   parameters already exists. Otherwise, it constructs and records a new
///   scatterer and returns a pointer to that one.  The model-building code
///   should use this member rather than the 'new' keyword to link volume
///   cells to Scatterers.
///
/// ## CACHEING PROCEDURE
///
///   To further ameliorate the memory demands of the scattering class,
///   the class will maintain the ability to flush the Probability
///   arrays if the net memory demands become too large.  Each
///   scattering object will keep a use-count so that more frequently
///   used Probability arrays can be preferentially kept compared to
///   less frequently used arrays.  In essence, the most frequently
///   used arrays will be kept in a cache, whereas the less frequently
///   used arrays will be re-generated as needed.
///
class Scatterer : public PhononSource {
protected:

  // ::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Enums and Constants  (Scatterer Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::

  enum in_types_e {IN_P, IN_S, NUM_INTYPES};
                                // IN_P and IN_S agree with raytype
                                // enum and can be used interchangably
                                // where P/S raytypes needed.

  enum out_types_e {GPP, GPS, GSP, GSS, NUM_OUTTYPES};
                                // Four basic conversions denote the
                                // four scattering distributions
                                // (shapes) used. (Note that Sato
                                // distinguishes between G_S->SH and
                                // G_S->SV, but we keep just a net
                                // distribution for S->S and
                                // distinguish polarization via a
                                // separate m_spol array.)


  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Class-Static Infrastructure  (Scatterer Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::


  static Scatterer  * cm_ll_first;  // Points to first element of the
                                    // linked list of Scatterer objects.

  static Scatterer   * cm_ll_last;  // Points to last element of the
                                    // linked list of Scatterer objects.

  static Count cm_ll_count;         // Keeps tally of how many elements 
                                    // are on the linked list.

  static bool cm_MFPOverride_b;         // True to overide mean-free paths
  static Real cm_MFPOverrides[RAY_NBT]; // with explicit values. (Applies
                                        // to all Scatterer objects.)

  static bool cm_NoDeflect_b;       // True to suppress deflection (phonon
                                    // direction, raytype, polarization
                                    // will not be changed by scattering
                                    // event). (Applies to all objects.)

                // NOTE: if MFP and deflection are BOTH overridden, then
                // ScatterParams have no remaining effect, and no more than
                // one object need be created. GetScattererMatchingParams()
                // should detect this and squash par to limit list to one
                // scatterer.

public:
  ;
  static Scatterer * 
  GetScattererMatchingParams(ScatterParams par);
                      // Retrieve or create a Scatterer object
                      // matching to within tolerances the provided
                      // parameters.

  static Count GetCollectionSize() {return cm_ll_count;}
                      // How many Scatterers are in the collection of
                      // Scatterer objects known to the class.

  static void OverrideMFP(Real mfp_P, Real mfp_S) {
    cm_MFPOverride_b = true;
    cm_MFPOverrides[RAY_P] = mfp_P;
    cm_MFPOverrides[RAY_S] = mfp_S;
  }

  static void SetNoDeflect() {
    cm_NoDeflect_b = true;
  }


protected:

  // :::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Private Member Variables  (Scatterer Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::

  Scatterer        * mpllPrev;  // Linked-list pointers
  Scatterer        * mpllNext;  //

  Real  mMeanFreeP[RAY_NBT];    // Mean Free Path before scattering.
                                // Indexed by raytype. 

  Real    mDipoles[RAY_NBT];    // Dipole moments of the distributions.
  Real  mQuadpoles[RAY_NBT];    // Quadrupole moments.
                                //
                                // These values are computed because
                                // they make an interesting way to
                                // quantitatively summarize the
                                // scattering shapes, but they are not
                                // needed or used in the actual
                                // simulation. They are just computed
                                // so that they can be reported to the
                                // user if desired.

  ScatterParams     mParams;    // Parameters used to generate the
                                // probability distributions.
                                // Retained incase the distribution
                                // needs to be recomputed after a
                                // cache flush. (The ability to dump
                                // the distribution for memory
                                // conservation is a planned future
                                // feature.)

  std::vector<Real>  m_spol;    // Polarization angles for S->S
                                // conversions


public:
  
  // :::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Scatterer Class) :::
  // :::::::::::::::::::::::::::::::::::::::

  Scatterer(ScatterParams par);


protected:
;
  // :::::::::::::::::::::::::::::::::::::::::::
  // ::: Helper Functions  (Scatterer Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::

  void PopulateProbDists(ScatterParams par);
  void PopulateWholeProbs();
  void ComputeMFPs();
  void ComputeDipoles();
  void ComputeQuads();


public:
;
  // :::::::::::::::::::::::::::::::::::::::::::::
  // ::: Manipulate Methods  (Scatterer Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::

  void SetMeanFreePathsPS(Real, Real);
  //            Though the MFPs are computed automatically by the
  //            constructor, this method can be used to change them
  //            after the fact.  This is mainly for diagnostic and
  //            development/debugging purposes.


  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Compute Methods  (Scatterer Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::

  Real GetRandomPathLength(raytype intype);
  Phonon GetRandomScatteredRelativePhonon(raytype intype);


  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Testing Methods  (Scatterer Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::

  void test_random_rayset(raytype intype, // Output to stdout a
                          int count);     // randomly generated set
                                          // of rays; For testing
                                          // scatter patterns.

  static void PrintAllScatteringStats();  // Walks the linked-list of
                                          // all scatterers and
                                          // outputs their stats.

  void PrintStats();    // Dump stats for particular scatterer.


}; // class Scatterer
////

///
#endif //#ifndef SCATTERERS_H_
//
