// sources.hpp
//
// This file defines a base class for a whole hierarchy of classes to
// represent seismic sources. The basic "family" of these sources will
// divide principally over those sources that are "events" in the
// sense of an earthquake or explosion event, and those that are
// "scatterers" in the sense of treating a scattering event as a
// source of re-radiated phonons.
//
// The principle deriviative files to this file are: events.hpp, and
// scatterers.hpp.
//
// Generally, a code module that needs "sources" will not include this
// file directly, but will include either (or both) of events.hpp or
// scatterers.hpp.  Those files, in turn, include this file.
//
//
#ifndef SOURCES_H_
#define SOURCES_H_
//
#include <iostream> // TODO: drop after prototyping ends
#include <vector>
#include "geom.hpp"
#include "raytype.hpp"
#include "probability.hpp"


//////
// CLASSES: Forward Declarations
//
//   FROM OTHER HEADERS:  (Referenced here by pointer only -
//                         no need to include full header.)

class Phonon;   /* Defined in phonon.hpp */


//////
// TYPEDEFS:
//

//////
// CLASSES: Definitions
//

//////
// CLASS:   PhononSource
///@brief
///
///  Generic base class for Phonon sprayers.
///
///  A very generic base class of a Phonon Source.  Doesn't even include
///  the probability function.  (because at this point the dimensionality
///  of the funciton isn't known.  Could be single-type, or P-S, or
///  scattering PP, SS, PS, SP, etc.)
///
///  ## FEATURES:
///
///  * Static member to maintain take-off angle list.
///
class PhononSource {
protected:

  // :::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Internal Typedefs  (PhononSource Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::

  /*
    // DEPRECATED - REMOVE
  typedef Real                            Probability;
  typedef std::vector<Probability>        ProbabilityArray;
  typedef std::vector<ProbabilityArray>   ProbabilityMatrix;
  */

protected:

  // ::::::::::::::::::::::::::::::::::::::::::::
  // ::: Static Members  (PhononSource Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::

  static int nTOA;           ///< Number of Take-Off-Angles in the TOA set.
  static S2::S2Set * pTOA;   ///< Pointer to class-level take-off-angle
                             ///  array. Defines set of take-off angles to be
                             ///  used by all PhononSource objects.  Needs to
                             ///  be initialized by Set_TOA_Array() at some
                             ///  point prior to construction of first object.

public:

  // ::::::::::::::::::::::::::::::::::::::::::::
  // ::: Static Methods  (PhononSource Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::

  static void Set_TOA_Array(S2::S2Set * toa) {
    pTOA = toa;
    nTOA = pTOA->size();
  }


protected:

  // :::::::::::::::::::::::::::::::::::::::::
  // ::: Member Data  (PhononSource Class) :::
  // :::::::::::::::::::::::::::::::::::::::::

  // PROBABILITIES: 
  //
  // A Probability Distribution (ProbDist object) is an array of
  // values corresponding to the directions in the TOA array that
  // encode the probability of a phonon taking flight at that
  // particular Take-Off Angle.  The PhononSource class maintans a set
  // (vector) of these distributions, corresponding to a plurality of
  // possible output conditions.  (E.g., an event-source generated
  // phonon might begin life as either a P or as an S phonon, and the
  // distributions would be different for each case.)
  //
  // In addition to the set of probability distributions, an
  // additional set of probability arrays is maintained, called "Whole
  // Probabilities" (they are still "distributions," but I use the
  // terminology "arrays" here because they are defined over a set of
  // conditions rather than over the range of ToA's).  The "whole"
  // probabilities define the probability of a phonon assuming each of
  // the states described by the separate Probability Distributions.
  // (E.g., for an event source, the WholeProbs arrays contain
  // elements for the "whole" probabilities of generating a P phonon
  // vs. an S phonon, etc.)  The WholeProbs array is vectorised (ie,
  // we maintain a set of WholeProbs arrays, not just a single array)
  // because the probabilities of generation in each condition might
  // be dependent on some input state.  (This is not the case for
  // event-source generation, so the vector is length one, but for
  // scattering sources, the available output condtions depends on
  // input condition.  Ie, a phonon arriving as a P phonon has a
  // different set of output possibilities compared with one arriving
  // as an S phonon.)

  std::vector<ProbDist>  mPDists;       //!< The Distributions. These
                                        //!< encode the "shape," if you
                                        //!< will, of the radiation/spray
                                        //!< patterns.

  /** Top-tier probability distribution. This encodes the relative
      probability of each "shape" as defined in the mPDists
      distributions. */
  std::vector<ProbDist>  mWholeProbs;   // The Whole Probs
                                        // Array(s). These encode the
                                        // relative probabilities of
                                        // each "shape" defined in the
                                        // mPDists distributions.

  // The process of "spraying" a phonon is then a two-step process
  // where we first randomly choose which "shape," or mPDist
  // Distribution, to use.  This first choice is made based on the
  // relative probabilities contained in the mWholeProbs arrays.  And
  // then the next step is to randomly choose a Take-Off Angle, based
  // on the probabilities encoded in the chosen mPDist distribution.
  // This will give us both an output state (are we P? are we S?) and
  // an output direction, and thus a Phonon can be generated and
  // returned with those properties set.


public:

  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (PhononSource Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::

  PhononSource(int nraytypes_in, int nraytypes_out);

  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Produce-Something Methods  (PhononSource) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  Phonon GenerateRandomPhonon(raytype inray);

  // ::::::::::::::::::::::::::::::::::::::::::::
  // ::: Output Methods  (PhononSource Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::

  void output_wholespace_prob_matrix();
  void output_differential_probabilities(raytype inray, raytype outray, 
                                         const char* symbol = "c",
                                         Real scalefactor = 0.2);
  void output_random_rayset(int nrays, raytype inray = RAY_NA);


};


///
#endif //#ifndef SOURCES_H_
//
