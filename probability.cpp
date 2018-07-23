// probability.cpp
//
#include <iostream>
#include <cstdlib>      /* rand(), srand(), RAND_MAX, exit() */
#include "probability.hpp"
using namespace std;

//////
// METHOD:   ProbDist :: Integrate()
//
//   Switches internal representation of the probability distribution
//   from differential form to cumulative form.  This is done by
//   integrating (summing) the relative probabilities currently in the
//   array over the index of the array.  Turns the contents of the
//   array from an arbitrary non-negative function to a
//   monotonic-or-flatly increasing non-negative function whose final
//   value is the sum-total probability (which may be non-unitary as
//   no assumption of normalization is made).  In this form, the array
//   can be used for probabalistic random choosing of an index.
//
void ProbDist::Integrate() {

  if (mbIntegrated || !mbAllocated) {   // Sanity Check 
    AckSpit();          // Barf if we've already integrated or if
                        // nothing to integrate.
  }

  Count numel = mDist.size();
  for (Index i = 1; i < numel; i++) {
    mDist[i] += mDist[i-1];
  }

  mbIntegrated = true;

}//
//


//////
// METHOD:   ProbDist :: GetMagnitude()
//
//   Return magnitude of distribution (the non-normalized sum
//   total probability).
//
Real ProbDist::GetMagnitude() {

  if (!mbIntegrated) {Integrate();}
  return mDist.back();  // Return last element of cumulative sum

}//
//


//////
// METHOD:   ProbDist :: GetDiffProb()
//
//   Returns a (normalized) differential probability for a given
//   index.
//
//   Calculates the differential probability by taking the difference
//   between value at mDist[idx] and mDist[idx-1] and dividing by the
//   magnitude of the distribution.  Note that this manner of
//   calculation assumes that the distribution is in integrated
//   representation. This assumption is safe because when we request
//   the magnitude it will force a conversion to integrated form if
//   that conversion has not already happened.
//
Real ProbDist::GetDiffProb(Index idx) {

  Real mag = GetMagnitude();  /* (Guarantees Integrated Form) */
  Real prev = 0;
  Real diff;

  if (idx>0) {prev = mDist[idx-1];}

  diff = mDist[idx] - prev;

  if (mag==0) {
    return 0;
  } else {
    return (diff/mag);
  }

}//
//


//////
// METHOD:   ProbDist :: GetRandomIndex()
//
//   Return a randomly-generated index into the distribution.
//
//   RANDOMNESS: Uses rand(), srand(), and RAND_MAX from stdlib.h,
//   necessitating #include <cstdlib>.  TODO: Consider looking into
//   C++'s random library from #include <random> as a possible
//   alternative.  Particularly as regards giving better
//   fineness. (b/c even the current version of g++ on a 64-bit intel
//   machine still only gives 31 bits of randomness with a call to
//   rand(); RAND_MAX is 2^31-1)
//
//   RANDOM SEED: It is presumed that the random number generator has
//   been properly seeded before this function is called.
//
Index ProbDist::GetRandomIndex() {

  if (!mbIntegrated) {Integrate();}  // Needs integral representation

  Index k1 = 0;                 // Initial lower bound for search
  Index k2 = mDist.size() - 1;  // Initial upper bound
  Index k;
  Real r = mDist[k2] * ((Real) rand() / RAND_MAX); 
             // A value between 0.0 and the final (and presumably
             // maximum) value of mDist.

  while (k1 != k2) {
    k = (k1 + k2) >> 1; // Shift-right one bit. Equivalent of
                        // divide-by-two, but no ambiguity about
                        // integer rounding. (Rounds down.)
    if (r <= mDist[k]) {
      k2 = k;           // Then exclude indices greater than k from
                        // the search.
    } else {
      k1 = k+1;         // Otherwise exclude indices less-than or
                        // equal-to k from the search.
    }
  }
  return k2;

}//
//


//////
// METHOD:   ProbDist :: AckSpit()
//
//   Raise a complaint if class used improperly in some way.
//
void ProbDist::AckSpit() {

  cerr << "ERROR: Improper use of ProbDist class.\n";

  exit(1);      // TODO: maybe raise an exception instead.

}
