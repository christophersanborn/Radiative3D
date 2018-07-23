// elastic.cpp
//
#include "elastic.hpp"

namespace Elastic {

// Helper Functions:
inline Real Q::HelperL(Velocity vpvs) {
  Real L;
  L = vpvs.Vs()/vpvs.Vp();  // beta/alpha
  L *= L;                   // (beta/alpha)^2
  L *= (4./3.);             // (4/3)(beta/alpha)^2
  return (L);
}
  

Real Q::Qp(Velocity vpvs) const {
  Real Qval;
  if (mUnknown==QUNK_QP) {
    Real L = HelperL(vpvs);
    Real L1 = 1. - L;
    if (L==0.) {        // Prevent a NaN if Vs = 0 and Qs = 0.
      Qval = 0;         // (Assumes Vs->0 dominates in the limit.)
    } else {            // This allows reasonable Qp values if user
      Qval = L/mQs;     // sets Qs=0 for water. 
    }                   //
    Qval += L1/mQk;
    Qval = 1./Qval;
  } else {
    Qval = mQp;
  }
  return (Qval);
}


Real Q::Qs(Velocity vpvs) const {
  Real Qval;
  if (mUnknown==QUNK_QS) {
    Real L = HelperL(vpvs);
    Real L1 = 1. - L;
    Qval =  1./mQp;     // Can return nan if mQp
    Qval -= L1/mQk;     // and mQk are both infinite
    Qval =  L/Qval;     // and S velocity is zero.
  } else {
    Qval = mQs;
  }
  return (Qval);
}


Real Q::Qk(Velocity vpvs) const {
  Real Qval;
  if (mUnknown==QUNK_QK) {
    Real L = HelperL(vpvs);
    Real L1 = 1. - L;
    Qval =  1./mQp;     // Can return nan if mQp and
    Qval -= L/mQs;      // mQs are both infinite and
    Qval =  L1/Qval;    // if vs = sqrt(3/4) * vp,
  } else {              // (or about 0.866 vp).
    Qval = mQk;
  }
  return (Qval);
}


} // End namespace Elastic
///
