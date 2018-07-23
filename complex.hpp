// complex.hpp
//
// This file develops a class to represent complex numbers.  It is
// based (at least at present) on the STL std::complex<T> class
// template, and does not seek to provide any new functionality, but
// rather just a convenient implementation.  In particular, I have
// implemented abs(), arg(), and norm() from the std namespace as
// member functions, which is what I wish the template designers had
// done.  All member functions are defined inline inside the class
// definition, and thus, with compiler optimizations, should hopefully
// not add any significant overhead over using the raw STL template.
// Also, everything needed is in this header file - i.e., there is no
// corresponding .cpp module.
// 
#ifndef COMPLEX_H_
#define COMPLEX_H_
//
#include <complex>
#include "typedefs.hpp"     /* typedef Real */

//////
// CLASS:   ::::  Complex  ::::
//
//   Provides a full-featured representation of complex numbers
//   consisting of a real and imaginary part utilizing the same
//   machine representation as defined by type Real (Which is presumed
//   to have been previously defined).
//
//   Implementation is as a derived class from the std::complex<T>
//   template from the STL library, which might be a horrendous idea,
//   which I will attempt to determine by trying it anyway and seeing
//   what problems might arise.  And even if all works, I'm not sure
//   what the long-term implications will be for portability,
//   platform-independence, and future-proof'ness.  Worst case
//   scenario, this class might at some point in the future need to be
//   implemented from scratch.  But I'm gonna cross my fingers and
//   leave that for another day that hopefully won't come...
//
class Complex : public std::complex<Real> {
public:
  Complex() : std::complex<Real>() {}
  Complex(Real x) : std::complex<Real>(x) {}
  Complex(Real x, Real y) : std::complex<Real>(x,y) {}
  Complex(const std::complex<Real> & other) : std::complex<Real>(other) {}
                    // (This last one allows the derived class to take
                    // assignment from objects of the base class.)

  Real abs() {return std::abs(*this);}
  Real arg() {return std::arg(*this);}
  Real norm() {return std::norm(*this);}  
                    // Norm here means magnitude-squared, which I'm
                    // not totally a fan of because it's ambiguous,
                    // but I'll stick with it because it is apparently
                    // the standard.
  
};

///
#endif //#infdef COMPLEX_H_
//
