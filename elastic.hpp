// elastic.hpp
//
// Defines classes to represent physical dimensional quantities in
// elastic media, and in some cases codifies or automates the
// relationship between these quantities (e.g. in cases where
// quantities are not independent).
//
// The ideal manifestation of this, from a theory standpoint, would be
// a class called, eg, Moduli, encoding the shear and bulk moduli as
// complex tensorial quantities, and material density as a scalar
// quantity.  From this information, all needed elastic attributes (P
// and S velocities (perhaps as functions of direction), density, and
// P and S attenuation factors (also perhaps directionally dependent))
// could be derived.  Ie, by some means we would initialize a Moduli
// object, and then we would ask it for needed info.
//
// A meta-class called Elastic would enclose the Moduli and also a set
// of heterogeneity specifications (ultimately to be turned into
// scattering parameters).  This Elastic object, or perhaps a
// flattening ro reduction of it, would then be suitable for use as a
// grid node.
//
// We will build out these classes gradually, and will not initially
// tackle tensorial representation.  Ie, from the outset, we are
// assuming isotropic media.  But we're gonna construct things in as
// generic a way as possible, so that interface changes will be
// minimized when we transition to non-isotropic modelling.
//
// Note, the purpose of this module is to simplify and clarify
// representation of complex multidimensional quantities for use in
// interpreting or translating physical data.  The purpose of this
// module is NOT to optimize either storage or speed.  These classes
// are intended to facilitate data input/output from the program,
// e.g., in the grid-building process.  But for code that needs to run
// FAST, (i.e., inner simulation loop), it may be wise to NOT use
// these classes.  Encoding just the quantities needed for computation
// as Reals will most likely be preferable.
//
// Classes in this module fall into two categories: "complete", and
// "contingent". A "complete" class is one where the class contains
// all information necessary to respond to all interrogations with no
// arguments to the accessor functions. A "contingent" class is one
// where arguments must be provided to the accessor functions.  An
// example of a complete class is the Veloc class, storing P and S
// velocities.  The accessor functions simply return the values.  An
// example of a contingent class is the Q class. Here, Qp, Qs, and Qk
// are all related to each other via a formula in which the velocities
// are a parameter.  Thus, the Qx() accessor functions take a Veloc
// argument so that the requested Q value can be solved for if it
// isn't cached in the Q object.  The reason to have a meta-class like
// Elastic is because Elastic is a "complete" class, which contains
// both Velocities and a contingent Q object.  The Qx() accessors to
// the Elastic class do not require arguments, since the member Veloc
// objet becomes the argument to the contingent Q accessor functions.
//
// HIERARCHIES:
//
// Base classes in this module are aschematic, which means they
// provide no public constructors allowing for their values to be
// initialized by arguments.  Rather, they must be initialized through
// special constructor-only derived classes called schemas. Use of
// schema classes lexical clarity when initializing objects, as the
// schema class name makes clear the expected arguments and their
// ordering.  (Ie, if I want to specify Q by giving Q_mu and Q_kappa
// (as opposed to Q_p and Q_s, e.g.), I would use the QmQk schema, as
// this makes the meaning of the arguments clear.)
//
// CLASSES:
//
// * Complete:
//
//   class Velocity
//   class Density
//   class Moduli (*Planned*)
//
//
// * Contingent:
//
//   class Q          // Attenuation Q - depends on Veloc
//   class HetSpec    // Heterogeneity spectrum
//
//
// * Meta:
//
//   class Elastic    // Contains Veloc, Dens, Q
//   class HElastic   // Contains Veloc, Dens, Q, and HetSpec
//
//
// * Schemas:  (Constructor-only derived classes)
//
//   class VpVs : Velocity  // Set velocities as a (Vs, Vp) tupple
//   class Qinf : Q         // Set Qs and Qp both infinite
//   class QpQs : Q         // Set Q's as a (Qp, Qs) tupple
//   class QmQk : Q         // Set Q's as a (Qmu, Qkappa) tupple
//   class QpQk : Q         // Set Q's as a (Qp, Qkappa) tupple
//
//
#ifndef ELASTIC_H_
#define ELASTIC_H_
//
#include "typedefs.hpp"

namespace Elastic {

class Velocity {
private:

  Real mVp;             // P and S velocity.  In future, we may encode
  Real mVs;             // tensor elements for non-isotropic velocities.

protected:

  Velocity(Real Vp, Real Vs) :  // Protected; for derived-class use only
    mVp (Vp),
    mVs (Vs)
  {}

public:

  Velocity() :          // Default constructor => 0,0
    mVp (0.),           // Use derived classes to initialize values
    mVs (0.)            //
  {}

  Real Vp() const {return mVp;}
  Real Vs() const {return mVs;}

};


class VpVs : public Velocity {
public:

  VpVs(Real vp, Real vs) :
    Velocity (vp,vs)
  {}

};


class Density {
private:
  Real mDens;
public:
  Density(Real dens) :
    mDens(dens)
  {}
  Real Value() const {return mDens;}
};


class Q {
protected:
  enum q_unknown_e {
    QUNK_QP,            // Given any two, we solve for the third.
    QUNK_QS,            // This enum identifies the unknown quantity.
    QUNK_QK             //
  };

private:
  q_unknown_e mUnknown; //
  Real mQp;             // P-wave Q
  Real mQs;             // S-wave Q  (Also called Q_mu or shear modulus Q)
  Real mQk;             // Bulk modulus (Kappa) Q
  
protected:
  Q(q_unknown_e U, Real Qp, Real Qs, Real Qk) : // Protected; for derived-
    mUnknown(U),                                // class use only.
    mQp(Qp),
    mQs(Qs),
    mQk(Qk)
  {}

public:
  Q() :                 // Default constructor, all Q's infinite.
    mUnknown(QUNK_QK),  // Note: no other base-class constructor provided;
    mQp (1./0.),        // use derived classes to initialize Q objects.
    mQs (1./0.),        //
    mQk (1./0.)         //
  {}

  Real Qp(Velocity vpvs) const;     // Return Q values, solving for unknown if
  Real Qs(Velocity vpvs) const;     // necessary (velocities needed if value
  Real Qk(Velocity vpvs) const;     // unknown).

private:
  static Real HelperL(Velocity vpvs);

};


class Qinf : public Q {         // Explicit way to choose infinte Q
public:                         // values (i.e., no attenuation.)
  Qinf() :
    Q (QUNK_QK, 1./0., 1./0., 0)
  {}
};

class QpQs : public Q {         // Specify Q_p and Q_s
public:
  QpQs(Real Qp, Real Qs) :
    Q (QUNK_QK, Qp, Qs, 0)
  {}
};


class QmQk : public Q {         // Specify Q_mu and Q_kappa
public:                         // (Q_mu is the same as Qs. (mu is
  QmQk(Real Qs, Real Qk=1./0.)  // shear modulus; kappa is bulk
    : Q(QUNK_QP, 0, Qs, Qk)     // modulus.)) AK135 specifies Q's this
  {}                            // way.
};


class QkQm : public Q {         // Similar to QmQk but reverses input
public:                         // order.  (As in AK135f)
  QkQm(Real Qk, Real Qm)        //
    : Q(QUNK_QP, 0, Qm, Qk)
  {}
};


class QpQk : public Q {         // Specify Q_p and Q_kappa
public:
  QpQk(Real Qp, Real Qk=1./0.)
    : Q (QUNK_QS, Qp, 0, Qk)
  {}
};


class HetSpec {
private:

  Real mNu;     // 
  Real mEps;    // Heterogeneity amplitude (fractional)
  Real mA;      // Spatial scale-length (corner)
  Real mKappa;  // Von Karman K

protected:

  HetSpec(Real nu, Real eps, Real a, Real k)  // protected: use schemas to
    : mNu(nu), mEps(eps), mA(a), mKappa(k)    // set values
  {}

public: 

  HetSpec()                         // A default option - (eps zero)
    : mNu(0.8), mEps(0.00),         //
      mA(1.0), mKappa(0.5)
  {}

  Real nu() const {return mNu;}
  Real eps() const {return mEps;}
  Real a() const {return mA;}
  Real kappa() const {return mKappa;}

};


class HSneak : public HetSpec {     // Schema explicitly denoting arg order
public:                             // as nu, eps, a, kappa (n,e,a,k)
  HSneak(Real nu, Real eps, Real a, Real k) :
    HetSpec(nu,eps,a,k)
  {}
};


class SElastic {     // Simple elastic, no heterogeneity
protected:

  Velocity mV;
  Density  mRho;
  Q        mQ;

public:;

  SElastic(Velocity vpvs, Density dens, Q q) :
    mV(vpvs), mRho(dens), mQ(q) // Note: Density can also be passed
  {}                            // as a Real. It will be promoted.

  SElastic(Density dens, Velocity vpvs, Q q)  // AK135 specifies Density be-
    : mV(vpvs), mRho(dens), mQ(q)             // fore velocity.  Typelocking
  {}                                          // disambiguates.

  Velocity getV() const {return mV;}
  Density getDens() const {return mRho;}
  Q getQ() const {return mQ;}

  Real Vp() const {return mV.Vp();}
  Real Vs() const {return mV.Vs();}
  Real Rho() const {return mRho.Value();}
  Real Qp() const {return mQ.Qp(mV);}
  Real Qs() const {return mQ.Qs(mV);}
  Real Qk() const {return mQ.Qk(mV);}

};

class HElastic : public SElastic {      // Adds in heterogeneity spectrum
protected:

  HetSpec mH;

public:;

  HElastic(Velocity vpvs, Density dens, Q q, HetSpec h)
    : SElastic(vpvs, dens, q), 
      mH(h)
  {}

  HElastic(Density dens, Velocity vpvs, Q q, HetSpec h)
    : SElastic(vpvs, dens, q), 
      mH(h)
  {}

  HetSpec getHS() const {return mH;}

};

}; // namespace Elastic
////

///
#endif //#ifndef ELASTIC_H_
//
