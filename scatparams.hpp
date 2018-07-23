// scatparams.hpp
//
// This file develops a class that represents scattering parameters in
// the manner used by Sato and Fehler and provides computation methods
// to get quantities derived from those parameters (E.g., Sato's G
// functions, X functions, and P functions)
//
//
#ifndef SCATPARAMS_H_
#define SCATPARAMS_H_
//
#include "geom.hpp"       /* For S2:: and Geometry::Pi */
#include "elastic.hpp"    /* Elastic::HetSpec and 
                             Elastic::Velocity */

//////
// CLASS:   ::: ScatterParams :::
// ENCAPS:  Scattering Parameters
//
//   Encapsulates the set of parameters that specify a particular
//   scattering regime.  This includes such parameters as the scale
//   length and scattering strength, etc.
//
// USAGE PREREQUISITES:
//
//   Class member cm_omega needs to be set via SetFrequencyHertz()
//   before a ScatterParams object can be constructed from a HetSpec
//   and a Velocity alone. If the user attempts otherwise, an
//   exception of type std::logic_error will be thrown.
//
class ScatterParams {
private:

  // :::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Member Variables  (ScatterParams Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::
  //

  // Scattering Parameters:
  //
  //                            (Form borrowed from PSPhonon) 
  //
  //                            (In the future we may parameterize
  //                            differently, but for now just do what
  //                            Shearer and Earle did.)
  //

  // These values describe the heterogeneity spectrum, and are
  // provided to the constructor via a HetSpec object:

  Real  nu;         // density vs. velocity pert. scaling (see 4.48)
  Real  eps;        // RMS velocity perturbation (fractional)
  Real  a;          // correlation distance
  Real  kappa;      // von Karman parameter for PSDF

  // Values for the following are computed during construction from
  // knowlege of the background velocities (provided as an argument to
  // the constructor):

  Real  el;         // S-wave wavenumber (=om/beta0)
  Real  gam0;       // Pvel/Svel (often assumed = sqrt(3) )


private:

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Class-Private Members  (ScatterParams Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::

  static Real cm_omega;       // Class needs to be informed of phonon
  static bool cm_omega_known; // angular frequency omega before el can
                              // be computed. Ideally, this should be
                              // set prior to constructing any
                              // ScatterParams object.

public:

  // :::::::::::::::::::::::::::::::::::::::::::::
  // ::: Set-Methods for Class-Private Members :::
  // :::::::::::::::::::::::::::::::::::::::::::::

  static void SetFrequencyHertz(Real freq) {    // Use this to set the
    cm_omega = 2.0 * freq * Geometry::Pi;       // cm_omega member
    cm_omega_known = true;
  }


public:

  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Custom Data Types  (ScatterParams Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::


  // Define a few one-element enums to use as type-locked
  // symbols. These willserve as initializer flags to trigger
  // special-purpose constructors:
  //
  //         For example, use as: 
  //
  //             ScatterParams SP(ScatterParams::SATO_TEST);
  //
  //         to initialize a ScatterParams object with
  //         special-case constructor call.
  //

  enum SATO_TEST_e {SATO_TEST};
  enum SHEARER_TEST_e {SHEARER_TEST};
  enum SP_FORWARD_200_e {SP_FORWARD_200}; // Forward Scattering with
                                          // MFP_S =~ 200,
                                          // MFP_P =~ 600.


  // :::::::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (ScatterParams Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::

public: 

  // Constructor takes background velocities and a heterogeneity
  // spectrum, and computes the parameters needed for Sato and Fehler
  // formulation:

  ScatterParams(Elastic::Velocity vpvs,
                Elastic::HetSpec hs) :
    nu    ( hs.nu()    ),
    eps   ( hs.eps()   ),
    a     ( hs.a()     ),
    kappa ( hs.kappa() ),
    el    ( cm_omega/vpvs.Vs() ),   // Depends on frequency
    gam0  ( vpvs.Vp()/vpvs.Vs()) {
    if (!cm_omega_known)
      {AckSpit();}  // Barf if cm_omega not yet set
  }

  // This alternate constructor takes the heterogeneity spectrum and
  // then directly takes the el and gam0 paramters.  User won't
  // typically use this constructor for general earth modelling, but
  // could be used for more diagnostic purposes or for exploring the
  // parameter space.

  ScatterParams(Elastic::HetSpec hs,
                Real el, Real gam0) :
    nu    ( hs.nu()    ),
    eps   ( hs.eps()   ),
    a     ( hs.a()     ),
    kappa ( hs.kappa() ),
    el    ( el   ),
    gam0  ( gam0 )
  {}

  // These ones represent special test cases:

  ScatterParams(SATO_TEST_e dummy) :  // Params intended to replicate
    nu(0.8), eps(1.0), a(0.1),        // figure 4.8 in Sato & Fehler.
    kappa(0.5), el(1.0), gam0(1.7321)
  {}

  ScatterParams(SHEARER_TEST_e dummy) : // Params intended to be
    nu(0.8), eps(1.0), a(1.0),          // similar to those used by
    kappa(0.5), el(1.0),                // Shearer and Earle, 2004
    gam0(1.7321)                        // (TODO: Find out what those
  {}                                    // parameters would be - these
                                        // are just placeholder
                                        // numbers.)

  ScatterParams(SP_FORWARD_200_e dummy) : // Results in strong
    nu(0.8), eps(0.01), a(4.0),           // forward scattering and
    kappa(0.8), el(2.17411),              // MFP's (P,S) of ~600,
    gam0(1.7301)                          // ~200.
  {}


public:

  // ::::::::::::::::::::::::::::::::::::::::::::
  // ::: Member Access  (ScatterParams Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::

  Real GetNu() const {return nu;}
  Real GetEps() const {return eps;}
  Real GetA() const {return a;}
  Real GetKappa() const {return kappa;}
  Real GetL() const {return el;}
  Real GetGam0() const {return gam0;}


  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Comparison Methods  (ScatterParams Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::
  //
  //                            Provides metrics for comparing one set
  //                            of scatterparams to another, to
  //                            facilitate parsimony in allocating the
  //                            relatively memory intensive Scatterer
  //                            objects at model build time.
  //

  Real CompareRoughly(const ScatterParams & other) const;
  //    Provides a VERY rough comparison between two ScatterParams
  //    which should be used only for the first design iteration of
  //    optimizing Scatterer allocation.  The comparison is
  //    unsophisticated in that it only looks at the squared magnitude
  //    change in values, and sums them to an aggregate (possibly a
  //    weighted sum), but has no particular smarts about the affect
  //    each parameter change has on the actual resultant scattering
  //    behavior.  This, it might in practice be over/under sensitive
  //    to particualar parameters, which may have performance or
  //    realism implications.


  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Computational Result Methods  (ScatterParams Class) ::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //                                  These constant methods return
  //                                  functional mathematical results
  //                                  based on the parameters stored
  //                                  in this class.
  //
  //                                  They do not change the internal
  //                                  state of the object.
  //

  void GSATO (S2::S2Point toa,          // (A`la PSPhonon)
              Real & gpp, Real & gps, 
              Real & gsp, Real & gss, 
              Real & spol) const;

  void XSATO (S2::S2Point toa,          // (A`la PSPhonon)
              Real & xpp, Real & xps,
              Real & xsp, Real & xss_psi, 
              Real & xss_zeta) const;

  Real PSATO (Real m) const;            // Sato P function, kinda like
                                        // EXPSATO in PSPhonon, but we
                                        // use von Karman autocor


private:

  // :::::::::::::::::::::::::::::::::::::::::::::
  // ::: Error Handling  (ScatterParams Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::

  void AckSpit(); // Register a complaint if class used improperly


}; // class ScatterParams
////

///
#endif //#ifndef SCATPARAMS_H_
//
