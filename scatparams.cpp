// scatparams.cpp
//
#include <iostream>
#include <cmath>      /* for pow(), sin(), cos(), tgamma() */
#include "scatparams.hpp"

using namespace std;

// in this file:
// CLASS IMPLEMENTATIONS FOR:
//
//   o  Class ScatterParams
//


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  ScatterParams                                       ****
// ****                                                              ****
//

//
// STATIC MEMBERS:
//

Real ScatterParams::cm_omega = 1.0;         // Meaningless default - Actual
bool ScatterParams::cm_omega_known = false; // value is to be set in the
                                            // Model constructor via a call
                                            // to SetFrequencyHertz()

//////
// METHOD:  ScatterParams :: CompareRoughly()
//
//   Very rough comparrison between two ScatterParams objects,
//   returned as a sum of squared differences in parameter values.
//   This is very unsophisticated.
//
Real ScatterParams::CompareRoughly(const ScatterParams & other) const {
  Real dnu = (other.nu - nu);
  Real deps = (other.eps - eps);
  Real da = (other.a - a);
  Real dkappa = (other.kappa - kappa);
  Real del = (other.el - el);
  Real dgam0 = (other.gam0 - gam0);
  Real sum = dnu*dnu + deps*deps 
           + da*da + dkappa*dkappa 
           + del*del + dgam0*dgam0;
  return sum;
}


//////
// METHOD:  ScatterParams :: GSATO()
//
// PURPOSE: Computes (4.52) from Sato and Fehler using von Karman
// autocorellation function (PSATO).  The implementation here mimics
// the implementation in PSPhonon.
//
// Comment block from PSPhonon: (Note: Parameters el, nu, gam0, eps,
// and a are member variables and do not need to be passed, and psi
// and zeta are packed in toa.)
//
//! GSATO computes (4.52) from Sato and Fehler (exponential autocor)
//!   Inputs:  psi  =  spherical coor. angle (radians)
//!            zeta =  sph. coor angle from x3 (radians)
//!            el   =  S-wave wavenumber (=om/beta0)
//!            nu   =  density vs. velocity pert. scaling (see 4.48)
//!            gam0 =  Pvel/Svel (often assumed = sqrt(3) )
//!            eps  =  RMS velocity perturbation
//!            a    =  correlation distance
//!   Returns: gpp,gps,gsp,gss  =  from eqn. (4.52)
//!            spol =  S-to-S scattered S polarization (radians)
//!                 =  0 for pure psi direction
//!
void ScatterParams::GSATO(S2::S2Point toa,
                          Real & gpp, Real & gps, 
                          Real & gsp, Real & gss, 
                          Real & spol) const {
  const Real pi4  = 4. * Geometry::Pi;
  const Real el4  = pow(el,4);
  const Real gam2 = pow(gam0,2);
  const Real psi = toa.Theta();

  Real xpp, xps, xsp, xss_psi, xss_zeta;
  Real arg;

  XSATO(toa, xpp, xps, xsp, xss_psi, xss_zeta); // Compute X-values

  Real xpp2 = xpp * xpp;                      // Square them all
  Real xps2 = xps * xps;
  Real xsp2 = xsp * xsp;
  Real xss_psi2  = xss_psi * xss_psi;
  Real xss_zeta2 = xss_zeta * xss_zeta;

  // GPP:
  arg = (2.*el/gam0)*sin(psi/2.);               // For PSATO()
  gpp = (el4/pi4) * xpp2 * PSATO(arg);
  if (gpp < 1.e-30) gpp=0.; // PSPhonon sets this lower-limit
                            // value. Not sure why. But I retain
                            // it here.
  // GPS:
  arg = (el/gam0) * sqrt(1. + gam2 - 2.*gam0*cos(psi));
  gps = (1./gam0) * (el4/pi4) * xps2 * PSATO(arg);
  if (gps < 1.e-30) gps=0.;

  // GSP:
  gsp = gam0 * (el4/pi4) * xsp2 * PSATO(arg);
  if (gsp < 1.e-30) gsp=0.;

  // GSS:
  arg = 2.*el*sin(psi/2.);
  gss = (el4/pi4) * (xss_psi2 + xss_zeta2) * PSATO(arg);
  if (gss < 1.e-30) gss=0.;
  spol=atan2(xss_zeta,xss_psi); // Allows the whole range of
                                // polarization angles (-Pi to Pi).
  return;

}


//////
// METHOD:  ScatterParams :: XSATO()
//
// PURPOSE: Computes (4.50) from Sato and Fehler
//
// COMMENT BLOCK FROM PSPhonon: (Note: Parameters nu and gam0 are
// implemented here as member variables and do not need to be passed
// in; and note also that psi and zeta are both packed into
// "take-off-angle" toa as toa.phi() and toa.theta(), respectively.)
//
//! XSATO computes (4.50) from Sato and Fehler
//!   Inputs:  psi  =  spherical coor. angle (radians)
//!            zeta =  sph. coor angle from x3 (radians)
//!            nu   =  density vs. velocity pert. scaling (see 4.48)
//!            gam0 =  Pvel/Svel (often assumed = sqrt(3) )
//!   Returns: xpp, xps, xsp, xss_psi, xss_zeta  =  from eqn. (4.50)
//!
void ScatterParams::XSATO(S2::S2Point toa,
                          Real & xpp, Real & xps, 
                          Real & xsp, Real & xss_psi, 
                          Real & xss_zeta) const {

  const Real gam2  = gam0 * gam0;
  const Real cpsi  = cos(toa.Theta());     // psi <--> toa.Theta()
  const Real c2psi = cos(2.*toa.Theta());
  const Real spsi  = sin(toa.Theta());
  const Real czeta = cos(toa.Phi());       // zeta <--> toa.Phi()
  const Real szeta = sin(toa.Phi());
  const Real spsi2 = spsi*spsi;

  xpp = (1./gam2) * (nu * (-1. + cpsi + (2./gam2) * spsi2) 
                     - 2. + (4./gam2) * spsi2);

  xps = -spsi*(nu*(1.-(2./gam0)*cpsi)-(4./gam0)*cpsi);

  xsp = (1./gam2) * spsi*czeta*(nu*(1.-(2./gam0)*cpsi) 
                                - (4./gam0)*cpsi);

  xss_psi = czeta*(nu*(cpsi-c2psi)-2.*c2psi);

  xss_zeta = szeta*(nu*(cpsi-1.)+2.*cpsi);

  return;

}


//////
// METHOD:  ScatterParams :: PSATO()
//
// PURPOSE: Compute Sato and Fehler's P(m) function using von Karman
// autocorrelation.  (Note that when k=0.5, von Karman is equivalent
// to exponential, which is what is used in PSPhonon.)
//
// COMMENT BLOCK from PSPhonon: (Note that eps and a are member
// variables of the class, and thus do not need to be passed in.)
//
//! function EXPSATO computes (2.10) from Sato and Fehler
//!    Inputs:  eps  =  RMS velocity perturbation
//!             a    =  correlation distance
//!             m    =  wavenumber
//!    Returns: P(m) =  PSDF (Power Spectral Density Function)
//!
Real ScatterParams::PSATO(Real m) const {
  const Real pi32  = pow(Geometry::Pi,1.5); // pi^(3/2)

  const Real numer = (8.*pi32*eps*eps*a*a*a)  
                   * tgamma(kappa+1.5)
                   / tgamma(kappa);

  const Real denom = pow((1.+a*a*m*m),(kappa+1.5));

  return numer/denom;
}


//////
// METHOD:   ScatterParams :: AckSpit()
//
//   Raise a complaint if class used improperly.
//
void ScatterParams::AckSpit() {

  throw(std::logic_error(
    "Can't construct ScatterParams before frequency is known."));

}
