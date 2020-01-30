// raypath.hpp
//
// Defines classes that define specific geometric "ray paths" that a phonon
// may follow. These are essentially simple geometric primatives, such as the
// center and radius of a circular arc.  These paths are dependent on the
// velocity profile in a medium cell, thus are generally computed by methods
// in the MediumCell hierarchy, but are passed to methods in the CellFace
// hierarchy to, e.g., find points on intersection with the boundaries of the
// cells.
//
// Note that for the moment it's only circular arcs computed in spherical
// velocity gradients that are using these classes.  But this seems to be a
// clean organization, so it may make some sence to go back and retrofit the
// linear ray path code and the circular arc Tetra code to make use of "ray
// path" classes here to describe paths.
//
#ifndef RAYPATH_H_
#define RAYPATH_H_
//
#include <string>
#include "geom.hpp"

//////
// CLASS:  cache_RD2_precompute
//
// When we compute the Ray Arc geometry in Radial-Quadratic (RD2)
// velocity profiles, we cache in this structure certain precomputed
// values that are going to be used in multiple subsequent
// calculations.
//
struct cache_RD2_precompute {
  Real S;           // Model center to Arc center Separation distance.
  Real S2;          // S^2
  Real TwoSQ;       // 2*S*Q; (Q being Arc Radius)
  Real CosZeta;     // Zeta is angle between line S and a tangent to the v==0
                    // isosurface through the Arc Center. Cosine of that.
  Real SinZeta;     // Sine of Zeta
  Real CotZetaBy2;  // Cotangent of Zeta/2
  Real timeCoef;    // Coefficient used in computing travel times: equal to
                    // (a*S*SinZeta)^-1 where 'a' is the quadratic coefficient
                    // of the velocity profile (i.e. from v=ar^2+c).
  cache_RD2_precompute() {}
  // Specify via just four parameters:
  cache_RD2_precompute(Real _S2, Real _Q, Real _zeroRad2, Real _a) :
    S          (  sqrt(_S2)  ),
    S2         (  _S2        ),
    TwoSQ      (  2*S*_Q     ),
    CosZeta    (  (S2 + _Q*_Q - _zeroRad2)/TwoSQ  ),
    SinZeta    (  sqrt(1 - CosZeta*CosZeta)       ),
    CotZetaBy2 (  (1+CosZeta)/SinZeta             ),
    timeCoef   (  1/(_a*S*SinZeta)                )
    {}
};

//////
// CLASS:  RayArcAttributes
//
// Defines a CIRCULAR ray arc with a known Center and Radius, and a coordinate
// system anchored on the arc center, wherein u3 points along the velocity
// gradient (towards what could be called a "ray bottoming"), and u1 points in
// the direction of the ray tangent when the ray is at the bottoming
// point. The coordinate system can be used to define positions along the ray
// in terms of a either theta-coordinate or an arc-length coordinate, where
// positive values are beyond the bottoming (and uptrending) and negative
// values are before the bottoming (and downtrending).
//
// This struct is suitable for ray arcs in Radial-Quadratic (RD2) velocities
// or cartesian Linear-gradient (LD1), and possibly others, as these produce
// circular ray arcs.  Additional profile-dependent precomputes are stored in
// a union cache.
//
struct RayArcAttributes {

  Real Radius;        // Radius of the ray arc (from arc center).
  Real Rad2;          // Radius-squared

  R3::XYZ Center;     // Centerpoint of the ray arc.

                      // REFERENCE BASIS anchored on ARC CENTER:
  R3::XYZ u3;         // Points towards ray "bottom" (i.e. in direction of
                      // velocity gradient).
  R3::XYZ u2;         // Out-of-plane.
  R3::XYZ u1;         // In-plane perpendicular. Aligns with ray tangent
                      // at arc bottom.

  union cache {
    cache_RD2_precompute RD2; // Used for radial-quadratic arcs.
    cache() {}
  } c;

  RayArcAttributes() {}

  // NOTES: If ray is pure vertical, then should hold Bottom=0,
  //        Radius=+inf, Center=(0,0,0), u3=u2=(0,0,0),
  //        u1 = direction of ray tangent.
  // TODO: Finalize these "special case" decisions.

  Real AngleOffsetFromBottom(const R3::XYZ & loc) const;
        // Given a location 'loc' assumed to be in the arc plane, return the
        // angle made between center-bottom and center-loc. Positive angles
        // mean loc is ahead of the bottom, negative behind, determined by
        // assuming u1 contains the ray tangent at bottom.

  R3::XYZ PositionFromAngle(Real angle) const;
        // Given angle relative to bottom, return location on arc path.

  R3::XYZ DirectionFromAngle(Real angle) const;
        // Given angle relative to bottom, return direction of the ray
        // tangent at that point.

  std::string str() const;

};

///
#endif //#ifndef RAYPATH_H_
//
