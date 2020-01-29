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
// CLASS:  RayArcAttributes
//
// Defines a CIRCULAR ray arc with a known Center and Radius, and a coordinate
// system anchored on the arc center, wherein u3 points along the velocity
// gradient (towards what could be called a "ray bottoming"), and u1 points in
// the direction of the ray tangent when the ray is at the bottoming
// point. The coordinate systme can be used to define positions along the ray
// in terms of a either theta-coordinate or an arc-length coordinate, where
// positive values are beyond the bottoming (and uptrending) and negative
// values are before the bottoming (and downtrending).
//
struct RayArcAttributes {

  Real Bottom;        // Bottoming radial coord of ray arc. (Dist from model
                      // center.)  DEPRECATE - (Not really relevant outside of
                      // spherical velocity gradients)

  Real Radius;        // Radius of the ray arc (from arc center).

  R3::XYZ Center;     // Centerpoint of the ray arc.

                      // REFERENCE BASIS anchored on ARC CENTER:
  R3::XYZ u3;         // Points towards ray "bottom" (i.e. in direction of
                      // velocity gradient).
  R3::XYZ u2;         // Out-of-plane.
  R3::XYZ u1;         // In-plane perpendicular. Aligns with ray tangent
                      // at arc bottom.

  // NOTES: If ray is pure vertical, then should hold Bottom=0,
  //        Radius=+inf, Center=(0,0,0), u3=u2=(0,0,0),
  //        u1 = direction of ray tangent.
  // TODO: Finalize these "special case" decisions.

  Real AngleOffsetFromBottom(const R3::XYZ & loc) const;
        // Given a location 'loc' assumed to be in the arc plane, return the
        // angle made between center-bottom and center-loc. Positive angles
        // mean loc is ahead of the bottom, negative behind, determined by
        // assuming u1 contains the ray tangent at bottom.

  std::string str() const;

};

///
#endif //#ifndef RAYPATH_H_
//
