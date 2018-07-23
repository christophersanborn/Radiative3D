// geom.hpp
///@file
///
/// This is a top-level include for the Geometry library of classes
/// developed to support Radiative3D.  Including this file gets you all
/// of the following:
///
///     *  typedefs.hpp
///     *  geom_base.hpp, defining the Geometry namespace
///     *  geom_r3.hpp, defining the R3 namespace
///     *  geom_s2.hpp, defining the S2 namespace
///
///   NOTE: geom_r4.hpp is *not* included by this header, as so few
///     modules actually need it. If the R4 namespace is needed, the
///     module should additionally include geom_r4.hpp after geom.hpp.
///
#ifndef GEOM_H_
#define GEOM_H_
//
#include "geom_r3.hpp"
#include "geom_s2.hpp"
//
//  No need to directly include "geom_base.hpp" since it is
//  included indirectly via geom_r3.hpp and geom_s2.hpp.
//

//////
// *** TYPEDEFS:
//

namespace S2 {
typedef ThetaPhi      S2Point;  /* Deprecated */
typedef ThetaPhi_Set  S2Set;    /* Deprecated */
                      // Deprecation Note: I don't like the names of
                      // these typedefs, or even necessarily that they
                      // exist.  I have removed them from the geom lib
                      // files themselves, and retain these only in
                      // this header to preserve application code
                      // already written, but that code SHOULD be
                      // re-written.
};


///
#endif //#ifndef GEOM_H_
//
