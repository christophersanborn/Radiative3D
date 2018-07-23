// global.cpp
//
// All global objects should be instantiated here. (Defined/declared
// elsewhere, but instantiated here.)  This serves two purposes. One:
// provides an easy reference showing everything that's in the global
// namespace, and Two: solves the initialization-order problem.
// Global objects are initialized in the order in which they appear in
// a translation unit, but when several global objects exist across
// multiple translation units, their initialization order is
// undefined. This isn't really a problem for this project yet, since
// none of the globals depend on each other (though they may in the
// future), but nonetheless, I sleep better knowing the initialization
// order is known and controlable.  So we just put them all in here.
// But in every other aspect, the hpp/cpp structure remains the same.
//
#include "ecs.hpp"
#include "dataout.hpp"

//////
// GLOBAL OBJECTS:
//
EarthCoords  ECS;       // Earth Coordinate System
DataReporter dataout;   // Communication line to outside world

///
//
