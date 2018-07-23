// typedefs.hpp
//
// This file defines elemental types to be used throughout the
// Radiative3D project.  The goal to is provide abstract numeric types
// that remove decisions of machine implementation from the core
// coding process.  For example, when coding the Radiative3D
// application-layer code, I do not want to have to be deciding
// between single- and double-precision floating point types, for
// example.  That choice is made here, in this file, and applies
// throughout. This also makes it easy to make an application-wide
// change later on, and assess the performance and precision
// implications easily, and then make the best choice.  Several
// abstract types are defined here, covering both real-number and
// integer use cases.
//
// This file has recently been expanded to additionally include string
// types and exception types.
//
#ifndef TYPEDEFS_H_
#define TYPEDEFS_H_
//
#include <string>
#include <sstream>
#include <stdexcept>
#include "config/opt-fptype.hpp"    // Choose float-type as
                                    // build-time option
#ifndef TYPEDEFS_H_FP_TYPE          //
#define TYPEDEFS_H_FP_TYPE double   // Use "double" as default
#endif                              //


//
// ***************************
// ***(  Indicial Types:  )***
// ***************************
//

typedef unsigned int             Index;
typedef unsigned long int     BigIndex;
typedef unsigned short int  SmallIndex;
//
//          Index Type: this type is intended to
//          represent sequentially enumerated indices
//          beginning at 0.
//
//          Use "Index" for situations in which the
//          machine-native integer type is considered
//          "big enough".
//
//          Use "BigIndex" For situations where we may
//          be indexing a very large set, and need the
//          biggest machine-supported integer type.
//
//          Use "SmallIndex" for situations where we
//          KNOW the indexed set will be not more than a
//          handful of elements.
//

typedef signed int            RelIndex;
typedef signed long int    BigRelIndex;
typedef signed short int SmallRelIndex;
//
//          Relative Index Type: this type is intended to
//          represent differences between sequentially
//          enumerated indices.  It is therefore a signed
//          type, to facilitate negative values.
//

//
// **************************
// ***(  Numeric Types:  )***
// **************************
//

typedef unsigned int             Count;  
typedef unsigned long int     BigCount;
typedef unsigned short int  SmallCount;
//
//          Count Type: for counting things.  Same
//          underlying type as Index, but using a
//          different name to preserve semantic
//          distinction. "Count" is more like an
//          integral "extent," whereas "Index" is more
//          about a set of identifiers.
//
//          Big and Small variants provided on the same
//          logic as for Index.
//

typedef  TYPEDEFS_H_FP_TYPE     Real;
//
//          Real Type: for representing "real"
//          (non-integral) values.  We use this typedef
//          to choose the floating-point type to use for
//          these values.
//
//
//          A typedef is also provided for
//          complex-values. See complex.hpp for this
//          implementation.
//

//
// ***********************
// ***(  Char Types:  )***
// ***********************
//

typedef  std::string                Text;
typedef  std::stringstream    TextStream;
//
//          Notational shorthand for string types.  Inclusion
//          of this header also guarantess that the headers
//          neaded for these type are included.
//

//
// ******************************
// ***(  Operational Types:  )***
// ******************************
//

typedef  std::runtime_error      Runtime;
typedef  std::logic_error        Invalid;
//
//          Here I define exactly TWO exception data types, and an
//          attempt will be made to use ONLY these two type throughout
//          the code. (Otherwise, choosing which 'standard' exception
//          type is appropriate to each situation is a heady and
//          brain-clenching mess.)  The basic schema will be to use
//          'Runtime' for situations in which the program throws
//          because the user has provided malformed input (e.g. bad
//          command line args, or input data that cannot be parsed or
//          handled, etc.), and 'Invalid' will be used when a function
//          throws because the args provided to the function are
//          malformed according to the needs of the function.
//
//          Said a different way: 'Runtime' implies that the user did
//          something wrong, whereas 'Invalid' implies that the
//          software author did something wrong... i.e. it suggests
//          that a bug has been discovered.
//
//          In the future, these may be implemented as derived types.
//          But typedefs suffice for now.
//

///
#endif //#infdef TYPEDEFS_H_
//
