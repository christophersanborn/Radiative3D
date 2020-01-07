// grid.cpp
//
#include <iostream>
#include <iomanip>    /* std::setw(), etc. */
#include <sstream>    /* std::ostringstream */
#include <cstdlib>    /* exit() */
#include <vector>
#include <stdexcept>
#include "grid.hpp"


// In this file:
// CLASS IMPLEMENTATIONS FOR:
//
//   o  Class GridData
//   o  Class GridNode
//   o  Class Grid
//
// Search on "&&&&" to jump between class implementations in this
// file.
//

//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  GridNode                                            ****
// ****                                                              ****
//

//////
// METHOD:  GridNode :: SetLocation()
//
//   Sets the value of mLocation.  Effective available only once per
//   grid node - subsequent calls are silently ignored.  (Detected
//   heuristically by comparing with default 0,0,0 value.)  This is
//   because the grid-definition file format permits a particular node
//   to be defined up to two times, to encode sharp jumps in seismic
//   parameters across layers. It is meaningful to specify seismic
//   parameters more than once, but not location.  This test gives
//   precedence to the FIRST defintion of location, and allows the
//   user to either repeat the location for the second-definition, or
//   to use 0,0,0 shorthand. (Alternate shorthands are technically
//   permissible, but could circumvent the test if a particular node
//   really is located at 0,0,0)
//
//   The tuple (x,y,z) can be in any coordinate system supported by the
//   EarthCoords ECS and will be on-the-fly converted to model space
//   coordinates upon readout with the Loc() method.
//
void GridNode::SetLocation(Real x, Real y, Real z) {

  if (mLocation.IsNull()) {     // Prevents subsequent redefinition

    mLocation.SetTriple(x,y,z);

  }

}

void GridNode::AdjustLocation(Real x, Real y, Real z) {
  mLocation.SetTriple(x + mLocation.x1(),
                      y + mLocation.x2(),
                      z + mLocation.x3());
}


//////
// METHOD:  GridNode :: SetAttributes()
//
//   Populates up to two GridData objects to represent elastic attributes
//   at nodal location (in case of single definition) or just above and
//   just below node location (in case of double definition).
//
void GridNode::SetAttributes(GridData gdata) {

  if (mNextData == GN_NLAY) {   // If both slots already full, then throw:
    throw(std::runtime_error(   //
      "GridNode: SetAttributes: Too many definitions for GridNode.\n"));
  } // Else slots are available, so assign:
  mData[mNextData] = gdata;                       // Assign attributes.
  mNextData = (mNextData == GN_ABOVE) ? GN_BELOW  // Increment mNextData.
                                      : GN_NLAY;  //

}


//////
// METHOD:  GridNode :: Loc()
//
//   Returns model space XYZ location of the node. (Implies conversion
//   from user ECS coordinate system.)
//
R3::XYZ GridNode::Loc() const {
  return ECS.Convert(mLocation);
}


//////
// METHOD:  GridNode :: Data()
//
//   Returns GridData object to represent elastic attributes at node
//   location.  'lay' argument applies when node is discontinuous
//   (i.e. when two distinct attributes have been specified).  In this
//   case, 'lay' specifies whether we want attributes just "above" node or
//   just "below" node.
//
GridData GridNode::Data(layers_e lay) const {

  if (mNextData == GN_ABOVE) {  // If slots are empty, then throw:
    throw(std::runtime_error(   //
      "GridNode: Data: No attributes set for GridNode.\n"));
  }

  GridData selected;            // Choose data object to return:
 
  if (mNextData == GN_BELOW) {  // If only one definition,
    selected = mData[GN_ABOVE]; // then select it, regardless of 'lay';
  } else {                      //
    selected = mData[lay];      // else select the requested object.
  }                             //

 return ECS.Convert(mLocation, selected);   // Convert() will apply Earth
                                            // Flattening, if appropriate.

}


//////
// METHOD:  GridNode :: OutputAsAscii()
//
//   Outputs to the grid node to an an ascii output stream.  Locations
//   and properties are presented in the user's choice of Output
//   Coordinate System (OCS). This involves a conversion from the
//   Earth Coordinate System (ECS) in which the node locations and
//   data are stored.  (See ecs.hpp for more.)
//
void GridNode::OutputAsAscii(std::ostream & out, 
                             const std::string & prefix) const {

  int num_lines = 2;
  if (mNextData == GN_BELOW) {num_lines = 1;}
  if (mNextData == GN_ABOVE) {num_lines = 0;}

  std::ios_base::fmtflags flags =   // Save ostream state
    out.flags();

  out.precision(5);
  out << std::fixed << std::right;

  R3::XYZ iloc = ECS.Convert(mLocation);            // Convert from ECS to ICS
  EarthCoords::Generic oloc = ECS.OutConvert(iloc); // ...and from ICS to OCS.

  if (num_lines == 0) {             // Then data is missing...
    out << prefix 
        << std::setw(11) << oloc.x1() << " "
        << std::setw(11) << oloc.x2() << " "
        << std::setw(11) << oloc.x3() << "  "
        << "        ***       ***       ***         ***       *** "
        << "        ***       ***       ***       ***\n";
  }
  else {
    for (int iside = 0; iside < num_lines; iside++) {
      GridData idat = ECS.Convert(mLocation,mData[iside]);  // ECS to ICS
      GridData odat = ECS.OutConvert(iloc,idat);            // ICS to OCS
      out << prefix
          << std::setw(11) << oloc.x1() << " "
          << std::setw(11) << oloc.x2() << " "
          << std::setw(11) << oloc.x3() << "    "
          << std::setw(9) << odat.Vp() << " "
          << std::setw(9) << odat.Vs() << " "
          << std::setw(9) << odat.Rho() << "   ";
      out.precision(1);
      out << std::setw(9) << odat.Qp() << " "
          << std::setw(9) << odat.Qs() << "   ";
      out.precision(5);
      out << std::setw(9) << odat.getHS().nu() << " "
          << std::setw(9) << odat.getHS().eps() << " "
          << std::setw(9) << odat.getHS().a() << " "
          << std::setw(9) << odat.getHS().kappa() << "\n" ;
    }
  }
                    
  out.flags(flags);                 // Restore ostream state
}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  Grid                                                ****
// ****                                                              ****
//
// IMPLEMENTATION NOTES:
//
// o  GridNode Storage:
//
//    The GridNode's are stored in a linear array (vector<GridNode>),
//    but are accessed outside of the class via three indexes
//    representing a 3D cartesian grid array. We pack the nodes in
//    contiguous memory such that the x-index varies fastest, and the
//    z-index slowest.
//

//////
// CONSTRUCTOR:
//
Grid::Grid() :
  mModelTarget (MOD_AUTO),
  mN (0), mNi (0), mNj (0), mNk (0),
  mIndexBase ( 0 ) 
{
} 

//////
// METHOD:   Grid :: SetSize()
//
void Grid::SetSize(Count ni, Count nj, Count nk) {

  if (mN != 0) {
    std::cerr << "Error: Attempt to resize an already-sized grid.\n";
    exit(1);  // Abort
  } 
  else {
    mNi = ni;
    mNj = nj;
    mNk = nk;
    mN  = ni * nj * nk;
    mGrid.clear();
    mGrid.resize(mN);
 }

}
  

//////
// METHOD:  Grid :: SetIndexBase()
//
void Grid::SetIndexBase(Index base) {
  mIndexBase = base;
}


//////
// METHOD:  Grid :: SetMapping()
//
//   Sets up the correct mapping scheme in the ECS based on the user's
//   choice of coordinate system and curvilinear mapping.
//
//   The Grid class uses two enums to select mapping, separating
//   coordinate choice and curvature into two questions.  The EarthCoords
//   class (provides the global ecs object), on the other hand, just has
//   one enum to select "mapping system".  In some sense, this function is
//   just a translation between the two schemes, providing an interface to
//   the user that saves them from needing to know any particulars about
//   the EarthCoords implementation.
//
//   Note that some combinations of coordinate choice and curvature are
//   not sensible and/or are not available as mappings in the ECS. These
//   combinations will result in a mapping choice of
//   EarthCoords::MAP_UNSUPPORTED which will trigger exceptions if any
//   coordinate conversions are requested from the ECS.
//
//
void Grid::SetMapping(gs_coords_e coords, curvature_e curve) {

  if (coords==GC_ENU && curve==GC_ORTHO) {
    ECS.SetMapping(EarthCoords::ENU_ORTHO);
    ECS.SetEarthFlattening(false);
  } 
  else if (coords==GC_ENU && curve==GC_FLATTENED) {
    ECS.SetMapping(EarthCoords::ENU_ORTHO);
    ECS.SetEarthFlattening(true);
  }
  else if (coords==GC_RAE && curve==GC_ORTHO) {
    ECS.SetMapping(EarthCoords::RAE_ORTHO);
    ECS.SetEarthFlattening(false);
  }
  else if (coords==GC_RAE && curve==GC_FLATTENED) {
    ECS.SetMapping(EarthCoords::RAE_ORTHO);
    ECS.SetEarthFlattening(true);
  }
  else if (coords==GC_RAE && curve==GC_CURVED) {
    ECS.SetMapping(EarthCoords::RAE_CURVED);
    ECS.SetEarthFlattening(false);
  }
  else {
    ECS.SetMapping(EarthCoords::MAP_UNSUPPORTED);
  }

}


//////
// METHOD:  Grid :: WNode()
//
//   Returns a writable reference to the node indicated by the
//   indices. Index-base adjustment is applied to the indices.
//
//   Bounds checking is performed on the indices.
//
GridNode & Grid::WNode(Index i, Index j, Index k) {
  i -= mIndexBase;
  j -= mIndexBase;
  k -= mIndexBase;
  if ((i >= mNi) || (j >= mNj) || (k >= mNk)) {
    std::cerr << "ERROR: Grid index out of bounds.\n";
    exit(1);
  }
  Index index = k*(mNj*mNi) + j*(mNi) + i;
  return mGrid[index];
}


//////
// METHOD:  Grid :: Node()
//
//   Returns a read-only reference to the node indicated by the
//   indices. Index-base adjustment is NOT applied to the
//   indices. (Assumes base 0 indexing.)
//
//   Bounds checking is NOT performed on the indices.
//
const GridNode & Grid::Node(Index i, Index j, Index k) const {
  Index index = k*(mNj*mNi) + j*(mNi) + i;
  return mGrid[index];
}


//////
// METHOD:  Grid :: RelNode()
//
//   Returns a read-only reference to the node indicated by the
//   indices, as adjusted by the relative indices. Index-base
//   adjustment is NOT applied to the indices. (Assumes base 0
//   indexing.)
//
//   Bounds checking is NOT performed on the indices.
//
const GridNode & Grid::RelNode(Index i, Index j, Index k,
                               RelIndex ri, RelIndex rj, RelIndex rk) const {
  Index index = (k+rk)*(mNj*mNi) 
              + (j+rj)*(mNi) 
              + (i+ri);
  return mGrid[index];
}


//////
// METHOD:  Grid :: GetModelType()
//
//   Returns a code that tells the Model() class what build-out method
//   to use to build the model from the grid.  The user can set this
//   code when building the grid, or this method can make a guess
//   based on the layout of the grid as defined.
//
Grid::model_target_e Grid::GetModelType() const {
  if (mModelTarget != MOD_AUTO) {
    return mModelTarget;
  } else {
    if ((mNi==3) && (mNj==1) && (mNk>1)) {
      return MOD_CYLINDER;
    } else if ((mNi==1) && (mNj==1) && (mNk>1)) {
      return MOD_SPHERESHELL;
    } else {
      return MOD_TETRAWCG;
    }
  }
}


//////
// METHOD:  Grid :: DumpGridToAscii()
//
void Grid::DumpGridToAscii() const {
  
  std::ostream & out = std::cout;
  std::ostringstream nodeprefix;

  out << "#  R3D_GRID:\n"
      << "#  Line 1:  ni nj nk\n"
      << "#  Line 2:  Index_Base\n"
      << "#  Lines 3 and up describe grid nodes:\n"
      << "#  i j k    x y x    vp vs rho qp qs    nu eps a k\n"
      << "#\n"
      << mNi << " " << mNj << " " << mNk << "\n"
      << 0 << "\n";
 
  for (Index k=0; k<mNk; k++) {
    for (Index j=0; j<mNj; j++) {
      for (Index i=0; i<mNi; i++) {
        nodeprefix.str(""); // Clear nodeprefix
        nodeprefix << std::right
                   << std::setw(3) << i << " " 
                   << std::setw(3) << j << " " 
                   << std::setw(3) << k << "  ";
        Node(i,j,k).OutputAsAscii(out, nodeprefix.str());
      }
    }
  }

  out << "#  END R3D_GRID\n";

}
