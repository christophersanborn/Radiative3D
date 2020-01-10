// grid.hpp
//
#ifndef GRID_H_
#define GRID_H_
//
#include <vector>
#include <ostream>  /* used in OutputAsAscii input args */
#include <string>
#include "ecs.hpp"      /* includes elastic.hpp and geom.hpp */


// In this file:
// CLASS DEFINITIONS FOR:
//
//   o  Class GridData
//   o  Class GridNode
//   o  Class Grid
//
// Search on "CLASS:" to jump between class definitions in this 
// file.
//

// *** FORWARD DECLARATIONS:
//

// *** TYPEDEFS:
//

// *** CLASS DEFINITIONS:
//

//////
// CLASS:   ::: GridData :::
//
//   Encapsulates the set of model data (velocities, scattering
//   parameters, etc.) known to exist at a given grid location (ie, at
//   a GridNode object).  The GridNode object will "own" one or more
//   of these GridData objects (ie, it will have two if a Node has a
//   different velocity above vs below a discontinuous boundary).
//
//   Derives from Elastic::HElastic, which encapsulates more-or-less
//   equivalent data.  Inherits member access functions Vp(), Vs(),
//   Rho(), Qp(), Qs(), getV(), getHS().
//
class GridData : public Elastic::HElastic {
public:

  GridData(Elastic::Velocity vel, Real density,
           Elastic::Q q, Elastic::HetSpec hs)
    : Elastic::HElastic(vel,density,q,hs) {
  }

  GridData(const Elastic::HElastic & he) :
    Elastic::HElastic(he)   // Copy construct from HElastic base
  {}

  GridData() : Elastic::HElastic(
                 Elastic::VpVs(0,0),0,
                 Elastic::Qinf(),Elastic::HetSpec())
  {}

};


//////
// CLASS:   ::: GridNode :::
//
//   Encapsulates a single node in the warped cartesian grid that
//   defines the model.  An array of these objects will be kept to
//   give the grid its regular (indicial) structure.  The GridNode
//   object will keep all needed info to define the location of the
//   node, the physical properties at that location, and any needed
//   relational information to facilitate building the model.
//
//   Grid "attributes" (referring to material properties) specified at
//   a particular grid node are meant to indicate knowledge of those
//   attributes at a particular point in space, and thus are expected
//   to affect all MediumCells that touch the node. (The MediumCells
//   may use such knowledge of their vertices to model, e.g.,
//   gradients of these properties, or they may define a stairstep
//   structure, giving precedence to certain nodes according to what
//   makes sence for a given modelling situation.  It is the
//   ultimately choice of the particular MediumCell.)
//
//   Because there may be a need to model discontinuous jumps in
//   material properties, (e.g. a Moho discontinuity), a mechanism is
//   provided to signal those discontinuities, even in the case where
//   MediumCells would interpolate the properties continuously between
//   nodes.  The mechanism is to allow up to TWO sets of attributes to
//   be defined on a given grid node, with the optional second
//   attribute set specifying properties immediately "below" the grid
//   node, and the first set then representing properties
//   "above". (Ordinarily the singular first set would apply both
//   above and below.)  This allows for discontinuities to be
//   specified along lateral interfaces within the Earth model.
//   (There is currently no mechanism to provide discontinuities along
//   vertical interfaces, but these would be rarely encountered in the
//   type of Earth modeling that we envision this software being used
//   for.)
//
class GridNode {
public: 

  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Class-Specific Enum Constants  (GridNode Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::

  enum layers_e { GN_ABOVE,     // ID's the "layer" above the node
                  GN_BELOW,     // ID's the "layer" below the node
                  GN_NLAY }; 

protected:
  ;
  // :::::::::::::::::::::::::::::::::::::
  // ::: Member Data  (GridNode Class) :::
  // :::::::::::::::::::::::::::::::::::::

  EarthCoords::Generic  mLocation;  // Location of node in generic (user
                                    // chosen) coordinates.  (ECS will used
                                    // later to map these onto model space
                                    // (XYZ) coords.)
                                    //
  layers_e  mNextData;        // Indicates next GridData slot to be populated
                              // by the SetAttributes function. Also serves
                              // as an indicator of how many attribute sets
                              // have been filled.
  GridData  mData[GN_NLAY];   // Placeholders for two "layers" of data. (If
                              // both have been filled, then this node sits
                              // on a first-order discontinuity.)

public:
  ;
  // ::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (GridNode Class) :::
  // ::::::::::::::::::::::::::::::::::::::

  GridNode() : 
    mLocation (  0,0,0   ),
    mNextData ( GN_ABOVE )
    // mData[] will be default constructed.
  {}
  
  //GridNode(const GridNode & obj);
  // NOTE: Copy-constructor not needed yet - object is flat, so
  // default will do. But there are places in the code where GridNodes
  // are copied, so if this object starts using dynamic allocation,
  // (and it probably will when we add handling for discontinuities
  // across more than just the horizontal interface), then a proper
  // copy-constructor WILL be needed.

  ~GridNode() {
  }


  // ::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property Set Methods  (GridNode Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::

  void SetLocation(Real x, Real y, Real z);     // Effective only once, ignored
  void SetLocation(EarthCoords::Generic xyz) {  // silently if called more than
    SetLocation(xyz.x1(), xyz.x2(), xyz.x3());  // once.
  }

  void AdjustLocation(Real x, Real y, Real z);  // Add x,y,z values to the
                                                // current location in user
                                                // coordinates.

  void SetAttributes(GridData gdata);
            // Can be called up to twice per grid node.  Second call differ-
            // entiates layer below fr/ layer above. This version takes full
            // constructed GridData object, but user will usually prefer one
            // of the forms below specifying the attributes explicitly.
  void SetAttributes(Elastic::Velocity vel, Real density,
                     Elastic::Q q, Elastic::HetSpec hs)
             {SetAttributes(GridData(vel,density,q,hs));}
  void SetAttributes(Real density, Elastic::Velocity vel,
                     Elastic::Q q, Elastic::HetSpec hs)
             {SetAttributes(GridData(vel,density,q,hs));}

  void ClearAttributes() {        // If we want to replace the attributes
    mNextData = GN_ABOVE;         // of a node, ClearAttributes() and set
  }                               // them again.

  // ::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (GridNode Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::

  R3::XYZ Loc() const;    // Returns location of node in XYZ model space
                          // (implies a conversion from user coordinate
                          // space, handled by ECS object.)

  EarthCoords::Generic GetRawLoc() const {  // Return location raw value
    return mLocation;                       // without ECS conversion.
  }

  GridData Data(layers_e lay) const;
  //                              Returns selected GridData object
  //                              (representing attributes "just
  //                              above" or "just below" the node).

  bool IsDiscontinuous() const {    // Return true/false as to
    return (mNextData == GN_NLAY);  // whether this node sits on a
  }                                 // first-order discontinuty.


  // :::::::::::::::::::::::::::::::::::::::::::
  // ::: Data-Dump Methods  (GridNode Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::

  void OutputAsAscii(std::ostream & out, const std::string & prefix) const;


};//
///
//


//////
// CLASS:   ::: Grid :::
//
//   Encapsulates an array of GridNode objects and supporting
//   metadata.  Grid nodes are currently indexed in three-dimensional
//   regular IJK fashion, though other indexing schemes may be
//   suported in the future to support more flexible modeling
//   schemes. Currently two types of modeling schemes are supported by
//   the IJK grids.  These are Warped Cartesian Grid models and
//   layered "Cylinder" models.
//
//   GRID INDEXING: Node access is by i,j,k indices. The relationship
//   between the ijk indices and spatial x,y,z directions is up to the
//   user, but in general the i and j indices will correspond in some
//   way to lateral extent and the k index with depth.  Correlating k
//   with depth is a somewhat firm constraint, insofar as whenever a
//   grid node has dual attribute specifications, (used to signify a
//   sharp change or discontinuity of material properties across an
//   interface), the interpretation is that the first value applies
//   just "above" the node, (i.e. in the smaller-k direction), and
//   that the second value set applies just "below" the node, (i.e. on
//   the larger-k side).  An illustration of node arrangement and
//   indexing and an example relationship to XYZ space is as follows:
//
//     ijk
//     000-------010-------020-------                                         .
//      |\         \         \               o----> (+j) Northwards (typical) .
//      | \         \         \              |\                               .
//      |  \         \         \             | \                              .
//      |  100-------110-------120-------    |  (+i) Eastwards (typical)      .
//      |   |\        |\        |\           v                                .
//     001  | \       | \       | \          (+k) Downwards (typical)         .
//      |\  |  \      |  \      |  \                                          .
//      | \ |         |         |                                             .
//      |  \|         |         |                   E.g.:  +i  -->  +x        .
//      |  101-------111-------121-------                  +j  -->  +y        .
//      |   |\        |\        |\                         +k  -->  -z        .
//          | \       | \       | \                                           .
//
//   Again, though, the relationship between ijk and xyz is COMPLETELY
//   a user choice (they will set an XYZ value on each node), and so
//   the matchup shown here is totally for example purposes.
//
class Grid {

public:
  ;
  // ::::::::::::::::::::::::::::::::::
  // ::: Member Types  (Grid Class) :::
  // ::::::::::::::::::::::::::::::::::


  enum model_target_e {   // Model Target - i.e., What type of model
                          // does the grid describe:
                          //
    MOD_AUTO,             //  Unknown - guess from grid grid dims.
    MOD_CYLINDER,         //  Model made up of stacked cylinders.
    MOD_TETRAWCG,         //  Model made of tetrahedra packed in regular
                          //   array following a warped cartesian grid (WCG).
    MOD_SPHERESHELL       //  Model made up of concentric spherical
                          //   shells.
  };

  enum gs_coords_e {      // Grid-Space Coordinates Choice:
                          //
    GC_ENU,               //   (East,North,Up)
    GC_RAE,               //   (Range, Azimuth, Elevation)
    GC_LLE,               //   (Lat, Lon, Elevation)
                          //
  };

  enum curvature_e {      // Curvature model to apply:
                          // (Compatible coords systems listed in parantheses)
                          //
    GC_ORTHO,             //   No curvature.  (XYZ, RAE)
    GC_FLATTENED,         //   Curvature simulated through depth
                          //    transformation.  (XYZ, RAE)
                          //
    GC_CURVED,            //   Curvilinear. Use GC_CURVED when it's important
    GC_SPHERICAL          //    that ICS origin be at ECS origin, e.g. when
                          //    doing upper-Earth models with WCG's.  Use
                          //    GC_SPHERICAL when it's important that ICS
                          //    origin is at Earth Center, like when doing
                          //    whole-Earth models with spherical
                          //    shells.  (RAE, LLE)
  };                      //


private:
  ;
  // :::::::::::::::::::::::::::::::::
  // ::: Member Data  (Grid Class) :::
  // :::::::::::::::::::::::::::::::::

  model_target_e            // Specifies the target model type that the
    mModelTarget;           // grid is supposed to define.  Used in
                            // Model() constructor to determine which
                            // build-out method to call.

  Count      mN;            // Total number of Nodes in the Grid
  Count     mNi;            // Number of nodes laterally, (i)
  Count     mNj;            // Number of nodes laterally, (j)
  Count     mNk;            // Number of nodes vertically (depth)
                            // (mN = mNi * mNj * mNk)

  Index    mIndexBase;      // Index base for addressing grid nodes (0 or
                            // 1). Affects grid building (the Set_
                            // functions) but not grid reading.  This is a
                            // convenience to allow modellers to use a base
                            // one index sytem if they choose.

  std::vector<GridNode> mGrid; // The node array


public:

  Grid();                   // Constructor.  Takes no args.  Grids are a
                            // "construct first, build later" kinda thing.

  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Set Methods  (Grid Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::

  void SetSize(Count ni, Count nj, Count nk);
  void SetIndexBase(Index base);
  void SetMapping(gs_coords_e, curvature_e);

  GridNode & WNode(Index i, Index j, Index k);
                    // Returns a writable reference to the node indicated
                    // by the indices. Index-base adjustment is applied to
                    // the indices.


  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (Grid Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::

  model_target_e GetModelType() const;
                    // Returns a value that the Model() class can use to
                    // decide which build-out method should be used to build
                    // the model from the grid.  NOTE: Returns a definitive
                    // value; will NOT return MOD_AUTO.  If mModelType is
                    // MOD_AUTO internally, then this method will guess the
                    // intended target by analyzing grid structure.

  const GridNode & Node(Index i, Index j, Index k) const;
                    // Returns a constant reference to the node indicated
                    // by the indices. Index-base adjustment is NOT applied
                    // to the indices. (Uses base 0 indexing.)

  const GridNode & RelNode(Index i, Index j, Index k,
                           RelIndex ri, RelIndex rj, RelIndex rk) const;
                    // Returns a constant reference to the node indicated by
                    // the indices, as adjusted by the relative indices.  Used
                    // to select adjacent nodes, as in when selecting eight
                    // nodes that form a "block" in a Warped Cartesian grid
                    // structure.  Uses base 0 indexing.
 
  Count N() const {return mN;}
  Count Ni() const {return mNi;}
  Count Nj() const {return mNj;}
  Count Nk() const {return mNk;}


  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Canned (Pre-Made) Grid Methods  (Grid Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  //
  //                    These provide a few built-in, ready-to-go grid
  //                    designs, mainy to use for testing purposes.
  //

  void ConstructGrid_UnitTest_Cylinder();
                        // A simplified cylinder model for testing
                        // purposes.  (NOT WRITTEN YET)

  void ConstructGrid_UnitTest_Tetra();
                        // A simplified tetra model for testing
                        // purposes.  (NOT WRITTEN YET)


  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: User-Written Methods  (Grid Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::

  void ConstructGridManual(int Selection, const std::vector<Real> &args);  
                               // Define this function, e.g. in main.cpp,
                               // if you want to use a manully-created
                               // compiled-in grid. (Used for testing.
                               // Final code will read grid from a file.)


  // :::::::::::::::::::::::::::::::::::::::
  // ::: Data-Dump Methods  (Grid Class) :::
  // :::::::::::::::::::::::::::::::::::::::

  void DumpGridToAscii() const;


};


////
///
#endif //#ifndef MODEL_H_
//
