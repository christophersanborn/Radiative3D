// user_ToySphere_inc.cpp
//
//   Defines a specific model.
//   #include this file in user.cpp
//

//////
// Toy Sphere Model
//
// Uses RAE (Range, Azimuth, Elevation) coordinate scheme in a spherical
// mapping.
//
void ToySphere(Grid & gr, const std::vector<Real> & args) {

  using Elastic::Velocity;
  using Elastic::VpVs;
  using Elastic::Q;
  using Elastic::QmQk;
  using Elastic::HetSpec;
  using Elastic::HSneak;

  Real sc_nu  = 0.8;    // Default scat params -
  Real sc_eps = 0.01;   // we'll scan args for user values
  Real sc_a   = 4.00;   //
  Real sc_k = 0.8;      //
  Real q = 1000;        //

  //HetSpec HSCr = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Crust
  HetSpec HSMa = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Mantle
  HetSpec HSCo = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Core Outer
  //HetSpec HSCi = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Core Inner


  //Q QCr = QmQk(q);          //
  Q QMa = QmQk(q);          //
  Q QCo = QmQk(q);          //
  //Q QCi = QmQk(q);          //

  gr.SetSize(1,1,3);            // Sets index bounds
  gr.SetIndexBase(0);           // When addressing nodes, use base 0
  gr.SetMapping(Grid::GC_RAE, Grid::GC_SPHERICAL);

  gr.WNode(0,0,0).SetLocation ( 0, 0, 0 );
  gr.WNode(0,0,1).SetLocation ( 0, 0, -4000.0 );  // "Core"
  gr.WNode(0,0,2).SetLocation ( 0, 0, -6371.0 );  // Center

  gr.WNode(0,0,0).SetAttributes( VpVs( 5.00, 2.60), 3.60, QMa, HSMa );
  gr.WNode(0,0,1).SetAttributes( VpVs( 9.00, 6.00), 4.00, QMa, HSMa );
  gr.WNode(0,0,1).SetAttributes( VpVs(10.00, 8.00), 4.20, QCo, HSCo );
  gr.WNode(0,0,2).SetAttributes( VpVs(14.00,12.00), 4.90, QCo, HSCo );


}
