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
  Real q = 200;         //

  HetSpec HSCr = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Crust
  HetSpec HSMa = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Mantle
  HetSpec HSCo = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Core Outer
  HetSpec HSCi = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Core Inner


  Q QCr = QmQk(q);          //
  Q QMa = QmQk(q);          //
  Q QCo = QmQk(q);          //
  Q QCi = QmQk(q);          //

  gr.SetSize(1,1,5);            // Sets index bounds
  gr.SetIndexBase(0);           // When addressing nodes, use base 0
  gr.SetMapping(Grid::GC_RAE, Grid::GC_SPHERICAL);

  gr.WNode(0,0,0).SetLocation ( 0, 0, 0 );
  gr.WNode(0,0,1).SetLocation ( 0, 0, -200.0 );
  gr.WNode(0,0,2).SetLocation ( 0, 0, -2890.0 );
  gr.WNode(0,0,3).SetLocation ( 0, 0, -5150.0 );
  gr.WNode(0,0,4).SetLocation ( 0, 0, -6371.0 );

  gr.WNode(0,0,0).SetAttributes( VpVs( 4.50, 2.60), 1.20, QCr, HSCr );
  gr.WNode(0,0,1).SetAttributes( VpVs( 9.00, 5.50), 2.20, QMa, HSMa );
  gr.WNode(0,0,2).SetAttributes( VpVs(14.00, 7.90), 5.60, QMa, HSMa );

  gr.WNode(0,0,2).SetAttributes( VpVs( 8.00, 0.00), 10.0, QCo, HSCo );
  gr.WNode(0,0,3).SetAttributes( VpVs(10.00, 0.00), 12.0, QCo, HSCo );

  gr.WNode(0,0,3).SetAttributes( VpVs(11.00, 3.50), 13.0, QCi, HSCi );
  gr.WNode(0,0,4).SetAttributes( VpVs(12.00, 3.60), 13.5, QCi, HSCi );


}
