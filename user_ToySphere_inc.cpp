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

  gr.SetSize(1,1,10);           // Sets index bounds
  gr.SetIndexBase(0);           // When addressing nodes, use base 0
  gr.SetMapping(Grid::GC_RAE, Grid::GC_SPHERICAL);

  gr.WNode(0,0,0).SetLocation ( 0, 0, 0 );
  gr.WNode(0,0,1).SetLocation ( 0, 0, -140.0 );
  gr.WNode(0,0,2).SetLocation ( 0, 0, -200.0 );   // Moho
  gr.WNode(0,0,3).SetLocation ( 0, 0, -600.0 );
  gr.WNode(0,0,4).SetLocation ( 0, 0, -2000.0 );
  gr.WNode(0,0,5).SetLocation ( 0, 0, -2890.0 );  // CMB
  gr.WNode(0,0,6).SetLocation ( 0, 0, -3500.0 );
  gr.WNode(0,0,7).SetLocation ( 0, 0, -4500.0 );
  gr.WNode(0,0,8).SetLocation ( 0, 0, -5150.0 );  // ICB
  gr.WNode(0,0,9).SetLocation ( 0, 0, -6371.0 );  // Center

  gr.WNode(0,0,0).SetAttributes( VpVs( 6.10, 2.60), 3.60, QCr, HSCr );
  gr.WNode(0,0,1).SetAttributes( VpVs( 7.20, 4.00), 4.00, QCr, HSCr );
  gr.WNode(0,0,2).SetAttributes( VpVs( 7.60, 4.40), 4.00, QCr, HSCr );

  gr.WNode(0,0,2).SetAttributes( VpVs( 8.00, 4.60), 4.20, QMa, HSMa );
  gr.WNode(0,0,3).SetAttributes( VpVs( 9.80, 6.00), 4.90, QMa, HSMa );
  gr.WNode(0,0,4).SetAttributes( VpVs(13.00, 7.20), 5.60, QMa, HSMa );
  gr.WNode(0,0,5).SetAttributes( VpVs(13.00, 7.20), 5.60, QMa, HSMa );

  gr.WNode(0,0,5).SetAttributes( VpVs( 8.00, 0.00), 10.0, QCo, HSCo );
  gr.WNode(0,0,6).SetAttributes( VpVs( 9.00, 0.00), 11.0, QCo, HSCo );
  gr.WNode(0,0,7).SetAttributes( VpVs(10.00, 0.00), 12.0, QCo, HSCo );
  gr.WNode(0,0,8).SetAttributes( VpVs(10.00, 0.00), 12.0, QCo, HSCo );

  gr.WNode(0,0,8).SetAttributes( VpVs(11.00, 3.50), 13.0, QCi, HSCi );
  gr.WNode(0,0,9).SetAttributes( VpVs(12.00, 3.60), 13.5, QCi, HSCi );


}
