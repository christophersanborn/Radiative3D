// user_SphereEarth_inc.cpp
//
//   Defines a specific model.
//   #include this file in user.cpp
//

//////
// Spherical Earth Model
//
// Uses RAE (Range, Azimuth, Elevation) coordinate scheme in a spherical
// mapping.
//
void SphereEarth(Grid & gr, const std::vector<Real> & args) {

  using Elastic::Velocity;
  using Elastic::VpVs;
  using Elastic::Q;
  using Elastic::QmQk;
  using Elastic::HetSpec;
  using Elastic::HSneak;

  Real sc_nu  = 0.8;    // Default scat params -
  Real sc_eps = 0.005;  // we'll scan args for user values
  Real sc_a   = 4.00;   //
  Real sc_k = 0.8;      //
  Real q = 2000;        //

  HetSpec HSCr = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Crust
  HetSpec HSMa = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Mantle
  HetSpec HSCo = HSneak(sc_nu, sc_eps, 2*sc_a, sc_k); // Core Outer
  HetSpec HSCi = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Core Inner


  Q QCr = QmQk(0.8*q);      //
  Q QMa = QmQk(q);          //
  Q QCo = QmQk(1,20*q);     //
  Q QCi = QmQk(q);          //

  gr.SetSize(1,1,16);           // Sets index bounds
  gr.SetIndexBase(0);           // When addressing nodes, use base 0
  gr.SetMapping(Grid::GC_RAE, Grid::GC_SPHERICAL);

  gr.WNode(0,0,0).SetLocation ( 0, 0, 0 );
  gr.WNode(0,0,1).SetLocation ( 0, 0, -100.0 );   // Moho
  gr.WNode(0,0,2).SetLocation ( 0, 0, -410.0 );
  gr.WNode(0,0,3).SetLocation ( 0, 0, -660.0 );
  gr.WNode(0,0,4).SetLocation ( 0, 0, -958.0 );
  gr.WNode(0,0,5).SetLocation ( 0, 0, -1354.0 );
  gr.WNode(0,0,6).SetLocation ( 0, 0, -2047.0 );
  gr.WNode(0,0,7).SetLocation ( 0, 0, -2789.0 );
  gr.WNode(0,0,8).SetLocation ( 0, 0, -2891.0 );  // CMB
  gr.WNode(0,0,9).SetLocation ( 0, 0, -3594.0 );
  gr.WNode(0,0,10).SetLocation ( 0, 0, -4298.0 );
  gr.WNode(0,0,11).SetLocation ( 0, 0, -4852.0 );
  gr.WNode(0,0,12).SetLocation ( 0, 0, -5153.0 );  // ICB
  gr.WNode(0,0,13).SetLocation ( 0, 0, -5661.0 );
  gr.WNode(0,0,14).SetLocation ( 0, 0, -6066.0 );
  gr.WNode(0,0,15).SetLocation ( 0, 0, -6371.0 );  // Center

  gr.WNode(0,0,0).SetAttributes( VpVs( 5.80, 3.20), 2.60, QCr, HSCr );
  gr.WNode(0,0,1).SetAttributes( VpVs( 6.80, 3.90), 2.92, QCr, HSCr );

  gr.WNode(0,0,1).SetAttributes( VpVs( 8.04, 4.48), 3.64, QMa, HSMa );
  gr.WNode(0,0,2).SetAttributes( VpVs( 9.03, 4.87), 3.51, QMa, HSMa );
  gr.WNode(0,0,2).SetAttributes( VpVs( 9.36, 5.08), 3.93, QMa, HSMa );
  gr.WNode(0,0,3).SetAttributes( VpVs(10.20, 5.61), 3.92, QMa, HSMa );
  gr.WNode(0,0,3).SetAttributes( VpVs(10.79, 5.96), 4.24, QMa, HSMa );
  gr.WNode(0,0,4).SetAttributes( VpVs(11.39, 6.35), 5.60, QMa, HSMa );
  gr.WNode(0,0,5).SetAttributes( VpVs(11.99, 6.60), 4.57, QMa, HSMa );
  gr.WNode(0,0,6).SetAttributes( VpVs(12.85, 6.94), 5.13, QMa, HSMa );
  gr.WNode(0,0,7).SetAttributes( VpVs(13.65, 7.26), 5.72, QMa, HSMa );
  gr.WNode(0,0,8).SetAttributes( VpVs(13.66, 7.28), 5.77, QMa, HSMa );

  gr.WNode(0,0,8).SetAttributes( VpVs( 8.00, 1e-5), 9.91, QCo, HSCo );
  gr.WNode(0,0,9).SetAttributes( VpVs( 9.08, 1e-5), 10.9, QCo, HSCo );
  gr.WNode(0,0,10).SetAttributes( VpVs(9.79, 1e-5), 11.6, QCo, HSCo );
  gr.WNode(0,0,11).SetAttributes( VpVs(10.17, 1e-5), 12.0, QCo, HSCo );
  gr.WNode(0,0,12).SetAttributes( VpVs(10.29, 1e-5), 12.1, QCo, HSCo );

  gr.WNode(0,0,12).SetAttributes( VpVs(11.04, 3.50), 12.7, QCi, HSCi );
  gr.WNode(0,0,13).SetAttributes( VpVs(11.18, 3.61), 12.9, QCi, HSCi );
  gr.WNode(0,0,14).SetAttributes( VpVs(11.25, 3.66), 13.0, QCi, HSCi );
  gr.WNode(0,0,15).SetAttributes( VpVs(11.26, 3.67), 13.0, QCi, HSCi );


}
