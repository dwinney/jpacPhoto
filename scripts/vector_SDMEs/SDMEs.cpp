// Model for high energy photoproduction of the vector mesons via Regge exchanges
//
// Adapted from models and codes originally written by V. Mathieu in [1]
// Plots SDMEs for rho, omega, and phi photoproduction. 
// Reproduces Fig. 8 in [1] into "sdmes.pdf"
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Univeritat Bonn (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// References:
// [1] - https://arxiv.org/abs/1802.09403
// ------------------------------------------------------------------------------


#include "plotter.hpp"
#include "regge/reggeon_exchange.hpp"
#include "amplitude.hpp"
#include "elementwise.hpp"

void SDMEs()
{
    using namespace jpacPhoto;
    using regge::reggeon_exchange;
    using complex = std::complex<double>;

    // Set up kinematics
    kinematics rho = new_kinematics(M_RHO);
    rho->set_meson_JP(VECTOR);

    kinematics omega = new_kinematics(M_OMEGA);
    omega->set_meson_JP(VECTOR);
    
    kinematics phi = new_kinematics(M_PHI);
    phi->set_meson_JP(VECTOR);

    //-------------------------------------------------
    // Unnatural exchanges
    
    double bU = 0.; // no form-factor 

    // Top couplings
    // Multiplied by masses and half to match normalization
    std::array<double,3> gT_pi  = {0.252*M_RHO*M_RHO/2, 0.696*M_OMEGA*M_OMEGA/2, 0.040*M_PHI*M_PHI/2};
    std::array<double,3> gT_eta = {0.136*M_RHO*M_RHO/2, 0.479*M_OMEGA*M_OMEGA/2, 0.210*M_PHI*M_PHI/2};

    // Bottom couplings
    double g_piNN = 13.26,  g_etaNN = 2.24;
    std::array<double, 2> betaB_pi  = {0., -2.*M_PROTON*g_piNN};  // Only flip allowed
    std::array<double, 2> betaB_eta = {0., -2.*M_PROTON*g_etaNN}; // Only flip allowed

    // Trajectory
    double alpha0U = -0.013, alphaPU = 0.7;

    // Amplitudes
    amplitude rho_pi  = new_amplitude<reggeon_exchange>(rho, -1, "#pi");
    rho_pi->set_parameters( {alpha0U, alphaPU, bU, 
                             gT_pi[0],   sqrt(2)*gT_pi[0],  gT_pi[0],  
                             betaB_pi[0],  betaB_pi[1]});

    amplitude rho_eta = new_amplitude<reggeon_exchange>(rho, -1, "#eta");
    rho_eta->set_parameters({alpha0U, alphaPU, bU, 
                             gT_eta[0],   sqrt(2)*gT_eta[0],  gT_eta[0],
                             betaB_eta[0], betaB_eta[1]});

    amplitude omega_pi  = new_amplitude<reggeon_exchange>(omega, -1, "#pi");
    omega_pi->set_parameters( {alpha0U, alphaPU, bU, 
                             gT_pi[1],   sqrt(2)*gT_pi[1],  gT_pi[1],  
                             betaB_pi[0],  betaB_pi[1]});

    amplitude omega_eta = new_amplitude<reggeon_exchange>(omega, -1, "#eta");
    omega_eta->set_parameters({alpha0U, alphaPU, bU, 
                             gT_eta[1],   sqrt(2)*gT_eta[1],  gT_eta[1],
                             betaB_eta[0], betaB_eta[1]});

    amplitude phi_pi  = new_amplitude<reggeon_exchange>(phi, -1, "#pi");
    phi_pi->set_parameters( {alpha0U, alphaPU, bU, 
                             gT_pi[2],   sqrt(2)*gT_pi[2],  gT_pi[2],  
                             betaB_pi[0],  betaB_pi[1]});

    amplitude phi_eta = new_amplitude<reggeon_exchange>(phi, -1, "#eta");
    phi_eta->set_parameters({alpha0U, alphaPU, bU, 
                             gT_eta[2],   sqrt(2)*gT_eta[2],  gT_eta[2],
                             betaB_eta[0], betaB_eta[1]});


    //-------------------------------------------------
    // Natural exchanges

    // Cutoffs
    double b_pom     = 3.60,  b_f2     = 0.55,  b_a2     = 0.53;

    // Trajectories
    double alp0_pom  = 1.08,  alp0_f2  = 0.5,   alp0_a2  = 0.5;
    double alpP_pom  = 0.2,   alpP_f2  = 0.9,   alpP_a2  = 0.9;

    // Overall normalization
    std::array<double,3> beta_pom = {2.506*2/alpP_pom/PI, 0.739*2/alpP_pom/PI, 0.932*2/alpP_pom/PI};
    std::array<double,3> beta_f2  = {2.476*2/alpP_f2/PI,  0.730*2/alpP_f2/PI,  0.                 }; 
    std::array<double,3> beta_a2  = {0.370*2/alpP_a2/PI,  1.256*2/alpP_a2/PI,  0.                 };

    // Relative couplings
    std::array<std::array<double,3>,3> rel;
    rel[0] = {1., 0.,             0.};
    rel[1] = {1., -0.95/sqrt(2.), -0.56};
    rel[2] = {1., -0.83/sqrt(2.),  0.} ;  

    // Bottom couplings 
    double kappa_pom = 0.,    kappa_f2 = 0,     kappa_a2 = -8.0;

    // Amplitudes
    amplitude rho_pomeron = new_amplitude<reggeon_exchange>(rho, +1, "Pomeron");
    rho_pomeron->set_parameters( {alp0_pom, alpP_pom, b_pom, 
                              beta_pom[0]*rel[0][0], beta_pom[0]*rel[0][1], beta_pom[0]*rel[0][2], 
                              1., kappa_pom} );

    amplitude rho_f2 = new_amplitude<reggeon_exchange>(rho, +1, "#it{f}_{2}");
    rho_f2->set_parameters( { alp0_f2, alpP_f2, b_f2, 
                              beta_f2[0]*rel[1][0], beta_f2[0]*rel[1][1], beta_f2[0]*rel[1][2], 
                              1., kappa_f2} );

    amplitude rho_a2 = new_amplitude<reggeon_exchange>(rho, +1, "#it{a}_{2}");
    rho_a2->set_parameters( { alp0_a2, alpP_a2, b_a2, 
                              beta_a2[0]*rel[2][0], beta_a2[0]*rel[2][1], beta_a2[0]*rel[2][2], 
                              1., kappa_a2} );

    amplitude omega_pomeron = new_amplitude<reggeon_exchange>(omega, +1, "Pomeron");
    omega_pomeron->set_parameters( {alp0_pom, alpP_pom, b_pom, 
                              beta_pom[1]*rel[0][0], beta_pom[1]*rel[0][1], beta_pom[1]*rel[0][2], 
                              1., kappa_pom} );

    amplitude omega_f2 = new_amplitude<reggeon_exchange>(omega, +1, "#it{f}_{2}");
    omega_f2->set_parameters( { alp0_f2, alpP_f2, b_f2, 
                              beta_f2[1]*rel[1][0], beta_f2[1]*rel[1][1], beta_f2[1]*rel[1][2], 
                              1., kappa_f2} );

    amplitude omega_a2 = new_amplitude<reggeon_exchange>(omega, +1, "#it{a}_{2}");
    omega_a2->set_parameters( { alp0_a2, alpP_a2, b_a2, 
                              beta_a2[1]*rel[2][0], beta_a2[1]*rel[2][1], beta_a2[1]*rel[2][2], 
                              1., kappa_a2} );

    amplitude phi_pomeron = new_amplitude<reggeon_exchange>(phi, +1, "Pomeron");
    phi_pomeron->set_parameters( {alp0_pom, alpP_pom, b_pom, 
                              beta_pom[2]*rel[0][0], beta_pom[2]*rel[0][1], beta_pom[2]*rel[0][2], 
                              1., kappa_pom} );    

    // //-------------------------------------------------
    // // Full models

    amplitude rho_total = rho_pomeron + rho_f2 + rho_a2 + rho_eta + rho_pi;
    rho_total->set_id("Total");

    amplitude omega_total = omega_pomeron + omega_f2 + omega_a2 + omega_eta + omega_pi;
    omega_total->set_id("Total");

    amplitude phi_total = phi_pomeron + phi_eta + phi_pi;
    phi_total->set_id("Total");

    //-------------------------------------------------
    // Plot results

    double s = s_cm(8.5);
    std::array<double,2> bounds = {0., 1.0};
    plotter plotter;

    plot p000 = plotter.new_plot();
    p000.set_ranges(bounds, {-0.35, 0.35});
    p000.set_labels("#minus #it{t}  [GeV^{2}]", "#rho^{0}_{00}");
    p000.add_curve( bounds, [&](double mt){ return  std::real( rho_total->mSDME_H(0, 0, 0, s, -mt) ); }, "#rho");
    p000.add_curve( bounds, [&](double mt){ return  std::real( phi_total->mSDME_H(0, 0, 0, s, -mt) ); }, "#phi");
    p000.add_curve( bounds, [&](double mt){ return  std::real( omega_total->mSDME_H(0, 0, 0, s, -mt) ); }, "#omega");
    p000.set_legend(0.7,0.3);

    plot p010 = plotter.new_plot();
    p010.set_ranges(bounds, {-0.35, 0.35});
    p010.set_labels("#minus #it{t}  [GeV^{2}]", "Re #rho^{0}_{10}");
    p010.add_curve( bounds, [&](double mt){ return  std::real( rho_total->mSDME_H(0, 1, 0, s, -mt) ); });
    p010.add_curve( bounds, [&](double mt){ return  std::real( phi_total->mSDME_H(0, 1, 0, s, -mt) ); });
    p010.add_curve( bounds, [&](double mt){ return  std::real( omega_total->mSDME_H(0, 1, 0, s, -mt) ); });

    plot p01m1 = plotter.new_plot();
    p01m1.set_ranges(bounds, {-0.35, 0.35});
    p01m1.set_labels("#minus #it{t}  [GeV^{2}]", "#rho^{0}_{1-1}");
    p01m1.add_curve( bounds, [&](double mt){ return  std::real( rho_total->mSDME_H(0, 1, -1, s, -mt) ); });
    p01m1.add_curve( bounds, [&](double mt){ return  std::real( phi_total->mSDME_H(0, 1, -1, s, -mt) ); });
    p01m1.add_curve( bounds, [&](double mt){ return  std::real( omega_total->mSDME_H(0, 1, -1, s, -mt) ); });

    plot p111 = plotter.new_plot();
    p111.set_ranges(bounds, {-0.35, 0.35});
    p111.set_labels("#minus #it{t}  [GeV^{2}]", "#rho^{1}_{11}");
    p111.add_curve( bounds, [&](double mt){ return  std::real( rho_total->mSDME_H(1, 1, 1, s, -mt) ); });
    p111.add_curve( bounds, [&](double mt){ return  std::real( phi_total->mSDME_H(1, 1, 1, s, -mt) ); });
    p111.add_curve( bounds, [&](double mt){ return  std::real( omega_total->mSDME_H(1, 1, 1, s, -mt) ); });

    plot p100 = plotter.new_plot();
    p100.set_ranges(bounds, {-0.35, 0.35});
    p100.set_labels("#minus #it{t}  [GeV^{2}]", "#rho^{1}_{00}");
    p100.add_curve( bounds, [&](double mt){ return  std::real( rho_total->mSDME_H(1, 0, 0, s, -mt) ); });
    p100.add_curve( bounds, [&](double mt){ return  std::real( phi_total->mSDME_H(1, 0, 0, s, -mt) ); });
    p100.add_curve( bounds, [&](double mt){ return  std::real( omega_total->mSDME_H(1, 0, 0, s, -mt) ); });
    
    plot p110 = plotter.new_plot();
    p110.set_ranges(bounds, {-0.35, 0.35});
    p110.set_labels("#minus #it{t}  [GeV^{2}]", "Re #rho^{1}_{10}");
    p110.add_curve( bounds, [&](double mt){ return  std::real( rho_total->mSDME_H(1, 1, 0, s, -mt) ); });
    p110.add_curve( bounds, [&](double mt){ return  std::real( phi_total->mSDME_H(1, 1, 0, s, -mt) ); });
    p110.add_curve( bounds, [&](double mt){ return  std::real( omega_total->mSDME_H(1, 1, 0, s, -mt) ); });

    plot p11m1 = plotter.new_plot();
    p11m1.set_ranges(bounds, {-0.6, 0.6});
    p11m1.set_labels("#minus #it{t}  [GeV^{2}]", "#rho^{1}_{1-1}");
    p11m1.add_curve( bounds, [&](double mt){ return  std::real( rho_total->mSDME_H(1, 1, -1, s, -mt) ); });
    p11m1.add_curve( bounds, [&](double mt){ return  std::real( phi_total->mSDME_H(1, 1, -1, s, -mt) ); });
    p11m1.add_curve( bounds, [&](double mt){ return  std::real( omega_total->mSDME_H(1, 1, -1, s, -mt) ); });

    plot p210 = plotter.new_plot();
    p210.set_ranges(bounds, {-0.6, 0.6});
    p210.set_labels("#minus #it{t}  [GeV^{2}]", "Im #rho^{2}_{10}");
    p210.add_curve( bounds, [&](double mt){ return  -std::imag( rho_total->mSDME_H(2, 1, 0, s, -mt) ); });
    p210.add_curve( bounds, [&](double mt){ return  -std::imag( phi_total->mSDME_H(2, 1, 0, s, -mt) ); });
    p210.add_curve( bounds, [&](double mt){ return  -std::imag( omega_total->mSDME_H(2, 1, 0, s, -mt) ); });

    plot p21m1 = plotter.new_plot();
    p21m1.set_ranges(bounds, {-0.6, 0.6});
    p21m1.set_labels("#minus #it{t}  [GeV^{2}]", "Im #rho^{2}_{1-1}");
    p21m1.add_curve( bounds, [&](double mt){ return  -std::imag( rho_total->mSDME_H(2, 1, -1, s, -mt) ); });
    p21m1.add_curve( bounds, [&](double mt){ return  -std::imag( phi_total->mSDME_H(2, 1, -1, s, -mt) ); });
    p21m1.add_curve( bounds, [&](double mt){ return  -std::imag( omega_total->mSDME_H(2, 1, -1, s, -mt) ); });

    plotter.combine({3,3}, {p000, p010, p01m1, p111, p100, p110, p11m1, p210, p21m1}, "sdmes.pdf");
};