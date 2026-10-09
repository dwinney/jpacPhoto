// Calculates the integrated cross section for both inclusive and exclusive 
// axial-vector production at near threshold and high energies
// Reproduces figs. 8 and 8 of [1]
//
// OUTPUT: NT.pdf (fig 8)
//         HE.pdf (fig 9)
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@iu.alumni.edu
//               winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------
// REFERENCES:
//
// [1] - https://arxiv.org/abs/2404.05326
// ------------------------------------------------------------------------------

#include "plotter.hpp"
#include "inclusive/vector_exchange.hpp"
#include "regge/vector_exchange.hpp"
#include "covariant/photon_exchange.hpp"

void inclusives()
{
    using namespace jpacPhoto;
    using namespace jpacPhoto::inclusive;
    using covariant::photon_exchange;

    //----------------------------------------------------------------------------
    // INPUTS

    // VMD proportionality couplings
    double gamma_omega = 56.34;
    double gamma_rho   = 16.37;
    double gamma_psi   = 36.85;

    // Actual top coupling
    double gC_gamma = 3.6E-2, gC_rho = 18.87E-3,  gC_omega = 10.46E-3;
    double gX_gamma = 3.2E-3, gX_rho = 0.0879857, gX_omega = 0.199228;    
   
    // Form factor cutoffs
    double lamRho   = 1.4, lamOmega = 1.2;

    //----------------------------------------------------------------------------
    // chi_c1

    kinematics kC = new_kinematics(M_CHIC1);
    kC->set_meson_JP(AXIALVECTOR);

    std::vector<double> C_omega_pars = {M_OMEGA, gC_omega, gamma_omega/2., lamOmega};
    std::vector<double> C_rho_pars   = {M_RHO,   gC_rho, gamma_rho/2.,   lamRho};
    std::vector<double> C_gamma_pars = {0.,      gC_gamma,   1,              0.};

    amplitude eC_omega = new_amplitude<photon_exchange>(kC);
    eC_omega->set_option(photon_exchange::kVMD);
    eC_omega->set_parameters(C_omega_pars);

    amplitude eC_rho   = new_amplitude<photon_exchange>(kC);
    eC_rho->set_option(photon_exchange::kVMD);
    eC_rho->set_parameters(C_rho_pars);

    amplitude eC_gam = new_amplitude<photon_exchange>(kC);
    eC_gam->set_option(photon_exchange::kVMD);
    eC_gam->set_parameters(C_gamma_pars);

    semi_inclusive iC = new_semi_inclusive<inclusive::vector_exchange>(kC);
    iC->set_parameters({gC_gamma, gC_rho, gC_omega});

    amplitude eC = eC_omega + eC_rho + eC_gam;
    semi_inclusive intC = iC + eC;

    // Reggeized exclusive amplitude
    amplitude rC_omega = new_amplitude<regge::vector_exchange>(kC);
    rC_omega->set_parameters({0.5, 0.9, 5.2E-4, 16., 0., 1.2});

    amplitude rC_rho = new_amplitude<regge::vector_exchange>(kC);
    rC_rho->set_parameters({0.5, 0.9, 9.2E-4, 2.4, 14.6, 1.4});

    amplitude rC = rC_omega + rC_rho;
    semi_inclusive irC = iC + rC;

    // //---------------------------------------------------------------------------
    // // X(3872)

    kinematics kX = new_kinematics(M_X3872);
    kX->set_meson_JP(AXIALVECTOR);

    std::vector<double> X_omega_pars = {M_OMEGA, gX_omega,  gamma_omega/2., lamOmega};
    std::vector<double> X_rho_pars   = {M_RHO,   gX_rho, gamma_rho/2.,   lamRho};
    std::vector<double> X_gamma_pars = {0., gX_gamma, 1, 0.};

    amplitude eX_omega = new_amplitude<photon_exchange>(kX);
    eX_omega->set_option(photon_exchange::kVMD);
    eX_omega->set_parameters(X_omega_pars);

    amplitude eX_rho   = new_amplitude<photon_exchange>(kX);
    eX_rho->set_option(photon_exchange::kVMD);
    eX_rho->set_parameters(X_rho_pars);

    amplitude eX_gam   = new_amplitude<photon_exchange>(kX);
    eX_gam->set_option(photon_exchange::kVMD);
    eX_gam->set_parameters(X_gamma_pars);

    semi_inclusive iX = new_semi_inclusive<inclusive::vector_exchange>(kX);
    iX->set_parameters({gX_gamma, gX_rho, gX_omega});

    amplitude eX = eX_omega + eX_rho;
    semi_inclusive intX = iX + eX;

    // Reggeized exclusive amplitude 
    amplitude rX_omega = new_amplitude<regge::vector_exchange>(kX);
    rX_omega->set_parameters({8.2E-3, 16., 0., lamOmega, 0.5, 0.9});

    amplitude rX_rho   = new_amplitude<regge::vector_exchange>(kX);
    rX_rho->set_parameters({3.6E-3, 2.4, 14.6, lamRho, 0.5, 0.9});

    amplitude rX = rX_omega + rX_rho;
    semi_inclusive irX = iX + rX;

    // --------------------------------------------------------------------------
    // Plot results

    // Bounds to plot
    std::array<double,2> X_NT = {kX->Wth() + EPS, 7.0};
    std::array<double,2> C_NT = {kC->Wth() + EPS, 7.0};
    std::array<double,2> HE   = {20, 60};

    plotter plotter;

    // Near threshold production plot
    plot p1 = plotter.new_plot();
    p1.set_curve_points(30);
    p1.set_logscale(false, true);
    p1.set_ranges({4.1, 7}, {1E-2, 2E3});
    p1.set_legend(0.27, 0.72);
    p1.set_labels( "#it{W}_{#gamma#it{p}}  [GeV]", "#sigma  [nb]");
    p1.print_to_terminal(true);
    p1.shade_region({W_cm(22), 10}, {kBlack, 1001});
    print("chic1 (inclusive)"); divider(2);
    p1.add_curve( C_NT, [&](double W){ return intC->integrated_xsection(W*W, 0.7); }, "#chi_{#it{c}1}");
    print("chic1 (exclusive)"); divider(2);
    p1.add_dashed(C_NT, [&](double W){ return eC->integrated_xsection(W*W); });
    print("X(3872) (inclusive)"); divider(2);
    p1.add_curve( X_NT, [&](double W){ return intX->integrated_xsection(W*W, 0.7); }, "#it{X}(3872)");
    print("X(3872) (exclusive)"); divider(2);
    p1.add_dashed(X_NT, [&](double W){ return eX->integrated_xsection(W*W); });
    p1.save("NT.pdf");

    // // Plot the breakdown of contributions for the chic1
    plot p3 = plotter.new_plot();
    p3.set_curve_points(20);
    p3.set_logscale(false, true);
    p3.set_ranges(HE, {5E-4, 3E2});
    p3.set_labels( "#it{W}_{#gamma#it{p}}  [GeV]", "#sigma  [pb]");
    p3.set_legend(0.20, 0.17);
    p3.add_header("#chi_{c1}(1#it{P})");
    iC->set_option(vector_exchange::kReggeized); 
    p3.print_to_terminal(true);
    p3.add_curve( HE, [&](double W){ return (irC->integrated_xsection(W*W)+eC_gam->integrated_xsection(W*W)) * 1E3; }, "Total");
    iC->set_parameters({0, gC_rho, gC_omega});
    p3.add_curve(HE, [&](double W){  return iC->integrated_xsection(W*W)     * 1E3; }, "Inclusive #it{V} exchange");
    p3.add_curve(HE, [&](double W){  return rC->integrated_xsection(W*W)     * 1E3; }, "Exclusive #it{V} exchange");
    iC->set_parameters({gC_gamma, 0., 0.});
    p3.add_curve(HE, [&](double W){  return iC->integrated_xsection(W*W)     * 1E3; }, "Inclusive #gamma exchange");
    p3.add_curve(HE, [&](double W){  return eC_gam->integrated_xsection(W*W) * 1E3; }, "Exclusive #gamma exchange");

    // Plot the breakdown of contributions for the X(3872)
    plot p2 = plotter.new_plot();
    p2.set_curve_points(20);
    p2.set_logscale(false, true);
    p2.set_ranges(HE, {5E-4, 3E2});
    p2.set_labels( "#it{W}_{#gamma#it{p}}  [GeV]", "#sigma  [pb]");
    p2.set_legend(0.80, 0.17);
    p2.add_header("#it{X}(3872)");
    p2.print_to_terminal(true);
    iX->set_option(vector_exchange::kReggeized);
    p2.add_curve( HE, [&](double W){ return (irX->integrated_xsection(W*W)+eX_gam->integrated_xsection(W*W)) * 1E3; });
    iX->set_parameters({0, gX_rho, gX_omega});
    p2.add_curve( HE, [&](double W){ return iX->integrated_xsection(W*W)     * 1E3; });
    p2.add_curve( HE, [&](double W){ return rX->integrated_xsection(W*W)     * 1E3; });
    iX->set_parameters({gX_gamma, 0., 0.});
    p2.add_curve( HE, [&](double W){ return iX->integrated_xsection(W*W)     * 1E3; });
    p2.add_curve( HE, [&](double W){ return eX_gam->integrated_xsection(W*W) * 1E3; });

    plotter.combine({2,1}, {p3,p2}, "HE.pdf");
};