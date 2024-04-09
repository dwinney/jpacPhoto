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
#include "semi_inclusive/vector_exchange.hpp"
#include "semi_inclusive/photon_exchange.hpp"
#include "regge/vector_exchange.hpp"
#include "covariant/photon_exchange.hpp"

void inclusives()
{
    using namespace jpacPhoto;
    using inclusive::vector_exchange;
    using covariant::photon_exchange;

    //----------------------------------------------------------------------------
    // INPUTS

    // VMD couplings
    double gamma_omega = 56.34;
    double gamma_rho   = 16.37;
    double gamma_psi   = 36.85;
   
    // Form factor cutoffs
    double lamRho   = 1.4, lamOmega = 1.2;

    //----------------------------------------------------------------------------
    // chi_c1

    kinematics kC = new_kinematics(M_CHIC1);
    kC->set_meson_JP(AXIALVECTOR);

    std::vector<double> C_omega_pars = {10.46E-3, gamma_omega/2., lamOmega};
    std::vector<double> C_rho_pars   = {18.87E-3, gamma_rho/2.,   lamRho};
    std::vector<double> C_gamma_pars = {3.6E-2, 1, 0.};

    amplitude eC_omega = new_amplitude<photon_exchange>(kC, M_OMEGA, "Omega Exchange");
    eC_omega->set_parameters(C_omega_pars);

    amplitude eC_rho   = new_amplitude<photon_exchange>(kC, M_RHO, "Rho Exchange");
    eC_rho->set_parameters(C_rho_pars);

    amplitude eC_gam = new_amplitude<photon_exchange>(kC, 0., "#gamma exchange");
    eC_gam->set_parameters(C_gamma_pars);

    semi_inclusive iC_V = new_semi_inclusive<inclusive::vector_exchange>(kC, "Inclusive");
    iC_V->set_parameters({C_rho_pars[0], C_omega_pars[0]});

    semi_inclusive iC_gam   = new_semi_inclusive<inclusive::photon_exchange>(kC, 0, "Inclusive");
    iC_gam->reggeized(true);
    iC_gam->set_parameters(C_gamma_pars);

    amplitude      eC = eC_omega + eC_rho;
    semi_inclusive iC = iC_V + iC_gam + eC;

    // Reggeized exclusive amplitude
    amplitude rC_omega = new_amplitude<regge::vector_exchange>(kC, "#omega exchange");
    rC_omega->set_parameters({5.2E-4, 16., 0., 1.2, 0.5, 0.9});

    amplitude rC_rho = new_amplitude<regge::vector_exchange>(kC, "#rho exchange");
    rC_rho->set_parameters({9.2E-4, 2.4, 14.6, 1.4, 0.5, 0.9});

    amplitude rC = rC_omega + rC_rho;
    semi_inclusive irC = iC_V + iC_gam + rC;

    //---------------------------------------------------------------------------
    // X(3872)

    kinematics kX = new_kinematics(M_X3872);
    kX->set_meson_JP(AXIALVECTOR);

    std::vector<double> X_omega_pars = {0.199228,  gamma_omega/2., lamOmega};
    std::vector<double> X_rho_pars   = {0.0879857, gamma_rho/2.,   lamRho};
    std::vector<double> X_gamma_pars = {3.2E-3, 1, 0.};

    amplitude eX_omega = new_amplitude<photon_exchange>(kX, M_OMEGA, "Omega Exchange");
    eX_omega->set_parameters(X_omega_pars);

    amplitude eX_rho   = new_amplitude<photon_exchange>(kX, M_RHO, "Rho Exchange");
    eX_rho->set_parameters(X_rho_pars);

    amplitude eX_gam   = new_amplitude<photon_exchange>(kX, 0., "#gamma Exchange");
    eX_gam->set_parameters(X_gamma_pars);

    semi_inclusive iX_V = new_semi_inclusive<inclusive::vector_exchange>(kX, "Inclusive");
    iX_V->set_parameters({X_rho_pars[0], X_omega_pars[0]});

    semi_inclusive iX_gam = new_semi_inclusive<inclusive::photon_exchange>(kX, 0, "Inclusive");
    iX_gam->reggeized(true);
    iX_gam->set_parameters(X_gamma_pars);

    amplitude      eX = eX_omega + eX_rho;
    semi_inclusive iX = iX_V + eX;

    // Reggeized exclusive amplitude 
    amplitude rX_omega = new_amplitude<regge::vector_exchange>(kX, "#omega exchange");
    rX_omega->set_parameters({8.2E-3, 16., 0., lamOmega, 0.5, 0.9});

    amplitude rX_rho   = new_amplitude<regge::vector_exchange>(kX, "#rho exchange");
    rX_rho->set_parameters({3.6E-3, 2.4, 14.6, lamRho, 0.5, 0.9});

    amplitude rX = rX_omega + rX_rho;
    semi_inclusive irX = iX_V + rX;

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
    p1.add_curve( C_NT, [&](double W){ return iC->integrated_xsection(W*W, 0.7); }, "#chi_{#it{c}1}");
    print("chic1 (exclusive)"); divider(2);
    p1.add_dashed(C_NT, [&](double W){ return eC->integrated_xsection(W*W); });
    print("X(3872) (inclusive)"); divider(2);
    p1.add_curve( X_NT, [&](double W){ return iX->integrated_xsection(W*W, 0.7); }, "#it{X}(3872)");
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
    iC_V->reggeized(true); 
    p3.print_to_terminal(true);
    p3.add_curve( HE, [&](double W){ return (irC->integrated_xsection(W*W)+eC_gam->integrated_xsection(W*W)) * 1E3; }, "Total");
    p3.add_curve(HE, [&](double W){  return iC_V->integrated_xsection(W*W)  * 1E3; },   "Inclusive #it{V} exchange");
    p3.add_curve(HE, [&](double W){  return rC->integrated_xsection(W*W)  * 1E3; },     "Exclusive #it{V} exchange");
    p3.add_curve(HE, [&](double W){  return iC_gam->integrated_xsection(W*W)  * 1E3; }, "Inclusive #gamma exchange");
    p3.add_curve(HE, [&](double W){  return eC_gam->integrated_xsection(W*W)  * 1E3; }, "Exclusive #gamma exchange");

    // Plot the breakdown of contributions for the X(3872)
    plot p2 = plotter.new_plot();
    p2.set_curve_points(20);
    p2.set_logscale(false, true);
    p2.set_ranges(HE, {5E-4, 3E2});
    p2.set_labels( "#it{W}_{#gamma#it{p}}  [GeV]", "#sigma  [pb]");
    p2.set_legend(0.80, 0.17);
    p2.add_header("#it{X}(3872)");
    p2.print_to_terminal(true);
    iX_V->reggeized(true);
    p2.add_curve( HE, [&](double W){ return irX->integrated_xsection(W*W) * 1E3; });
    p2.add_curve( HE, [&](double W){ return iX_V->integrated_xsection(W*W) * 1E3; });
    p2.add_curve( HE, [&](double W){ return rX->integrated_xsection(W*W) * 1E3; });
    p2.add_curve( HE, [&](double W){ return iX_gam->integrated_xsection(W*W) * 1E3; });
    p2.add_curve( HE, [&](double W){ return eX_gam->integrated_xsection(W*W) * 1E3; });

    plotter.combine({2,1}, {p3,p2}, "HE.pdf");
};