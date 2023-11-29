// Calculates the integrated cross section for both inclusive and exclusive 
// X(3872) via photon and vector meson exchanges
//
// OUTPUT: inclusive_X.pdf
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#include "plotter.hpp"
#include "inclusive/vector_exchange.hpp"
#include "analytic/vector_exchange.hpp"
#include "covariant/photon_exchange.hpp"

void NT_X3872()
{
    using namespace jpacPhoto;

    const int p = inclusive::vector_exchange::kProton;
    const int n = inclusive::vector_exchange::kNeutron;

    double Wth = M_X3872 + M_PROTON;
    
    //----------------------------------------------------------------------------
    // INPUTS

    // Couplings
    double gGamma  = 3.20E-3;
    double gRho    = 5.38E-3;
    double gOmega  = 3.54E-3;

    // VMD factors
    double etaRho   = 16.37;
    double etaOmega = 56.34;

    // Form factor cutoffs
    double lamRho   = 1.4;
    double lamOmega = 1.2;

    std::vector<double> parsGamma, parsRho, parsOmega;
    parsGamma = {gGamma, 1,        0       };
    parsRho   = {gRho,   etaRho,   lamRho  };
    parsOmega = {gOmega, etaOmega, lamOmega};

    //----------------------------------------------------------------------------
    // Set up inclusive amplitudes

    inclusive_process inc_gamma = new_inclusive_process<inclusive::vector_exchange>(M_X3872, "Inclusive");
    inc_gamma->set_parameters(parsGamma);

    inclusive_process inc_omega = new_inclusive_process<inclusive::vector_exchange>(M_X3872, M_OMEGA, "Inclusive");
    inc_omega->set_parameters(parsOmega);

    inclusive_process inc_rho   = new_inclusive_process<inclusive::vector_exchange>(M_X3872, M_RHO, "Inclusive");
    inc_rho->set_parameters(parsRho);

    //----------------------------------------------------------------------------
    // Set up exclusives amplitudes

    kinematics kX = new_kinematics(M_X3872);
    kX->set_meson_JP(AXIALVECTOR);

    amplitude exc_gamma = new_amplitude<covariant::photon_exchange>(kX, "Exclusive");
    exc_gamma->set_parameters(parsGamma);

    amplitude exc_omega = new_amplitude<covariant::photon_exchange>(kX, M_OMEGA, "Omega Exchange");
    exc_omega->set_parameters(parsOmega);

    amplitude exc_rho   = new_amplitude<covariant::photon_exchange>(kX, M_RHO, "Rho Exchange");
    exc_rho->set_parameters(parsRho);

    // For the neutron target we need to flip the sign of the coupling for the rho
    amplitude exc_rho_m = new_amplitude<covariant::photon_exchange>(kX, M_RHO, "#minus Rho Exhange");
    exc_rho->set_parameters({-gRho,   etaRho,   lamRho  });

    amplitude exc_mesons_p = exc_omega + exc_rho;
    exc_mesons_p->set_id("Exclusive");
    amplitude exc_mesons_n = exc_omega + exc_rho_m;
    exc_mesons_n->set_id("Exclusive");

    // --------------------------------------------------------------------------
    // Aux functions to help plotting easier

    // Bounds to plot
    std::array<double,2> NT = {Wth, 7};
    
    // Convert cross sections for primakoff into fb
    auto inc_gamma_fb   = [&]  (double W)
    {
        double sig = inc_gamma->integrated_xsection(W*W) * 1E6;
        print(W, sig);
        return sig;
    };
    auto exc_gamma_fb   = [&]  (double W)
    {
        double sig = exc_gamma->integrated_xsection(W*W) * 1E6;
        return sig;
    };

    // For meson exchanges we include both rho and omega contributions
    // However, this depends on isospin of the target
    // For protons they sum, for neutrons they subtract
    auto inc_mesons = [&]  (double W)
    {
        double sig = inc_omega->integrated_xsection(W*W) + inc_rho->integrated_xsection(W*W);
        print(W, sig);
        return sig;
    };    

    // --------------------------------------------------------------------------
    // Plot results

    plotter plotter;

    //---------------------------------------------------
    // Plot of pure primakoff in fb
    plot p1 = plotter.new_plot();

    p1.set_curve_points(30);
    p1.set_logscale(false, true);
    p1.set_ranges({4.6, 7}, {1E-3, 4E1});
    p1.set_legend(0.27, 0.72);
    p1.add_header("Primakoff Effect");
    p1.set_labels( "#it{W}_{#gamma#it{N}}  [GeV]", "#sigma(#gamma#it{N} #rightarrow #it{X}#it{N})  [fb]");

    p1.add_curve( NT, inc_gamma_fb, "Semi-inclusive");
    inc_gamma->set_option(n);
    p1.add_dashed(NT, inc_gamma_fb);

    p1.add_curve( NT, exc_gamma_fb, "Exclusive");
    exc_gamma->set_option(n);
    p1.add_dashed(NT, exc_gamma_fb);

    //---------------------------------------------------
    // Plot of vector meson exchange in nb
    plot p2 = plotter.new_plot();

    p2.set_curve_points(30);
    p2.set_logscale(false, true);
    p2.set_ranges({4.6, 7}, {10E-2, 2E3});
    p2.set_legend(0.27, 0.72);
    p2.add_header("#rho / #omega Exchange");
    p2.set_labels( "#it{W}_{#gamma#it{N}}  [GeV]", "#sigma(#gamma#it{N} #rightarrow #it{X}#it{N})  [nb]");
    p2.color_offset(2);

    p2.add_curve( NT, inc_mesons, "Semi-inclusive");
    inc_rho->set_option(n); inc_omega->set_option(n);
    p2.add_dashed(NT, inc_mesons);

    p2.add_curve( NT, sigma_w, exc_mesons_p);
    exc_mesons_n->set_option(n);
    p2.add_dashed(NT, sigma_w, exc_mesons_n);

    plotter.combine({2,1}, {p1, p2}, "inclusive_X.pdf");
};