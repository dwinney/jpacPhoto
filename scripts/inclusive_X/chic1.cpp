// Calculates the integrated cross section for both inclusive and exclusive 
// chi_c1 via photon and vector meson exchanges
//
// OUTPUT: inclusive_chic1.pdf
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
#include "analytic/vector_exchange.hpp"

void chic1()
{
    using namespace jpacPhoto;

    const int p = inclusive::vector_exchange::kProton;
    const int n = inclusive::vector_exchange::kNeutron;

    double Wth = M_CHIC1 + M_PROTON;
    
    //----------------------------------------------------------------------------
    // INPUTS

    // Couplings
    double gRho    = 1.140E-3;
    double gOmega  = 0.190E-3;
    double gPhi    = 0.110E-3;
    double gPsi    = 34.64E-3;

    // VMD factors
    double etaRho   = 16.37;
    double etaOmega = 56.34;
    double etaPhi   = 44.37;
    double etaPsi   = 36.85;

    // Form factor cutoffs
    double lamRho   = 1.4;
    double lamOmega = 1.2;

    std::vector<double> parsRho, parsOmega, parsPhi, parsPsi;
    parsRho   = {gRho,   etaRho,   lamRho  };
    parsOmega = {gOmega, etaOmega, lamOmega};
    parsPhi   = {gPhi,   etaPhi,   lamOmega};
    parsPsi   = {gPsi,   etaPsi,   0.};

    //----------------------------------------------------------------------------
    // Set up inclusive amplitudes

    inclusive_process inc_omega = new_inclusive_process<inclusive::vector_exchange>(M_CHIC1, M_OMEGA, "#omega");
    inc_omega->set_parameters(parsOmega);

    inclusive_process inc_rho   = new_inclusive_process<inclusive::vector_exchange>(M_CHIC1, M_RHO, "#rho");
    inc_rho->set_parameters(parsRho);

    inclusive_process inc_phi   = new_inclusive_process<inclusive::vector_exchange>(M_CHIC1, M_PHI, "#phi");
    inc_phi->set_parameters(parsPhi);

    std::vector<inclusive_process> exchanges = {inc_rho, inc_omega, inc_phi};

    //----------------------------------------------------------------------------
    // Set up exclusives amplitudes

    kinematics kChi = new_kinematics(M_CHIC1);
    kChi->set_meson_JP(AXIALVECTOR);

    amplitude exc_omega = new_amplitude<covariant::photon_exchange>(kChi, M_OMEGA, "#omega");
    exc_omega->set_parameters(parsOmega);

    amplitude exc_rho   = new_amplitude<covariant::photon_exchange>(kChi, M_RHO, "#rho");
    exc_rho->set_parameters(parsRho);

    amplitude exc_phi   = new_amplitude<covariant::photon_exchange>(kChi, M_PHI, "#phi");
    exc_phi->set_parameters(parsPhi);

    amplitude exc_psi   = new_amplitude<covariant::photon_exchange>(kChi, M_JPSI, "J/#psi");
    exc_psi->set_parameters(parsPsi);

    // For the neutron target we need to flip the sign of the coupling for the rho
    amplitude exc_rho_m = new_amplitude<covariant::photon_exchange>(kChi, M_RHO, "#minus #rho Exhange");
    exc_rho->set_parameters({-gRho,   etaRho,  lamRho  });

    amplitude exc_mesons_p = exc_rho + exc_omega  + exc_phi + exc_psi;
    exc_mesons_p->set_id("Exclusive");
    amplitude exc_mesons_n = exc_rho_m + exc_omega + exc_phi + exc_psi;
    exc_mesons_n->set_id("Exclusive");

    //----------------------------------------------------------------------------
    // Compare with explicitly 

    // Nucleon couplings 
    double gV_omega = 16.,  gT_omega = 0.;
    double gV_rho = 2.4,    gT_rho = 14.6;
    double gV_phi = -6.2,   gT_phi = 2.1;
    double gV_psi = 1.6E-3, gT_psi = 0.;
    
    // Photon couplings
    double gChi_omega   = 5.2E-4;
    double gChi_rho     = 9.2E-4;
    double gChi_phi     = 4.2E-4;
    double gChi_psi     = 1.;
    double gX_omega     = 8.2E-3;
    double gX_rho       = 3.6E-3;

    // chi_c1
    amplitude ChiC1_omegaL = new_amplitude<analytic::vector_exchange>(kChi, M_OMEGA, "#omega exchange");
    ChiC1_omegaL->set_parameters({gChi_omega, gV_omega, gT_omega, lamOmega});

    amplitude ChiC1_rhoL = new_amplitude<analytic::vector_exchange>(kChi, M_RHO, "#rho exchange");
    ChiC1_rhoL->set_parameters({gChi_rho, gV_rho, gT_rho, lamRho});

    amplitude ChiC1_phiL = new_amplitude<analytic::vector_exchange>(kChi, M_PHI, "#phi exchange");
    ChiC1_phiL->set_option(analytic::vector_exchange::kNoFF);
    ChiC1_phiL->set_parameters({gChi_phi, gV_phi, gT_phi});

    amplitude ChiC1_psiL = new_amplitude<analytic::vector_exchange>(kChi, M_JPSI, "#it{J}/#psi exchange");
    ChiC1_psiL->set_option(analytic::vector_exchange::kNoFF);

    ChiC1_psiL->set_parameters({gChi_psi, gV_psi, gT_psi});
    
    amplitude ChiC1_L = ChiC1_rhoL + ChiC1_omegaL + ChiC1_phiL + ChiC1_psiL;
    ChiC1_L->set_id("Sum");

    // --------------------------------------------------------------------------
    // Aux functions to help plotting easier

    // Bounds to plot
    std::array<double,2> NT = {Wth + EPS, 7};
    
    // For meson exchanges we include both rho and omega contributions
    auto inc_mesons = [&]  (double W)
    {
        double sig = 0;
        for (auto ex : exchanges)
        {
            sig += ex->integrated_xsection(W*W);
        };
        print(W, sig);
        return sig;
    };   

    int ex = 0;
    auto inc_single = [&] (double W)
    {
        return exchanges[ex]->integrated_xsection(W*W);
    };

    // --------------------------------------------------------------------------
    // Plot results

    plotter plotter;

    //---------------------------------------------------
    // 
    plot p1 = plotter.new_plot();
    
    p1.set_curve_points(100);
    p1.set_logscale(false, true);
    p1.set_ranges({4.35, 7}, {1E-4, 300});
    p1.set_legend(0.22, 0.7);
    p1.add_header("Exclusive");
    p1.set_labels( "#it{W}_{#gamma#it{p}}  [GeV]", "#sigma(#gamma#it{p} #rightarrow #chi#it{p})  [nb]");

    std::vector<amplitude> exc_exchanges = extract_subamplitudes(exc_mesons_p);
    std::vector<amplitude> hadron_exc    = extract_subamplitudes(ChiC1_L);
    for (int i = 0; i < exc_exchanges.size(); i++)
    {
        p1.add_curve( NT,  sigma_w, exc_exchanges[i]);
        p1.add_dashed(NT, sigma_w,  hadron_exc[i]);
    };

    // Plot of vector meson exchange in nb
    plot p2 = plotter.new_plot();

    p2.set_curve_points(40);
    p2.set_logscale(false, true);
    p2.set_ranges({4.3, 7}, {1E-3, 4E1});
    p2.set_legend(0.22, 0.70);
    p2.add_header("Semi-inclusive");
    p2.set_labels( "#it{W}_{#gamma#it{p}}  [GeV]", "#sigma(#gamma#it{p} #rightarrow #chi#it{p})  [nb]");

    for (int i = 0; i < exchanges.size(); i++)
    {
        ex = i;
        p2.add_curve( NT, inc_single, exchanges[ex]->id());
    }
    p2.color_offset(1);
    p2.add_curve( NT, inc_mesons, "Sum");

    plotter.combine({2,1}, {p1, p2}, "chic1.pdf");
};