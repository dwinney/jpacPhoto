// Plots the total inclusive cross-sections for the Z mesons 
// Reproduces Fig. 8 of [1]
//
// OUTPUT: Z_totals.pdf
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] 	arXiv:2209.05882 [hep-ph]
// ------------------------------------------------------------------------------

#include "constants.hpp"
#include "kinematics.hpp"
#include "plotter.hpp"
#include "semi_inclusive.hpp"
#include "semi_inclusive/pion_exchange.hpp"
#include "semi_inclusive/phase_space.hpp"
#include "analytic/pseudoscalar_exchange.hpp"

void Z_totals()
{
    using namespace jpacPhoto;

    plotter plotter;

    // ---------------------------------------------------------------------------
    // Couplings and constants
    // ---------------------------------------------------------------------------

    // Bottom vertex coupling (pi - nucleon - nucleon)
    double g_piNN = sqrt(2) * sqrt(4*PI*13.81); 

    // Cutoff for exponential form factor
    double lambda_pi = .900;  // MeV 

    // Zc(3900) couplings 
    double gc_jpsi  = 1.91; // psi coupling before VMD scaling
    double gc_gamma = E * F_JPSI * gc_jpsi / M_JPSI;

    // Zb(10610) couplings
    double gb_upsilon1S = 0.49, gb_upsilon2S = 3.30, gb_upsilon3S = 9.22;
    double gb_gamma = E * (F_UPSILON1S * gb_upsilon1S / M_UPSILON1S 
                         + F_UPSILON2S * gb_upsilon2S / M_UPSILON2S
                         + F_UPSILON3S * gb_upsilon3S / M_UPSILON3S);  

    // Zb'(10650) couplings
    double gbp_upsilon1S = 0.21, gbp_upsilon2S = 1.47, gbp_upsilon3S = 4.8;
    double gbp_gamma = E * (F_UPSILON1S * gbp_upsilon1S / M_UPSILON1S 
                         +  F_UPSILON2S * gbp_upsilon2S / M_UPSILON2S
                         +  F_UPSILON3S * gbp_upsilon3S / M_UPSILON3S);  

    // ---------------------------------------------------------------------------
    // Exclusive reactions 
    // ---------------------------------------------------------------------------

    kinematics kZc  = new_kinematics( M_ZC3900);
    kinematics kZb  = new_kinematics(M_ZB10610);
    kinematics kZbp = new_kinematics(M_ZB10650);

    kZc->set_meson_JP( AXIALVECTOR);
    kZb->set_meson_JP( AXIALVECTOR);
    kZbp->set_meson_JP(AXIALVECTOR);

    // Zc(3900)
    amplitude   Zce = new_amplitude<analytic::pseudoscalar_exchange>(kZc,   M_PION, "#it{Z}_{c}(3900)");
    Zce->set_parameters({gc_gamma, g_piNN, lambda_pi});

    // Zb(10610)
    amplitude   Zbe = new_amplitude<analytic::pseudoscalar_exchange>(kZb,   M_PION, "#it{Z}_{b}(10610)");
    Zbe->set_parameters({gb_gamma, g_piNN, lambda_pi});

    // Zb'(10650)
    amplitude   Zbpe = new_amplitude<analytic::pseudoscalar_exchange>(kZbp, M_PION, "#it{Z}_{b}'(10650)");
    Zbpe->set_parameters({gbp_gamma, g_piNN, lambda_pi});

    // ---------------------------------------------------------------------------
    // Z minus production
    // ---------------------------------------------------------------------------

    // Zc(3900)-
    semi_inclusive Zcm = new_semi_inclusive<inclusive::pion_exchange>(kZc, -1, "#it{Z}_{c}(3900)^{#minus}");
    Zcm->set_parameters(gc_gamma);

    // Zb(10610)-
    semi_inclusive Zbm = new_semi_inclusive<inclusive::pion_exchange>(kZb, -1, "#it{Z}_{b}(10610)^{#minus}");
    Zbm->set_parameters(gb_gamma);

    // Zb'(10650)-
    semi_inclusive Zbpm = new_semi_inclusive<inclusive::pion_exchange>(kZbp, -1, "#it{Z}_{b}(10650)^{#minus}");
    Zbpm->set_parameters(gbp_gamma);

    // ---------------------------------------------------------------------------
    // Z plus production
    // ---------------------------------------------------------------------------

    // Zc(3900)+
    semi_inclusive Zcp = new_semi_inclusive<inclusive::pion_exchange>(kZc, +1, "#it{Z}_{c}(3900)^{#plus}");
    Zcp->set_parameters(gc_gamma);
    Zcp += Zce;

    // Zb(10610)+
    semi_inclusive Zbp = new_semi_inclusive<inclusive::pion_exchange>(kZb, +1, "#it{Z}_{b}(10610)^{#plus}");
    Zbp->set_parameters(gb_gamma);
    Zbp += Zbe;

    // Zb'(10650)+
    semi_inclusive Zbpp = new_semi_inclusive<inclusive::pion_exchange>(kZbp, +1, "#it{Z}_{b}(10650)^{#plus}");
    Zbpp->set_parameters(gbp_gamma);
    Zbpp += Zbpe;

    // ---------------------------------------------------------------------------
    // Make plot
    // ---------------------------------------------------------------------------

    int N = 50;
    std::array<double,2> bounds = {4.8, 20};

    // Z minus plot
    auto sig = [](semi_inclusive to_plot)
    {
        return [to_plot](double w){ return to_plot->integrated_xsection(w*w); };
    };
    auto sige = [](amplitude to_plot)
    {
        return [to_plot](double w){ return to_plot->integrated_xsection(w*w); };
    };

    plot m = plotter.new_plot();
    m.set_curve_points(N);
    m.set_legend(0.6, 0.7);
    m.set_logscale(false, true);
    m.set_ranges({4.6,20}, {2E-1, 1E2});
    m.set_labels("#it{W}_{#gammap}  [GeV]", "#sigma [nb]");

    m.add_curve( bounds, sig(Zcm), Zcm->id());
    Zcm->set_option( inclusive::pion_exchange::kPwave );
    m.add_dashed(bounds, sig(Zcm));
    m.add_curve( bounds, sig(Zbm), Zbm->id());
    Zbm->set_option( inclusive::pion_exchange::kPwave );
    m.add_dashed(bounds, sig(Zbm));
    m.add_curve( bounds, sig(Zbpm), Zbpm->id());
    Zbpm->set_option( inclusive::pion_exchange::kPwave );
    m.add_dashed(bounds, sig(Zbpm));

    plot p = plotter.new_plot();
    p.set_curve_points(N);
    p.set_legend(0.6, 0.7);
    p.set_logscale(false, true);
    p.set_ranges({4.6,20}, {2E-1, 1E2});
    p.set_labels("#it{W}_{#gammap}  [GeV]", "#sigma [nb]");

    p.add_curve( bounds, sig(Zcp),  Zcp->id());
    p.add_dashed(bounds, sige(Zce));
    p.add_curve( bounds, sig(Zbp),  Zbp->id());
    p.add_dashed(bounds, sige(Zbe));
    p.add_curve( bounds, sig(Zbpp), Zbpp->id());
    p.add_dashed(bounds, sige(Zbpe));

    // Combine plots together
    plotter.combine({2,1}, {m, p}, "Z_totals.pdf");
};