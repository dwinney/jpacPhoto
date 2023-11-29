// Print out Sach's form factors for the proton and neutron 
// based on the parameterization of [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------
// References:
// [1] - https://arxiv.org/abs/1707.09063
// ------------------------------------------------------------------------------

#include "plotter.hpp"
#include "covariant/photon_exchange.hpp"

void sachs_FFs()
{
    using namespace jpacPhoto;
    
    kinematics kX = new_kinematics(M_X3872);
    kX->set_meson_JP(AXIALVECTOR);

    auto primakoff = new_raw_amplitude<covariant::photon_exchange>(kX);

    // Standard dipole form factor
    auto GD = [&] (double Q2)
    {
        return 1/pow(1 + Q2/0.71, 2.);
    };

    auto GE_over_GD = [&] (double Q2)
    {
        return primakoff->G_E(Q2) / GD(Q2);
    };

    auto GM_over_GD = [&] (double Q2)
    {
        return primakoff->G_M(Q2) / GD(Q2);
    };

    plotter plotter;

    plot pE = plotter.new_plot();
    pE.set_curve_points(1000);
    pE.set_logscale( true, false);
    pE.set_ranges({ 1E-2, 1E2}, {-0.6, 1.1} );
    pE.set_labels("#it{Q}^{2}  [GeV^{2}]", "G_{#it{E}} / G_{#it{D}}");

    primakoff->set_option( photon_exchange::kProton );
    pE.add_curve({1E-2, 1E2}, GE_over_GD , "#it{p}");
    primakoff->set_option( photon_exchange::kNeutron );
    pE.add_curve({1E-2, 1E2}, GE_over_GD , "#it{n}");

    plot pM = plotter.new_plot();
    pM.color_offset(2);
    pM.set_curve_points(1000);
    pM.set_logscale( true, false);
    pM.set_ranges({ 1E-2, 1E2}, {-0.2, 1.2} );
    pM.set_legend(0.3, 0.4);
    pM.set_labels("#it{Q}^{2}  [GeV^{2}]", "G_{#it{M}} / #mu_{#it{p}}G_{#it{D}}");

    primakoff->set_option( photon_exchange::kProton );
    pM.add_curve({1E-2, 1E2}, GM_over_GD , "#it{p}");
    primakoff->set_option( photon_exchange::kNeutron );
    pM.add_curve({1E-2, 1E2}, GM_over_GD , "#it{n}");

    plotter.combine({2,1}, {pE, pM}, "GEp.pdf");
};
