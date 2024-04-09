// Print out Sach's form factors for the proton and neutron 
// based on the parameterization of [1]. Reproduces fig. 2 in [2].
//
// Output: Gp.pdf 
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
//               Universitat Bonn, HISKP
// Email:        daniel.winney@iu.alumni.edu
//               winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------
// References:
// [1] - https://arxiv.org/abs/1707.09063
// [2] - https://arxiv.org/abs/2404.05326
// ------------------------------------------------------------------------------

#include "plotter.hpp"
#include "covariant/photon_exchange.hpp"

void sachs_FFs()
{
    using namespace jpacPhoto;
    using covariant::photon_exchange;

    double mu = 1.79;

    kinematics kX = new_kinematics(M_X3872);
    kX->set_meson_JP(AXIALVECTOR);

    auto primakoff = new_raw_amplitude<covariant::photon_exchange>(kX);

    auto GD = [&] (double Q2) { return 1/pow(1 + Q2/0.71, 2.); };
    auto GE_over_GD = [&] (double Q2) { return primakoff->G_E(Q2) / GD(Q2); };
    auto GM_over_GD = [&] (double Q2) { return primakoff->G_M(Q2) / GD(Q2) / mu; };

    plotter plotter;
    plot p = plotter.new_plot();
    p.set_curve_points(1000);
    p.set_logscale( true, false);
    p.set_ranges({ 1E-2, 1E2}, {-0.6, 1.2} );
    p.set_legend(0.3, 0.3);
    p.set_labels("#minus #it{t}  [GeV^{2}]", "#it{G} / #it{G}_{#it{D}}");

    p.add_curve({1E-2, 1E2}, GE_over_GD , "#it{G}_{#it{E}}");
    p.add_curve({1E-2, 1E2}, GM_over_GD , "#it{G}_{#it{M}} / #mu");
    p.save("Gp.pdf");
};
