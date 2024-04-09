// Draws the Chew-Low plot (physical region) for the b1 inclusive kinematics
// Reproduces Fig. 11 of [1]
//
// OUTPUT: chew_low.pdf 
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

#include "semi_inclusive.hpp"
#include "kinematics.hpp"
#include "semi_inclusive/pion_exchange.hpp"
#include "semi_inclusive/phase_space.hpp"
#include "plotter.hpp"

void chew_low()
{
    using namespace jpacPhoto;

    semi_inclusive b1 = new_semi_inclusive<inclusive::phase_space>(new_kinematics(M_B1), M_PROTON + M_PION, "b_{1} Phase space");

    double s = 9., M = M_PROTON + M_PION;
    std::array<double,2> bounds = {-b1->TMINfromM2(s, M*M) , -b1->TMAXfromM2(s, M*M)};

    plotter plotter;
    plot p = plotter.new_plot();
    p.set_ranges({0, 5.5}, {1, 3.3});
    p.add_curve(bounds, [&](double mt){ return b1->M2MAXfromT(s, -mt); });
    p.add_curve(bounds, [&](double mt){ return b1->M2MINfromT(s, -mt); });
    // p.set_legend(false);
    p.set_labels( "#minus#it{t} [GeV^{2}]", "#it{M}^{2}  [GeV^{2}]");
    p.save("chew_low.pdf");
};