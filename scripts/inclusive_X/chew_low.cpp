// Draws the Chew-Low plot (physical region) for the X(3872) inclusive kinematics
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

#include "inclusive_process.hpp"
#include "inclusive/phase_space.hpp"
#include "plotter.hpp"
#include "print.hpp"

using namespace jpacPhoto;

void chew_low()
{
    using namespace jpacPhoto;
    using complex = std::complex<double>;

    inclusive_process X = new_inclusive_process<inclusive::phase_space>(M_X3872, M_PROTON + M_PION, "X(3872) Phase space");

    double Wth = M_PROTON + M_PION;
    std::array<double,2> bounds;

    plotter plotter;
    plot p = plotter.new_plot();
    p.set_ranges({-0.1, 35}, {0.8, 3.5});

    auto M2max = [&](double mt)
    {
        return sqrt(X->M2MAXfromT(-mt));
    };
    auto M2min = [&](double mt)
    {
        return sqrt(X->M2MINfromT(-mt));
    };
    
    auto add_curve = [&](double W)
    {
        X->set_total_energy(W*W);
        bounds = {-X->TMINfromM2(Wth*Wth) , -X->TMAXfromM2(Wth*Wth)};

        p.add_curve(bounds, M2max, var_def("#sqrt{#it{s}}", W, "GeV"));
        p.color_offset(-1);
        p.add_curve(bounds, M2min);
    };

    add_curve(7);
    add_curve(6);
    add_curve(5);

    p.set_legend(0.73, 0.65);
    p.set_labels( "#it{Q}^{2} [GeV^{2}]", "#it{W}  [GeV]");
    p.save("chew_low.pdf");
};