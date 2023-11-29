// Plots F_1 and F_2 at fixed Q2 as a function of W using the low-energy
// empirical formulae of [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] - https://arxiv.org/abs/0712.3731
// ------------------------------------------------------------------------------

#include "sigma_tots/CB_F.hpp"
#include "plotter.hpp"

void CB_Fs()
{
    using namespace jpacPhoto;

    int p = CB_F::kProton;
    int n = CB_F::kNeutron;

    auto pF1 = new_inclusive_function<CB_F>(1, p); 
    auto pF2 = new_inclusive_function<CB_F>(2, p); 

    auto nF1 = new_inclusive_function<CB_F>(1, n); 
    auto nF2 = new_inclusive_function<CB_F>(2, n); 

    double Q2;     
    int iso;

    // Lambdas for all the different curves we want to plot

    auto F1 = [&](double W)
    {
        return (iso == p) ? pF1->evaluate(W*W, -Q2)
                          : nF1->evaluate(W*W, -Q2);
    };

    auto F2 = [&](double W)
    {
        return (iso == p) ? pF2->evaluate(W*W, -Q2)
                          : nF2->evaluate(W*W, -Q2);
    };

    plotter plotter;

    // F1 plot
    plot p1 = plotter.new_plot();
    p1.set_curve_points(1000);
    p1.set_ranges({1, 3}, {0, 4});
    p1.set_labels("#it{W} [GeV]", "#it{F}_{1}(#it{W}, #it{Q}^{2})");
    p1.set_legend(0.25,0.7);

    Q2 = 0.1;
    iso = p;
    p1.add_curve( {M_PROTON+M_PION, 3}, F1, "#it{Q}^{2} = 0.1 GeV^{2}");
    iso = n;
    p1.add_dashed({M_PROTON+M_PION, 3}, F1);
    
    Q2 = 0.5;
    p1.add_curve( {M_PROTON+M_PION, 3}, F1, "#it{Q}^{2} = 0.5 GeV^{2}");
    iso = n;
    p1.add_dashed({M_PROTON+M_PION, 3}, F1);

    Q2 = 1.0;
    p1.add_curve( {M_PROTON+M_PION, 3}, F1, "#it{Q}^{2} = 1.0 GeV^{2}");
    iso = n;
    p1.add_dashed({M_PROTON+M_PION, 3}, F1);
    
    plot p2 = plotter.new_plot();
    p2.set_curve_points(1000);
    p2.set_ranges({1, 3}, {0, 0.5});
    p2.set_labels("#it{W} [GeV]", "#it{F}_{2}(#it{W}, #it{Q}^{2})");
    p2.set_legend(0.65,0.35);

    Q2 = 0.1;
    iso = p;
    p2.add_curve( {M_PROTON+M_PION, 3}, F2, "#it{Q}^{2} = 0.1 GeV^{2}");
    iso = n;
    p2.add_dashed({M_PROTON+M_PION, 3}, F2);
    
    Q2 = 0.5;
    p2.add_curve( {M_PROTON+M_PION, 3}, F2, "#it{Q}^{2} = 0.5 GeV^{2}");
    iso = n;
    p2.add_dashed({M_PROTON+M_PION, 3}, F2);

    Q2 = 1.0;
    p2.add_curve( {M_PROTON+M_PION, 3}, F2, "#it{Q}^{2} = 1.0 GeV^{2}");
    iso = n;
    p2.add_dashed({M_PROTON+M_PION, 3}, F2);

    plotter.combine({2,1}, {p1,p2}, "CB_Fs.pdf");
};