// Plots inclusive unpolarized structure functions F_1 and F_2 at a variety of
// Q2 values as a function of Bjorken x.
// The high-energy parameterization of [1] is used.
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] - https://arxiv.org/abs/hep-ph/0402081
// ------------------------------------------------------------------------------

#include "sigma_tots/DL_F.hpp"
#include "plotter.hpp"

void DL_Fs()
{
    using namespace jpacPhoto;
    using complex = std::complex<double>;

    inclusive_function F1 = new_inclusive_function<DL_F>(1);
    inclusive_function F2 = new_inclusive_function<DL_F>(2);

    std::vector<std::array<double,2>> Q2xs = 
    {
        // ZEUS
        {0.045, 6.21E-7},
        {0.065, 1.02E-6},
        {0.085, 1.17E-6},
        {0.11, 1.51E-6},
        {0.15, 2.07E-6},
        {0.2,  3.15E-6},
        {0.25, 4E-6},
        {0.3, 5.52E-6},
        {0.4, 8.83E-6},
        {0.5, 1.67E-5},
        {0.65, 3.59E-5},
        // H1
        {1.5,  0.000041}, 
        {2.0,  0.000065},
        {2.5,  0.000065},
        {3.5,  0.000105},
        {5.0,  0.000165},
        {6.5,  0.000165},
        {8.5,  0.000260},
        {12.0, 0.000410},
        {15.0, 0.000410},
        {20.0, 0.000650},
        {25.0, 0.000650},
        {35.0, 0.00105},
        {45.0, 0.00165} 
    };

    double Q2;
    auto Q4F2 = [&](double x)
    {
        double s = Q2/x + M2_PROTON - Q2; 
        return Q2*Q2*F2->evaluate(s, -Q2);
    };
    auto Q4F1 = [&](double x)
    {
        double s = Q2/x + M2_PROTON - Q2; 
        return Q2*Q2*F1->evaluate(s, -Q2);
    };

    plotter plotter;

    // F2 plot
    plot p2 = plotter.new_plot();
    p2.set_ranges({5E-7, 1E-3}, {1E-4, 5E3});
    p2.set_logscale(true, true);
    p2.set_labels("#it{x}", "#it{Q}^{4} #it{F}_{2}(#it{x}, #it{Q}^{2})  [GeV^{4}]");

    p2.color_offset(11);
    for (auto Q2x : Q2xs)
    {
        p2.color_offset(-1);
        Q2 = Q2x[0];
        p2.add_curve({Q2x[1], 0.001}, Q4F2);
    }

    // F1 plot
    plot p1 = plotter.new_plot();
    p1.set_ranges({5E-7, 1E-3}, {1E-2, 5E6});
    p1.set_logscale(true, true);
    p1.set_labels("#it{x}", "#it{Q}^{4} #it{F}_{1}(#it{x}, #it{Q}^{2})  [GeV^{4}]");

    p1.color_offset(11);
    for (auto Q2x : Q2xs)
    {
        p1.color_offset(-1);
        Q2 = Q2x[0];
        p1.add_curve({Q2x[1], 0.001}, Q4F1);
    }

    plotter.combine({2,1}, {p1,p2}, "DL_Fs.pdf");
};