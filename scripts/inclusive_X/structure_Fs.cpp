// Print out Sach's form factors for the proton and neutron 
// based on the parameterization of [1] and [2]. Reproduces fig. 2 in [3].
//
// Output: CB_pdf (unpublished)
//         Fs_compare.pdf (fig. 2)
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
//               Universitat Bonn, HISKP
// Email:        daniel.winney@iu.alumni.edu
//               winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------
// References:
// [1] - https://arxiv.org/abs/hep-ph/0402081
// [2] - https://arxiv.org/abs/1512.09113
// [3] - https://arxiv.org/abs/2404.05326
// ------------------------------------------------------------------------------

#include "inclusive/structure_functions.hpp"
#include "plotter.hpp"

void structure_Fs()
{
    using namespace jpacPhoto;
    using complex = std::complex<double>;

    auto structure = structure_functions(); 

    double Wth = M_PROTON+M_PION + EPS;
    std::array<double,2> range = {M_PROTON+M_PION, 3};
    auto xB = [&](double W, double t){ return -t / (W*W - M2_PROTON - t); };

    // C&B plot near threshold
    plotter plotter;
    plot p1 = plotter.new_plot();
    p1.set_curve_points(1000);
    p1.set_ranges({1, 3.0}, {0, 0.5});
    p1.set_labels("#it{M}_{#it{X}} [GeV]", "#it{F}_{2}(#it{x}_{B}, #it{t})");
    p1.set_legend(0.35,0.75);
    p1.add_curve( {Wth, 3}, [&](double w){ return structure.F2(w*w, -0.1);}, "#it{t} = #minus 0.1 GeV^{2}");
    p1.add_dashed({Wth, 3}, [&](double w){ return 2*xB(w, -0.1)*structure.F1( w*w, -0.1);});
    p1.add_curve( {Wth, 3}, [&](double w){ return structure.F2(w*w, -2.0);}, "#it{t} = #minus 2.0 GeV^{2}");
    p1.add_dashed({Wth, 3}, [&](double w){ return 2*xB(w, -2.0)*structure.F1(w*w, -2.0);});
    p1.add_curve( {Wth, 3}, [&](double w){ return structure.F2(w*w, -10);}, "#it{t} = #minus 10 GeV^{2}");
    p1.add_dashed({Wth, 3}, [&](double w){ return 2*xB(w, -10 )*structure.F1(w*w, -10);});
    p1.save("CB_Fs.pdf");

    // D&L F1 plot
    plot p2 = plotter.new_plot();
    p2.set_curve_points(100);
    p2.set_ranges({1, 30}, {0, 0.8});
    p2.set_labels("#it{M}_{#it{X}} [GeV]", "#it{F}_{2}(#it{x}_{B}, #it{t})");
    p2.set_legend(0.25,0.75);

    structure.set_option(structure_functions::kDL);
    p2.add_curve(  {Wth, 30}, [&](double w){ return structure.F2(w*w, -0.1);}, "#it{t} = #minus 0.1 GeV^{2}");
    p2.add_dashed( {Wth, 30}, [&](double w){ return structure.F2(w*w, -2);});
    p2.add_dashed( {Wth, 30}, [&](double w){ return structure.F2(w*w, -10);});
    p2.save("DL_Fs.pdf");

};