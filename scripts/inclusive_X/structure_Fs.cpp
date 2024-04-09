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

#include "semi_inclusive/CB_F.hpp"
#include "semi_inclusive/DL_F.hpp"
#include "semi_inclusive/pdf_F.hpp"

#include "plotter.hpp"

void structure_Fs()
{
    using namespace jpacPhoto;
    using complex = std::complex<double>;

    auto F1_CnB = new_inclusive_function<CB_F>(1),    F2_CnB = new_inclusive_function<CB_F>(2); 
    auto F1_DnL = new_inclusive_function<DL_F>(1),    F2_DnL = new_inclusive_function<DL_F>(2);
    auto F1_PDF = new_inclusive_function<pdf_F>(1),   F2_PDF = new_inclusive_function<pdf_F>(2);

    double Wth = M_PROTON+M_PION + EPS;
    std::array<double,2> range = {M_PROTON+M_PION, 3};
    auto xB = [&](double W, double t){ return -t / (W*W - M2_PROTON - t); };

    // C&B plot near threshold
    plotter plotter;
    plot p1 = plotter.new_plot();
    p1.set_curve_points(100);
    p1.set_ranges({1, 3.0}, {0, 0.5});
    p1.set_labels("#it{M}_{#it{X}} [GeV]", "#it{F}_{2}(#it{x}_{B}, #it{t})");
    p1.set_legend(0.35,0.75);
    p1.add_curve( {Wth, 3}, [&](double w){ return F2_CnB->evaluate(w*w, -0.1);}, "#it{t} = #minus 0.1 GeV^{2}");
    p1.add_dashed({Wth, 3}, [&](double w){ return 2*xB(w, -0.1)*F1_CnB->evaluate(w*w, -0.1);});
    p1.add_curve( {Wth, 3}, [&](double w){ return F2_CnB->evaluate(w*w, -2.0);}, "#it{t} = #minus 2.0 GeV^{2}");
    p1.add_dashed({Wth, 3}, [&](double w){ return 2*xB(w, -2.0)*F1_CnB->evaluate(w*w, -2.0);});
    p1.add_curve( {Wth, 3}, [&](double w){ return F2_CnB->evaluate(w*w, -10);}, "#it{t} = #minus 10 GeV^{2}");
    p1.add_dashed({Wth, 3}, [&](double w){ return 2*xB(w, -10 )*F1_CnB->evaluate(w*w, -10);});
    p1.save("CB_Fs.pdf");

    // D&L F1 plot
    plot p2 = plotter.new_plot();
    p2.set_curve_points(100);
    p2.set_ranges({1, 30}, {0, 0.8});
    p2.set_labels("#it{M}_{#it{X}} [GeV]", "#it{F}_{2}(#it{x}_{B}, #it{t})");
    p2.set_legend(0.25,0.75);

    p2.add_curve(  {Wth, 30}, [&](double w){ return F2_DnL->evaluate(w*w, -0.1);}, "#it{t} = #minus 0.1 GeV^{2}");
    p2.add_curve(  {Wth, 30}, [&](double w){ return F2_PDF->evaluate(w*w, -2);},   "#it{t} = #minus 2 GeV^{2}");
    p2.add_dashed( {Wth, 30}, [&](double w){ return F2_DnL->evaluate(w*w, -2);});
    p2.add_curve(  {Wth, 30}, [&](double w){ return F2_PDF->evaluate(w*w, -10);},  "#it{t} = #minus 10 GeV^{2}");
    p2.add_dashed( {Wth, 30}, [&](double w){ return F2_DnL->evaluate(w*w, -10);});

    plot p3 = plotter.new_plot();
    p3.set_labels("#it{M}_{#it{X}} [GeV]", "#it{F}_{2}(#it{x}_{B}, #it{t})");
    p3.add_style_legend({"D&L", "CTEQ-TEA", "B&C"});
    p3.set_style_legend(0.25,0.75);
    p3.set_curve_points(100);
    p3.set_ranges({1, 3}, {0, 0.45});
    p3.add_curve( {1, 3},  [&](double M){ return F2_DnL->evaluate(M*M, -0.1) ; });
    p3.add_dotted( {1, 3},  [&](double M){ return F2_CnB->evaluate(M*M, -0.1) ; });
    p3.add_curve( {1, 3},  [&](double M){ return F2_DnL->evaluate(M*M, -2) ; });
    p3.add_dashed( {1, 3},  [&](double M){ return F2_PDF->evaluate(M*M, -2) ; });
    p3.add_dotted( {1, 3},  [&](double M){ return F2_CnB->evaluate(M*M, -2) ; });
    p3.add_curve( {1, 3},  [&](double M){ return F2_DnL->evaluate(M*M, -10) ; });
    p3.add_dashed( {1, 3},  [&](double M){ return F2_PDF->evaluate(M*M, -10) ; });
    p3.add_dotted( {1, 3},  [&](double M){ return F2_CnB->evaluate(M*M, -10) ; });
    plotter.combine({2,1}, {p3,p2}, "Fs_compare.pdf");
};