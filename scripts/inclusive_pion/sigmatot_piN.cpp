// Draws the total pion-nucleon cross-section using JPAC amplitudes & PDG parameterization
// Reproduces Fig. 3 of [1]
//
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

#include "inclusive_pion/PDG_piN.hpp"
#include "inclusive_pion/JPAC_piN.hpp"
#include "plotter.hpp"

void sigmatot_piN()
{
    using namespace jpacPhoto;

    JPAC_piN sigma;
    PDG_piN  sigma_pdg;

    plotter plotter;
    plot p = plotter.new_plot();

    double Wth = M_PION + M_PROTON + 1E-3;
    p.add_curve({ Wth, 2.5}, [&](double w){ return sigma(+1, w*w, M2_PION); }, "#pi^{#plus} #it{p}");
    p.add_dashed({Wth, 2.5}, [&](double w){ return sigma_pdg(+1, w*w); });
    p.add_curve( {Wth, 2.5}, [&](double w){ return sigma(-1, w*w, M2_PION); }, "#pi^{#minus} #it{p}");
    p.add_dashed({Wth, 2.5}, [&](double w){ return sigma_pdg(-1, w*w); });

    p.set_ranges({1, 2.5}, {6, 400});
    p.set_logscale(false, true);
    p.set_labels("#it{W}_{#gamma#it{p}}  [GeV]", "#sigma_{tot}^{#pi#it{p}} [mb]");
    p.set_legend(0.6, 0.7);
    p.save("sigmatot.pdf");
};