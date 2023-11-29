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

#include "sigma_tots/PDG_piN.hpp"
#include "sigma_tots/JPAC_piN.hpp"
#include "plotter.hpp"

void sigmatot_piN()
{
    using namespace jpacPhoto;

    inclusive_function PDG_pimp  = new_PDG_sigmatot(pimp);
    inclusive_function JPAC_pimp = new_inclusive_function<JPAC_piN>(-1);

    inclusive_function PDG_pipp  = new_PDG_sigmatot(pipp);
    inclusive_function JPAC_pipp = new_inclusive_function<JPAC_piN>(+1);

    inclusive_function to_plot;
    auto sig = [&](double w)
    {
        return to_plot->evaluate(w*w, M2_PION);
    };

    plotter plotter;
    plot p = plotter.new_plot();

    double Wth = M_PION + M_PROTON + 1E-3;
    to_plot = JPAC_pipp;
    p.add_curve({Wth, 2.5}, sig, "#pi^{#plus} #it{p}");
    to_plot = PDG_pipp;
    p.add_dashed({Wth,2.5}, sig);

    to_plot = JPAC_pimp;
    p.add_curve({Wth, 2.5}, sig, "#pi^{#minus} #it{p}");
    to_plot = PDG_pimp;
    p.add_dashed({Wth,2.5}, sig);

    p.set_ranges({1, 2.5}, {6, 400});
    p.set_logscale(false, true);
    p.set_labels("#it{W}_{#gamma#it{p}}  [GeV]", "#sigma_{tot}^{#pi#it{p}} [mb]");
    p.set_legend(0.6, 0.7);
    p.save("sigmatot.pdf");
};