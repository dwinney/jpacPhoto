// ---------------------------------------------------------------------------
// 
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] arXiv:2009.08345v1
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "pomeron_exchange.hpp"
#include "amplitude_fitter.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

void fit_test()
{

    // ---------------------------------------------------------------------------
     // ---------------------------------------------------------------------------
    // Preliminaries
    // ---------------------------------------------------------------------------

    // Low - energy trajectory and couplings
    linear_trajectory alpha (1, 0.94, 0.36, "LE");

    reaction_kinematics kjpsi (M_JPSI);
    kjpsi.set_meson_JP(1, -1);

    pomeron_exchange pomeron(&kjpsi, &alpha, false, "Pomeron exchange");

    // ---------------------------------------------------------------------------
    // Fit to 2019 GlueX data
    // ---------------------------------------------------------------------------

    amplitude_fitter fitter(&pomeron);

    // Raw integrated
    std::vector<double> egam      = {8.38,  8.74,  9.1, 9.46, 9.82, 10.18, 10.54, 10.9, 11.26, 11.62};
    std::vector<double> err_egam  = {0.18,  0.18, 0.18, 0.18, 0.18,  0.18,  0.18, 0.18,  0.18,  0.18};
    std::vector<double> sigma     = {0.116, 0.343, 0.313, 0.835, 0.868, 0.949, 1.383, 1.274, 2.158, 3.245};
    std::vector<double> err_sigma = {0.033, 0.105, 0.137, 0.278, 0.224, 0.213, 0.430, 0.276, 0.780, 1.00};

    // Convert Egamma to s
    std::vector<double> s, err_s;
    for (int i = 0; i < egam.size(); i++)
    {
        double W = W_cm(egam[i]);
        s.push_back(W*W);

        double eW = W_cm(err_egam[i]);
        err_s.push_back(eW*eW);
    };

    // Raw differential
    double egam_avg = 10.7;
    std::vector<double> mtp    = {(0.0+0.15)/2., (0.15+0.3)/2., (0.3+0.45)/2., (0.45+0.6)/2., (0.6+0.75)/2., (0.75+0.9)/2., (0.9+1.05)/2.};
    std::vector<double> dsigma = {1.642, 1.249, 1.088, 0.627, 0.599, 0.47, 0.4};
    std::vector<double> derror = {0.338, 0.265, 0.248, 0.182, 0.163, 0.145, 0.134};

    // Convert -t' to t
    std::vector<double> t;
    double savg = pow(W_cm(egam_avg), 2.);
    for (double mtp : mtp)
    {
        double tt = -mtp + kjpsi.t_man(savg, 0.);
        t.push_back(tt);
    };

    std::vector<double> initial_guess = {0.5, 0.4};
    fitter.set_parameter_labels({"A", "b"});

    fitter.add_integrated_data(        s,  sigma, err_sigma, "GlueX 2019 (integrated)");
    fitter.add_differential_data(savg, t, dsigma,    derror, "GlueX 2019 (differential)");
    fitter.do_fit(initial_guess);

    // ---------------------------------------------------------------------------
    // Plot the results
    // ---------------------------------------------------------------------------

    // Options
    int N = 100;
    double  xmin = 8.;
    double  xmax = 12.;

    double  ymin = 7.E-2;
    double  ymax = 20.;

    std::string filename = "fit_results.pdf";
    std::string ylabel  = "#sigma(#gamma #it{p} #rightarrow #it{J}/#psi #it{p})  [nb]";
    bool PRINT = false;


    auto F = [&](double x)
    {
        double w = W_cm(x);
        return pomeron.integrated_xsection(w*w);
    };

    jpacGraph1D* plotter = new jpacGraph1D();
    plotter->AddEntry(N, F, {xmin, xmax}, pomeron.get_id(), PRINT);
    plotter->AddDataPoints(egam, sigma, err_egam, err_sigma, "GlueX 2019");

    plotter->SetXaxis("#it{E}_{#gamma}  [GeV]", xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(1);
    plotter->SetLegend(0.3, 0.67);
    plotter->SetLegendOffset(0.3, 0.1);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

    return 0;
};