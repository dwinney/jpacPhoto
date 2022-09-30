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
#include "DLambdac_box.hpp"
#include "DstarLambdac_box.hpp"
#include "pomeron_exchange.hpp"
#include "amplitude_sum.hpp"

#include "amplitude_fitter.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

void box_fit()
{

    // ---------------------------------------------------------------------------
    // Preliminaries
    // ---------------------------------------------------------------------------

    reaction_kinematics kJpsi(M_JPSI, M_PROTON);
    kJpsi.set_meson_JP( {1, -1} );

    DLambdac_box      dlam_box( &kJpsi, "#bar{D}#Lambda_{c} box");
    DstarLambdac_box  dslam_box(&kJpsi, "#bar{D}*#Lambda_{c} box");

    // Generate the interpolation grids
    // dlam_box.set_verbose(1); dslam_box.set_verbose(1);
    // dlam_box.generate_grid( 6., {0., 1.5}, {200, 50}, "./grids/");
    // dslam_box.generate_grid(6., {0., 1.5}, {200, 50}, "./grids/");

    // Else you can just import one
    dlam_box.import_grid( "./grids/");
    dslam_box.import_grid("./grids/");

    // Combined the two boxes together
    amplitude_sum boxes(&kJpsi, {&dlam_box, &dslam_box}, "Boxes only");

    // Pomeron
    pomeron_exchange pomeron(&kJpsi, "Pomeron exchange");


    // Fit 1 (destructive)
    dlam_box.set_params({1.5,  0.8409382758});
    dslam_box.set_params({1.5, 0.8876040639});
    pomeron.set_params({ 0.3760808165,  0.2133018688, 1.1, 0.347252063});

    // Complete sum of all 
    amplitude_sum full_sum(&kJpsi, {&dlam_box, &dslam_box, &pomeron}, "Combined Fit");

    // ---------------------------------------------------------------------------
    // Fit to 2019 GlueX data
    // ---------------------------------------------------------------------------

    amplitude_fitter fitter(&full_sum);

    // Raw integrated data

    // GlueX 2019
    std::vector<double> egam      = {8.38 , 8.74,  9.1  , 9.46 , 9.82 , 10.18, 10.54, 10.9 , 11.26, 11.62};
    std::vector<double> err_egam  = {0.18 , 0.18,  0.18 , 0.18 , 0.18 , 0.18 , 0.18 , 0.18 , 0.18 , 0.18 };
    std::vector<double> sigma     = {0.116, 0.343, 0.313, 0.835, 0.868, 0.949, 1.383, 1.274, 2.158, 3.245};
    std::vector<double> err_sigma = {0.033, 0.105, 0.137, 0.278, 0.224, 0.213, 0.430, 0.276, 0.780, 1.00 };

    // Raw differential data
    double egam_avg = 10.7;
    std::vector<double> mtp      = {0.075, 0.225, 0.375, 0.525, 0.675, 0.825, 0.975};
    std::vector<double> err_mtp  = {0.075, 0.075, 0.075, 0.075, 0.075, 0.075, 0.075};
    std::vector<double> dsdt     = {1.642, 1.249, 1.088, 0.627, 0.599, 0.47 , 0.4  };
    std::vector<double> err_dsdt = {0.338, 0.265, 0.248, 0.182, 0.163, 0.145, 0.134};


    fitter.add_integrated_data(           egam, sigma, err_sigma, "GlueX 2019 (integrated)");
    fitter.add_differential_data(egam_avg, mtp,  dsdt,  err_dsdt, "GlueX 2019 (differential)", true); // last bool means we use t' not t
    
    fitter.use_beam_energy(); // Globally use Egamma instead of s for energy variables
    fitter.set_parameter_labels({"qmax D", "eta D", "qmax D*", "eta D*", "A", "b", "alpha(0)", "alpha'"});
    
    // Fix the dispersion cut-offs since they dont affect much
    fitter.fix_parameter("qmax D");
    fitter.fix_parameter("qmax D*");
    
    // Fix the pomeron trajectory to be that of the super high-energy data
    fitter.fix_parameter("alpha(0)");
    fitter.fix_parameter("alpha'");

    fitter.set_parameter_limits("eta D",    {0.0, 1.5}); // Limits here based on the size of our interpolation
    fitter.set_parameter_limits("eta D*",   {0.0, 1.5});
    fitter.set_parameter_limits("b", {0., 5.});          // b-slope should be positive 

    // // Start fit with out initial guess
    fitter.do_fit({1.5, 1., 1.5, 1., +0.5, 0.4, 1.151, 0.112}); 

    // ---------------------------------------------------------------------------
    // Plot the results
    // ---------------------------------------------------------------------------

    std::vector<amplitude*> amps;
    amps.push_back(&full_sum);
    amps.push_back(&dlam_box);
    amps.push_back(&dslam_box);
    amps.push_back(&pomeron);

    jpacGraph1D* plotter = new jpacGraph1D();

    // Options
    int N = 100;
    bool PRINT = false;
    double xmin, xmax, ymin, ymax;
    
    // Plot the integrated differential cross section in data region
    xmin = 8.; xmax = 11.8;
    ymin = 2.E-2; ymax = 20.;
    for (int i = 0; i < amps.size(); i++)
    {
        auto F = [&](double x)
        {
            double w = W_cm(x);
            return amps[i]->integrated_xsection(w*w);
        };

        plotter->AddEntry(N, F, {xmin, xmax}, amps[i]->get_id(), PRINT);
    }
    plotter->AddDataPoints(egam, sigma, err_egam, err_sigma, "GlueX 2019");

    plotter->SetXaxis("#it{E}_{#gamma}  [GeV]", xmin, xmax);
    plotter->SetYaxis("#sigma(#gamma #it{p} #rightarrow #it{J}/#psi #it{p})  [nb]", ymin, ymax);
    plotter->SetYlogscale(1);
    plotter->SetLegend(0.2, 0.71);
    plotter->Plot("integrated_fit.pdf");

    // Clear and plot the differential data also
    plotter->ClearData();

    xmin = 0.; xmax = 1.1;
    ymin = 5.E-3; ymax = 2.E1;
    for (int i = 0; i < amps.size(); i++)
    {
        auto G = [&](double mtp)
        {  
            double savg = pow(W_cm(egam_avg), 2.);
            double t = -mtp + kJpsi.t_man(savg, 0.);
            return amps[i]->differential_xsection(savg, t);
        };

        plotter->AddEntry(N, G, {xmin, xmax}, amps[i]->get_id(), PRINT);
    };
    plotter->AddDataPoints(mtp, dsdt, err_mtp, err_dsdt, "GlueX 2019");

    plotter->SetXaxis("#minus#it{t'}  [GeV^{2}]", xmin, xmax);
    plotter->SetYaxis("#it{d}#sigma/#it{dt} (#gamma #it{p} #rightarrow #it{J}/#psi #it{p})  [nb / GeV^{-2}]", ymin, ymax);
    plotter->SetYlogscale(1);
    plotter->SetLegend(0.5, 0.72);
    plotter->Plot("differential_fit.pdf");

    delete plotter;

    return 0;
};