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

void combined_fit()
{
    // ---------------------------------------------------------------------------
    // Toggles to make choosing scenarios easier

    // Amplitude selection
    bool POM       = 1;
    bool DLAM      = 1;
    bool DSLAM     = 1;

    // Data selection
    bool GX2019DIF = 1;
    bool GX2022DIF = 1;
    bool GX2019INT = 1;
    bool GX2022INT = 0;

    std::string fit_curve_label = "Sum";

    // Minuit print level and strategy
    int Nprint = 0;
    std::string strategy;
    strategy = "Simplex";
    // strategy = "Combined";

    // Starting guesses for each amplitude if they're being used
    std::vector<double> pom_guess   = {0.5, 0.3, 1.08, 0.25};
    std::vector<double> dlam_guess  = {1.5, 1.};
    std::vector<double> dslam_guess = {1.5, 1.};

    // ---------------------------------------------------------------------------
    // Amplitude setup
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

    // Pomeron
    // We dont predefine a pomeron trajectory so that they get treated as free parameters in the fit
    pomeron_exchange pomeron(&kJpsi, "Pomeron exchange");

    // Push all the amplitudes we want to include in the fit into a vector
    std::vector<amplitude*> amps;
    if (POM)   amps.push_back(&pomeron);
    if (DLAM)  amps.push_back(&dlam_box);
    if (DSLAM) amps.push_back(&dslam_box);

    // This object is what actually goes into the fitter
    amplitude_sum amp_to_fit(&kJpsi, amps, fit_curve_label);

    // ---------------------------------------------------------------------------
    //  GlueX data set selection
    // ---------------------------------------------------------------------------

    // Raw integrated data

    // GlueX 2019
    std::string int_id = "GlueX 2019 (integrated)";
    std::vector<double> egam      = {8.38 , 8.74,  9.1  , 9.46 , 9.82 , 10.18, 10.54, 10.9 , 11.26, 11.62};
    std::vector<double> err_egam  = {0.18 , 0.18,  0.18 , 0.18 , 0.18 , 0.18 , 0.18 , 0.18 , 0.18 , 0.18 };
    std::vector<double> sigma     = {0.116, 0.343, 0.313, 0.835, 0.868, 0.949, 1.383, 1.274, 2.158, 3.245};
    std::vector<double> err_sigma = {0.033, 0.105, 0.137, 0.278, 0.224, 0.213, 0.430, 0.276, 0.780, 1.00 };

    // Raw differential data
    
    // Gluex 2019 
    double egam_avg = 10.7;
    std::string dif_id = "GlueX 2019 (E = 10.7)";
    std::vector<double> mtp      = {0.075, 0.225, 0.375, 0.525, 0.675, 0.825, 0.975};
    std::vector<double> err_mtp  = {0.075, 0.075, 0.075, 0.075, 0.075, 0.075, 0.075};
    std::vector<double> dsdt     = {1.642, 1.249, 1.088, 0.627, 0.599, 0.47 , 0.4  };
    std::vector<double> err_dsdt = {0.338, 0.265, 0.248, 0.182, 0.163, 0.145, 0.134};

    // ---------------------------------------------------------------------------
    //  Fitter setup
    // ---------------------------------------------------------------------------

    amplitude_fitter fitter(&amp_to_fit, strategy);
    fitter.use_beam_energy(); // Globally use Egamma instead of s for energy variables
    fitter.set_print_level(Nprint);

    // -----------------------------------------
    // Choose which data we want to fit against

    // 2019
    if (GX2019INT) fitter.add_integrated_data(egam, sigma, err_sigma, int_id);
    if (GX2019DIF) fitter.add_differential_data(egam_avg,  mtp,   dsdt,   err_dsdt,  dif_id, true ); // last bool means we use t' not t
    
    // // 2022
    // if (GX2022INT) fitter.add_integrated_data(egam2, sigma2, err_sigma2, int_id2);
    // if (GX2022DIF) fitter.add_differential_data(  egams2,  mt2,  dsdt2,  err_dsdt2, dif_id2, false);    // this one uses raw t values

    // -----------------------------------------
    // Set up our parameter labels and options

    std::vector<std::string> pom_labels, dlam_labels, dslam_labels, labels;
    pom_labels   = {"A", "b",  "alpha(0)", "alpha'"};
    dlam_labels  = {"qmax D",  "eta D" };
    dslam_labels = {"qmax D*", "eta D*"};
    std::vector<double> initial_guess;

    if (POM)
    {
        labels.insert(  labels.end(), pom_labels.begin(),  pom_labels.end());
        initial_guess.insert( initial_guess.end(), pom_guess.begin(), pom_guess.end());
    } 
    if (DLAM)
    {
        labels.insert( labels.end(), dlam_labels.begin(),  dlam_labels.end());
        initial_guess.insert( initial_guess.end(), dlam_guess.begin(), dlam_guess.end());
    }
    if (DSLAM)
    { 
        labels.insert(labels.end(), dslam_labels.begin(), dslam_labels.end());
        initial_guess.insert( initial_guess.end(), dslam_guess.begin(), dslam_guess.end());
    };
    fitter.set_parameter_labels(labels);

    // Fix the dispersion cut-offs since they dont affect much
    if (DLAM)  
    {
        fitter.fix_parameter("qmax D");
        fitter.set_parameter_limits("eta D",    { 0.0, 1.5}); // Limits here based on the size of our interpolation
    }
    if (DSLAM)
    {
        fitter.fix_parameter("qmax D*");      
        fitter.set_parameter_limits("eta D*",   { 0.0, 1.5});
    } 
    if (POM)
    {
        // Pomeron parameters should all be positive
        fitter.set_parameter_limits("b",        {  0., 5.}); 
        fitter.set_parameter_limits("alpha(0)", {  0., 5.}); 
        fitter.set_parameter_limits("alpha'",   {  0., 5.}); 
    }  

    // Do fit and add out fit result to amplitudes to plot
    // If all amplitudes are turned off, skip this step
    if (POM + DLAM + DSLAM != 0) 
    {
        fitter.do_fit(initial_guess); 
     
         // If we fit a combination of multiple amplitudes add the fit result to the vector to be plotted
         // if we only fit a single amplitude dont do this to avoid double counting the curve
        if (POM + DLAM + DSLAM != 1) amps.push_back(&amp_to_fit);
    };

    // ---------------------------------------------------------------------------
    // Plot the results
    // ---------------------------------------------------------------------------

    jpacGraph1D* plotter = new jpacGraph1D();

    // Plotting Options
    int N = 100;
    bool PRINT = false;
    double xmin, xmax, ymin, ymax;
    std::ostringstream ss;
    
    // -----------------------------------------
    // Integrated cross-section plots

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
    
    if (GX2022INT) plotter->AddDataPoints(egam2, sigma2, err_egam2, err_sigma2, "GlueX 2022");
    if (GX2019INT) plotter->AddDataPoints(egam,   sigma,  err_egam,  err_sigma, "GlueX 2019");

    plotter->SetXaxis("#it{E}_{#gamma}  [GeV]", xmin, xmax);
    plotter->SetYaxis("#sigma(#gamma #it{p} #rightarrow #it{J}/#psi #it{p})  [nb]", ymin, ymax);
    plotter->SetYlogscale(1);
    plotter->SetLegend(0.2, 0.71);
    plotter->Plot("integrated_results.pdf");

    // Clear and plot the differential data also
    plotter->ClearData();

    if (GX2019DIF)
    {
        xmin = 0.; xmax = 1.1;
        ymin = 5.E-3; ymax = 5.E1;
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

        ss << std::setprecision(3) << "#it{E}_{#gamma} = " << egam_avg << " GeV";
        plotter->SetLegend(0.5, 0.69, ss.str());
        plotter->Plot("differential_2019.pdf");

        // Clear and plot the second differential data set
        plotter->ClearData();
    }
    if (GX2022DIF)
    {
        double savg2 = pow(W_cm(egam_avg2), 2.);
        xmin = -kJpsi.t_man(savg2, PI); xmax = -kJpsi.t_man(savg2, 0.);
        ymin = 1.E-3; ymax = 10.;
        for (int i = 0; i < amps.size(); i++)
        {
            auto G = [&](double mt)
            {  
                return amps[i]->differential_xsection(savg2, -mt);
            };

            plotter->AddEntry(N, G, {xmin, xmax}, amps[i]->get_id(), PRINT);
        };
        plotter->AddDataPoints(mt2, dsdt2, err_mt2, err_dsdt2, "GlueX 2022");
        plotter->SetXaxis("#minus#it{t}  [GeV^{2}]", 0.6, 5.6);
        plotter->SetYaxis("#it{d}#sigma/#it{dt} (#gamma #it{p} #rightarrow #it{J}/#psi #it{p})  [nb / GeV^{-2}]", ymin, ymax);
        plotter->SetYlogscale(1);

        ss.str(""); ss.clear();
        ss << std::setprecision(3) << "#it{E}_{#gamma} = " << egam_avg2 << " GeV";
        plotter->SetLegend(0.5, 0.69, ss.str());
        plotter->Plot("differential_2022.pdf");
    }
    
    delete plotter;
};
