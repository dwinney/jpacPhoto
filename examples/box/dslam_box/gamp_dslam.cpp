// ---------------------------------------------------------------------------
// Total photoproduction cross-sections of the Lambda_c D final state
// Amplitude constructed from considering open-charm exchanges in t and u channels
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
#include "pseudoscalar_exchange.hpp"
#include "vector_exchange.hpp"
#include "dirac_exchange.hpp"
#include "amplitude_sum.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

void gamp_dslam()
{
    // Form factor parameter
    double eta = 1.;
    double lambdaQCD = 0.25;

    double e     = sqrt(4.* PI * ALPHA);
    double gDsDs = 0.641;
    double gDDs  = 0.134;
    double gNDsL = -4.3;
    double gNDL  = -13.2;

    // ---------------------------------------------------------------------------
    // D phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for Dbar LambdaC in final state
    reaction_kinematics kDs (M_DSTAR, M_LAMBDAC);
    kDs.set_meson_JP(1, -1);

    dirac_exchange ds_lamcEx (&kDs, M_LAMBDAC, "#Lambda_{c} exchange");
    ds_lamcEx.set_params({e, gNDsL, 0.});
    ds_lamcEx.set_formfactor(2, M_LAMBDAC + eta * lambdaQCD);
    ds_lamcEx.force_covariant(true);

    vector_exchange ds_dstarEx (&kDs, M_DSTAR, "D* exchange");
    ds_dstarEx.set_params({gDsDs, gNDsL, 0.});
    ds_dstarEx.set_formfactor(2, M_DSTAR + eta * lambdaQCD);
    ds_dstarEx.force_covariant(true);

    pseudoscalar_exchange ds_dEx (&kDs, M_D, "D exchange");
    ds_dEx.set_params({gDDs, gNDL});
    ds_dEx.set_formfactor(2, M_D + eta * lambdaQCD);
    ds_dEx.force_covariant(true);
    
    amplitude_sum ds_sum (&kDs,  {&ds_dEx, &ds_dstarEx, &ds_lamcEx}, "Sum");

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;

    amps.push_back(&ds_dEx);
    amps.push_back(&ds_dstarEx);
    amps.push_back(&ds_lamcEx);
    amps.push_back(&ds_sum);

    int N = 100;
    double PRINT = true;

    double xmin = 4.;
    double xmax = 5.;

    double ymin = 0.;
    double ymax = 300.;

    std::string filename  = "gamp.pdf";
    std::string ylabel    = "#sigma(#gamma #it{p} #rightarrow #bar{#it{D}}* #Lambda_{c}^{+})   [nb]";
    std::string xlabel    = "#it{W}  [GeV]";

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    jpacGraph1D * plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    for (int n = 0; n < amps.size(); n++)
    {
        auto F = [&](double w)
        {
            return amps[n]->integrated_xsection(w*w);
        };

        plotter->AddEntry(N, F, {xmin, xmax}, amps[n]->get_id(), PRINT);
    };

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.2, 0.75);
    plotter->SetLegendOffset(0.5, 0.15);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

    return 0;
};