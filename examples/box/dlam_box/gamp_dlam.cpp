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

void gamp_dlam()
{
    // Form factor parameter
    double eta = 1.;
    double lambdaQCD = 0.25;

    // ---------------------------------------------------------------------------
    // D phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for Dbar LambdaC in final state
    reaction_kinematics kD (M_D, M_LAMBDAC);
    kD.set_meson_JP(0, -1);

    vector_exchange d_dstarEx (&kD, M_DSTAR, "D* exchange");
    d_dstarEx.set_params({0.134, -4.3, 0.});
    d_dstarEx.set_formfactor(2, M_DSTAR + eta * lambdaQCD);
    d_dstarEx.force_covariant(true);

    dirac_exchange d_lamcEx (&kD, M_LAMBDAC, "#Lambda_{c} exchange");
    d_lamcEx.set_params({sqrt(4.* PI * ALPHA), -13.2, 0.});
    d_lamcEx.set_formfactor(2, M_LAMBDAC + eta * lambdaQCD);
    d_lamcEx.force_covariant(true);

    amplitude_sum d_sum (&kD,  {&d_dstarEx, &d_lamcEx}, "Sum");

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;

    amps.push_back(&d_dstarEx);
    amps.push_back(&d_lamcEx);
    amps.push_back(&d_sum);

    int N = 100;
    double PRINT = true;

    double xmin = 4.;
    double xmax = 5.;

    double ymin = 0.;
    double ymax = 200.;

    std::string filename  = "gamp.pdf";
    std::string ylabel    = "#sigma(#gamma #it{p} #rightarrow #bar{#it{D}} #Lambda_{c}^{+})   [nb]";
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
    plotter->SetLegendOffset(0.5, 0.1);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

    return 0;
};