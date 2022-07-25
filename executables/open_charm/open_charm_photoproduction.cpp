// ---------------------------------------------------------------------------
// Photoproduction cross-sections of the Lambda_c Dbar/D* final state
// by open charm exchanges 
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] arXiv:2009.08345v1
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/dirac_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // Form factor parameter
    double eta = 1.;
    double lambdaQCD = 0.25;
    bool du_result = false;

    // ---------------------------------------------------------------------------
    // D phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for Dbar LambdaC in final state
    reaction_kinematics kD (M_D, M_LAMBDAC);
    kD.set_meson_JP(0, -1);

    vector_exchange d_dstarEx (&kD, M_DSTAR, "D^{*} exchange");
    d_dstarEx.set_params({0.134, -13.2, 0.});
    d_dstarEx.set_formfactor(2, M_DSTAR + eta * lambdaQCD);
    d_dstarEx.set_debug(1);

    dirac_exchange d_lamcEx (&kD, M_LAMBDAC, "#Lambda_{c} exchange");
    d_lamcEx.set_params({sqrt(4.* PI * ALPHA), -4.3});
    d_lamcEx.set_formfactor(2, M_LAMBDAC + eta * lambdaQCD);

    amplitude_sum d_sum (&kD,  {&d_dstarEx, &d_lamcEx}, "Sum");

    // ---------------------------------------------------------------------------
    // D* phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for D* LambdaC in final state
    reaction_kinematics kDstar (M_DSTAR, M_LAMBDAC);
    kDstar.set_meson_JP(1, -1);

    pseudoscalar_exchange dstar_dEx (&kDstar, M_D, "D exchange");
    dstar_dEx.set_params({0.134, -4.3});
    dstar_dEx.set_formfactor(2, M_D + eta * 0.250);
    dstar_dEx.set_debug(1);
    
    vector_exchange dstar_dstarEx (&kDstar, M_DSTAR, "D^{*} exchange");
    dstar_dstarEx.set_params({0.641, -13.2, 0.});
    dstar_dstarEx.set_formfactor(2, M_DSTAR + eta * 0.250);
    dstar_dstarEx.set_debug(1);

    dirac_exchange dstar_lamcEx (&kDstar, M_LAMBDAC, "#Lambda_{c} exchange");
    dstar_lamcEx.set_params({sqrt(4.* PI * ALPHA), -13.2});
    dstar_lamcEx.set_formfactor(2, M_LAMBDAC + eta * 0.250);
    
    amplitude_sum dstar_sum (&kDstar, {&dstar_dEx, &dstar_dstarEx, &dstar_lamcEx}, "D^{*} production");

    if (du_result)
    {
        dstar_dstarEx.set_debug(2);
        d_dstarEx.set_debug(2);
    }

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;

    // amps.push_back(d_dstarEx);
    // amps.push_back(d_lamcEx);
    // amps.push_back(dstar_dEx);
    // amps.push_back(dstar_dstarEx);
    // amps.push_back(dstar_lamcEx);
    amps.push_back(&d_sum);
    amps.push_back(&dstar_sum);

    int N = 50;
    double PRINT = true;

    double xmin = 8.3;
    double xmax = 10.5;

    double ymin = 0.;
    double ymax = 400.;

    std::string filename  = "open_charm.pdf";
    std::string ylabel    = "#sigma(#gamma p #rightarrow D^{(*)} #Lambda_{c}^{+})   [nb]";
    std::string xlabel    = "E_{#gamma}  [GeV]";

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    jpacGraph1D * plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    for (int n = 0; n < amps.size(); n++)
    {
        auto F = [&](double e)
        {
            double w = W_cm(e);
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