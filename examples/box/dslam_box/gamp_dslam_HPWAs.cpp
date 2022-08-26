// ---------------------------------------------------------------------------
// Photoproduction cross-sections of the Lambda_c D final state
// Constructed by considering individual s-channel partial wave projections
//
// Here we consider up to J = 5/2
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.edu
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "pseudoscalar_exchange.hpp"
#include "vector_exchange.hpp"
#include "dirac_exchange.hpp"
#include "amplitude_sum.hpp"
#include "projected_amplitude.hpp"
#include "jpacGraph1D.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

void gamp_dslam_HPWAs()
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
    // PW projection
    // ---------------------------------------------------------------------------
    
    // Take the sum ampitude and pass it to a projected_amplitude
    helicity_PWA hpwa(&ds_sum, 1, 3);

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    int N = 50;
    double PRINT = true;

    double xmin = sqrt(kDs.sth()) + 0.01;
    double xmax = 5.;

    double ymin = -0.8;
    double ymax =  1.0;

    std::string filename  = "gamp_dslam_hpwas.pdf";
    std::string ylabel    = "#it{b}_{#{}{#it{L}}}^{1/2} (#it{s})";
    std::string xlabel    = "#it{W}  [GeV]";

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    jpacGraph1D * plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    auto reF = [&](double w)
    {
        return hpwa.real_part(w*w);
    };
    auto imF = [&](double w)
    {
        return hpwa.imag_part(w*w);
    };

    hpwa.set_helicities({1, 1, +1, +1});
    plotter->AddEntry(N, reF, {xmin, xmax}, "#{}{ + +, + + }", PRINT);
    plotter->AddDashedEntry(N, imF, {xmin, xmax}, PRINT);

    hpwa.set_helicities({1, 1, 0, +1});
    plotter->AddEntry(N, reF, {xmin, xmax}, "#{}{ + +, 0 + }", PRINT);
    plotter->AddDashedEntry(N, imF, {xmin, xmax}, PRINT);

    hpwa.set_helicities({1, 1, 0, -1});
    plotter->AddEntry(N, reF, {xmin, xmax}, "#{}{ + +, 0 #minus }", PRINT);
    plotter->AddDashedEntry(N, imF, {xmin, xmax}, PRINT);

    hpwa.set_helicities({1, 1, -1, -1});
    plotter->AddEntry(N, reF, {xmin, xmax}, "#{}{ + +, #minus #minus }", PRINT);
    plotter->AddDashedEntry(N, imF, {xmin, xmax}, PRINT);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.2, 0.3);
    plotter->SetLegendOffset(0.5, 0.15);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

    return 0;
};