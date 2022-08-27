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

void psip_dslam_HPWAs()
{
    // Form factor parameter
    double eta = 1.;
    double lambdaQCD = 0.25;

    double gPsiDD   = 7.4;
    double gPsiDDs  = 3.83766;
    double gPsiDsDs = 7.999;
    double gDNL     = -13.2;
    double gDsNL    = -4.3;
    double gPsiLL   = -1.4;

    // ---------------------------------------------------------------------------
    // psi p -> D Lambda amplitudes
    // ---------------------------------------------------------------------------

    // Set up Kinematics for Dbar LambdaC in final state
    reaction_kinematics kDs (M_JPSI, M_PROTON, M_DSTAR, M_LAMBDAC);
    kDs.set_meson_JP(1, -1);

    pseudoscalar_exchange ds_dEx (&kDs, M_D, "D exchange");
    ds_dEx.set_params({gPsiDDs, gDNL});
    ds_dEx.set_formfactor(2, M_D + eta * lambdaQCD);
    ds_dEx.force_covariant(true);

    vector_exchange ds_dstarEx (&kDs, M_DSTAR, "D* exchange");
    ds_dstarEx.set_params({gPsiDsDs, gDsNL, 0.});
    ds_dstarEx.set_formfactor(2, M_DSTAR + eta * lambdaQCD);
    ds_dstarEx.force_covariant(true);

    dirac_exchange ds_lamcEx (&kDs, M_LAMBDAC, "#Lambda_{c} exchange");
    ds_lamcEx.set_params({gPsiLL, gDsNL, 0.});
    ds_lamcEx.set_formfactor(2, M_LAMBDAC + eta * lambdaQCD);
    ds_lamcEx.force_covariant(true);

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

    double ymin = -5.5;
    double ymax = +12.0;

    std::string filename  = "psip_dslam_hpwas.pdf";
    std::string ylabel    = "#it{c}_{#{}{#it{R}}}^{1/2} (#it{s})";
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

    hpwa.set_helicities({0, 1, -1, -1});
    plotter->AddEntry(N, reF, {xmin, xmax}, "#{}{ 0 +, #minus #minus }", PRINT);
    plotter->AddDashedEntry(N, imF, {xmin, xmax}, PRINT);

    hpwa.set_helicities({0, 1, 0, +1});
    plotter->AddEntry(N, reF, {xmin, xmax}, "#{}{ 0 +, 0 + }", PRINT);
    plotter->AddDashedEntry(N, imF, {xmin, xmax}, PRINT);

    hpwa.set_helicities({0, 1, 0, -1});
    plotter->AddEntry(N, reF, {xmin, xmax}, "#{}{ 0 +, 0 #minus }", PRINT);
    plotter->AddDashedEntry(N, imF, {xmin, xmax}, PRINT);

    hpwa.set_helicities({0, 1, +1, +1});
    plotter->AddEntry(N, reF, {xmin, xmax}, "#{}{ 0 +, + + }", PRINT);
    plotter->AddDashedEntry(N, imF, {xmin, xmax}, PRINT);


    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.3, 0.75);
    plotter->SetLegendOffset(0.4, 0.15);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

    return 0;
};