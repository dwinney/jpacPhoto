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

void psip_dlam_HPWAs()
{
    // Form factor parameter
    double eta = 1.;
    double lambdaQCD = 0.25;

    double gPsiDD  = 7.4;
    double gPsiDDs = gPsiDD / sqrt(M_DSTAR * M_D);
    double gDNL    = -13.2;
    double gDsNL   = -4.3;
    double gPsiLL  = -1.4;

    // ---------------------------------------------------------------------------
    // psi p -> D Lambda amplitudes
    // ---------------------------------------------------------------------------

    // Set up Kinematics for Dbar LambdaC in final state
    reaction_kinematics kD (M_JPSI, M_PROTON, M_D, M_LAMBDAC);
    kD.set_meson_JP(0, -1);

    pseudoscalar_exchange d_dEx (&kD, M_D, "D exchange");
    d_dEx.set_params({gPsiDD, gDNL});
    d_dEx.set_formfactor(2, M_D + eta * lambdaQCD);
    d_dEx.force_covariant(true);

    vector_exchange d_dstarEx (&kD, M_DSTAR, "D* exchange");
    d_dstarEx.set_params({gPsiDDs, gDsNL, 0.});
    d_dstarEx.set_formfactor(2, M_DSTAR + eta * lambdaQCD);
    d_dstarEx.force_covariant(true);

    dirac_exchange d_lamcEx (&kD, M_LAMBDAC, "#Lambda_{c} exchange");
    d_lamcEx.set_params({gPsiLL, gDNL, 0.});
    d_lamcEx.set_formfactor(2, M_LAMBDAC + eta * lambdaQCD);
    d_lamcEx.force_covariant(true);

    amplitude_sum d_sum (&kD,  {&d_dEx, &d_dstarEx, &d_lamcEx}, "Sum");

    // ---------------------------------------------------------------------------
    // PW projection
    // ---------------------------------------------------------------------------
    
    // Take the sum ampitude and pass it to a projected_amplitude
    helicity_PWA hpwa(&d_sum, 1);

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    int N = 50;
    double PRINT = true;

    double xmin = sqrt(kD.sth()) + 0.01;
    double xmax = 5.;

    double ymin = -4.;
    double ymax = +6.0;

    std::string filename  = "psip_dlam_hpwas.pdf";
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

    hpwa.set_helicities({1, 1, 0, +1});
    plotter->AddEntry(N, reF, {xmin, xmax}, "#{}{ + +, 0 + }", PRINT);
    plotter->AddDashedEntry(N, imF, {xmin, xmax}, PRINT);

    hpwa.set_helicities({1, 1, 0, -1});
    plotter->AddEntry(N, reF, {xmin, xmax}, "#{}{ + +, 0 #minus }", PRINT);
    plotter->AddDashedEntry(N, imF, {xmin, xmax}, PRINT);

    hpwa.set_helicities({0, 1, 0, +1});
    plotter->AddEntry(N, reF, {xmin, xmax}, "#{}{ 0 +, 0 + }", PRINT);
    plotter->AddDashedEntry(N, imF, {xmin, xmax}, PRINT);

    hpwa.set_helicities({0, 1, 0, -1});
    plotter->AddEntry(N, reF, {xmin, xmax}, "#{}{ 0 +, 0 #minus }", PRINT);
    plotter->AddDashedEntry(N, imF, {xmin, xmax}, PRINT);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.45, 0.74);
    plotter->SetLegendOffset(0.5, 0.15);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

    return 0;
};