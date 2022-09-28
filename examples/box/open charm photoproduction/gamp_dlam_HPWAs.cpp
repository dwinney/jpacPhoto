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

void gamp_dlam_HPWAs()
{
    // Form factor parameter
    double eta = 1.;
    double lambdaQCD = 0.25;

    double e     = sqrt(4.* PI * ALPHA);
    double gDDs  = 0.134;
    double gNDsL = -4.3;
    double gNDL  = -13.2;

    // ---------------------------------------------------------------------------
    // D phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for Dbar LambdaC in final state
    reaction_kinematics kD (M_D, M_LAMBDAC);
    kD.set_meson_JP(0, -1);

    vector_exchange d_dstarEx (&kD, M_DSTAR, "D* exchange");
    d_dstarEx.set_params({gDDs, gNDsL, 0.});
    d_dstarEx.set_formfactor(2, M_DSTAR + eta * lambdaQCD);
    d_dstarEx.force_covariant(true);

    dirac_exchange d_lamcEx (&kD, M_LAMBDAC, "#Lambda_{c} exchange");
    d_lamcEx.set_params({e, gNDL, 0.});
    d_lamcEx.set_formfactor(2, M_LAMBDAC + eta * lambdaQCD);
    d_lamcEx.force_covariant(true);

    amplitude_sum d_sum (&kD,  {&d_dstarEx, &d_lamcEx}, "Sum");

    // ---------------------------------------------------------------------------
    // PW projection
    // ---------------------------------------------------------------------------
    
    // Take the sum ampitude and pass it to a projected_amplitude
    int J = 1;
    helicity_PWA hpwa(&d_sum, J);

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    int N = 100;
    double PRINT = true;

    double xmin = sqrt(kD.sth() + 1.E-4);
    double xmax = 5.;

    double ymin = -1.;
    double ymax = 1.3;

    std::string filename  = "gamp_pwa.pdf";
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

    hpwa.set_helicities({1, 1, 0, +1});
    plotter->AddEntry(N, reF, {xmin, xmax}, "#{}{ + +, 0 + }", PRINT);
    hpwa.set_helicities({-1, -1, 0, -1});
    plotter->AddDashedEntry(N, reF, {xmin, xmax}, PRINT);

    hpwa.set_helicities({1, 1, 0, -1});
    plotter->AddEntry(N, reF, {xmin, xmax}, "#{}{ + +, 0 #minus }", PRINT);
    hpwa.set_helicities({-1, -1, 0, +1});
    plotter->AddDashedEntry(N, reF, {xmin, xmax}, PRINT);

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.2, 0.5);
    plotter->SetLegendOffset(0.5, 0.1);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

    return 0;
};