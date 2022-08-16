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

void psip_dlam_PWA()
{
    // Form factor parameter
    double eta = 1.;
    double lambdaQCD = 0.25;

    double gPsiDD  = 7.4;
    double gPsiDDs = gPsiDD / sqrt(M_D * M_DSTAR);
    double gDNL    = -4.3;
    double gDsNL   = -13.2;
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

    vector_exchange d_dstarEx (&kD, M_DSTAR, "D^{*} exchange");
    d_dstarEx.set_params({gPsiDDs, gDsNL, 0.});
    d_dstarEx.set_formfactor(2, M_DSTAR + eta * lambdaQCD);
    d_dstarEx.force_covariant(true);

    dirac_exchange d_lamcEx (&kD, M_LAMBDAC, "#Lambda_{c} exchange");
    d_lamcEx.set_params({sqrt(4.* PI * ALPHA), gPsiLL, 0.});
    d_lamcEx.set_formfactor(2, M_LAMBDAC + eta * lambdaQCD);
    d_lamcEx.force_covariant(true);

    amplitude_sum d_sum (&kD,  {&d_dEx, &d_dstarEx, &d_lamcEx}, "Sum");

    // ---------------------------------------------------------------------------
    // PW projection
    // ---------------------------------------------------------------------------
    
    // Take the sum ampitude and pass it to a projected_amplitude
    projected_amplitude d_sum1(&d_sum, 1, "#it{J} = 1/2");
    projected_amplitude d_sum3(&d_sum, 3, "#it{J} = 3/2");
    projected_amplitude d_sum5(&d_sum, 5, "#it{J} = 5/2");

    amplitude_sum d_sum135 (&kD,  {&d_sum1, &d_sum3, &d_sum5}, "Sum up to #it{J}_{max} = 5/2");

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(&d_sum1);
    amps.push_back(&d_sum3);
    amps.push_back(&d_sum5);
    amps.push_back(&d_sum135);
    amps.push_back(&d_sum);

    int N = 50;
    double PRINT = true;

    double xmin = 4.;
    double xmax = 6.;

    double ymin = 0.;
    double ymax = 200.;

    std::string filename  = "open_charm.pdf";
    std::string ylabel    = "#sigma(#psi #it{p} #rightarrow #bar{#it{D}} #Lambda_{c}^{+})   [#mub]";
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
            return amps[n]->integrated_xsection(w*w) * 1E-3;
        };

        plotter->AddEntry(N, F, {xmin, xmax}, amps[n]->get_id(), PRINT);
    };

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.2, 0.65);
    plotter->SetLegendOffset(0.5, 0.17);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

    return 0;
};