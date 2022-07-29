// ---------------------------------------------------------------------------
// Prediction for X(6900) photoproduction near threshold based on hypothetical
// omega exchange amplitudes.
// 
// Reproduces the plot in FIG 4 of [1].
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// REFERENCES:
// [1] arXiv:2008.01001 [hep-ph]
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "vector_exchange.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

void X6900()
{

    // ---------------------------------------------------------------------------
    // Amplitude
    // ---------------------------------------------------------------------------

    // X(6900)
    reaction_kinematics kX (6.900);
    kX.set_meson_JP(0, 1);

    double gV_psi = 1.6E-3, gT_psi = 0.;
    double gX_psi = 5.03;
    double gGam_psi = gX_psi * E * F_JPSI / M_JPSI;

    double gV_omega = 16., gT_omega = 0.;
    double gX_omega = 0.225;
    double gGam_omega = gX_omega * E * F_JPSI / M_JPSI;

    vector_exchange X_psi(&kX, M_JPSI, "J/#psi exchange, BR = 100%");
    X_psi.set_params({gGam_psi, gV_psi, gT_psi});

    vector_exchange X_omega(&kX, M_OMEGA, "#it{X}(6900) with BR(#it{X #rightarrow #psi#omega}) = 1%");
    X_omega.set_params({gGam_omega, gV_omega, gT_omega});
    X_omega.set_formfactor(true, 1.2);
   
    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    int N = 100;

    double  xmin = 7.5;
    double  xmax = 15.;

    double  ymin = 1.E-2;
    double  ymax = 40.;

    std::string filename = "omega_exchange.pdf";
    std::string xlabel   = "#it{W}_{#gammap} [GeV]";
    std::string ylabel   = "#sigma(#gamma #it{p} #rightarrow #it{X p})   [nb]";
    bool PRINT = true;

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude

    auto F = [&](double x)
    {
        return X_omega.integrated_xsection(x*x);
    };

    plotter->AddEntry(N, F, {xmin, xmax}, X_omega.get_id(), PRINT);
    
    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(true);
    
    plotter->SetLegend(0.2, 0.75);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
};