// ---------------------------------------------------------------------------
// Prediction for X(3872) and chi_c1(1P) photoproduction at high enegies
// Reproduces the plot in FIG 3b of [1].
// 
// USAGE:
// make X3872_high && ./X3872_high
//
// OUTPUT:
// X_regge.pdf
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
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

    // ---------------------------------------------------------------------------
    // Preliminaries
    // ---------------------------------------------------------------------------

    // Chi_c1(1P)
    reaction_kinematics kChi (M_CHIC1);
    kChi.set_meson_JP(1, +1);

    // X(3872)
    reaction_kinematics kX (M_X3872);
    kX.set_meson_JP(1, +1);

    // Nucleon couplings and cutoffs
    double gV_omega = 16., gT_omega = 0.;
    double LamOmega = 1.2;
    double gV_rho = 2.4, gT_rho = 14.6;
    double LamRho = 1.4;
    double gV_phi = -6.2, gT_phi = 2.1;
    double gV_psi = 1.6E-3, gT_psi = 0.;

    // Top couplings
    double gChi_omega = 5.2E-4;
    double gChi_rho = 9.2E-4;
    double gX_omega = 8.2E-3;
    double gX_rho = 3.6E-3;
    
    // Linear trajectory for the rho
    linear_trajectory alpha (-1, 0.5, 0.9, "#rho - #omega");

    // ---------------------------------------------------------------------------
    // High-Energy Amplitudes
    // ---------------------------------------------------------------------------

    // Chi_c1(1P)
    vector_exchange Chi_omega(&kChi, &alpha, "#omega");
    Chi_omega.set_params({gChi_omega, gV_omega, gT_omega});
    Chi_omega.set_formfactor(true, LamOmega);

    vector_exchange Chi_rho(&kChi, &alpha, "#rho");
    Chi_rho.set_params({gChi_rho, gV_rho, gT_rho});
    Chi_rho.set_formfactor(true, LamRho);

    // sum the amplitudes together
    amplitude_sum chi(&kChi, {&Chi_omega, &Chi_rho}, "#chi_{c1}(1P)");

    // X(3872)
    vector_exchange X_omega(&kX, &alpha, "#omega");
    X_omega.set_params({gX_omega, gV_omega, gT_omega});
    X_omega.set_formfactor(true, LamOmega);

    vector_exchange X_rho(&kX, &alpha, "#rho");
    X_rho.set_params({gX_rho, gV_rho, gT_rho});
    X_rho.set_formfactor(true, LamRho);

    // sum the amplitudes together
    amplitude_sum X(&kX, {&X_omega, &X_rho}, "#it{X}(3872)");

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(&chi);
    amps.push_back(&X);

    int N = 30;

    double  xmin = 20.;
    double  xmax = 60.;

    double  ymin = 1.E-5;
    double  ymax = 1.;

    std::string filename  = "X_regge.pdf";
    std::string ylabel    = "#sigma(#gamma #it{p} #rightarrow #it{X p})  [nb]";
    std::string xlabel    = "#it{W}_{#gammap}  [GeV]";

    bool PRINT = true;

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    for (int n = 0; n < amps.size(); n++)
    {
        auto F = [&](double x)
        {
            return amps[n]->integrated_xsection(x*x);
        };

        plotter->AddEntry(N, F, {xmin, xmax}, amps[n]->get_id(), PRINT);
    };

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(1);
    plotter->SetLegend(0.73, 0.65);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

    return 0;
}
