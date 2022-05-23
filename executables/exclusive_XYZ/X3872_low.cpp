// ---------------------------------------------------------------------------
// Prediction for X(3872) and chi_c1(1P) photoproduction at low enegies
// Reproduces the plot in FIG 3a of [1].
// 
// USAGE:
// make X3872_low && ./X3872_low
//
// OUTPUT:
// X_FS.pdf
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// REFERENCES:
// [1] arXiv:2008.01001 [hep-ph]
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "amplitudes/reaction_kinematics.hpp"
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
    reaction_kinematics * kChi = new reaction_kinematics(M_CHIC1);
    kChi->set_JP(1, 1);

    // X(3872)
    reaction_kinematics * kX = new reaction_kinematics(M_X3872);
    kX->set_JP(1, 1);

    // Nucleon couplings 
    double gV_omega = 16., gT_omega = 0.;
    double LamOmega = 1.2; // cutoff
    double gV_rho = 2.4, gT_rho = 14.6;
    double LamRho = 1.4; // cutoff
    double gV_phi = -6.2, gT_phi = 2.1;
    double gV_psi = 1.6E-3, gT_psi = 0.;

    // Photon couplings
    double gChi_omega   = 5.2E-4;
    double gChi_rho     = 9.2E-4;
    double gChi_phi     = 4.2E-4;
    double gChi_psi     = 1.;
    double gX_omega     = 8.2E-3;
    double gX_rho       = 3.6E-3;

    // ---------------------------------------------------------------------------
    // Low-Energy Amplitudes
    // ---------------------------------------------------------------------------

    // Chi_c1(1P)
    vector_exchange Chi_omega(kChi, M_OMEGA, "#omega");
    Chi_omega.set_params({gChi_omega, gV_omega, gT_omega});
    Chi_omega.set_formfactor(1, LamOmega);

    vector_exchange Chi_rho(kChi, M_RHO, "#rho");
    Chi_rho.set_params({gChi_rho, gV_rho, gT_rho});
    Chi_rho.set_formfactor(1, LamRho);

    vector_exchange Chi_phi(kChi, M_PHI, "#phi");
    Chi_phi.set_params({gChi_phi, gV_phi, gT_phi});

    vector_exchange Chi_psi(kChi, M_JPSI, "#psi");
    Chi_psi.set_params({gChi_psi, gV_psi, gT_psi});

    std::vector<amplitude*> chi_exchanges = {&Chi_omega, &Chi_rho, &Chi_phi, &Chi_psi};
    amplitude_sum chi(kChi, chi_exchanges, "#it{#chi_{c1}(1P)}");

    // X(3872)
    vector_exchange X_omega(kX, M_OMEGA, "#omega");
    X_omega.set_params({gX_omega, gV_omega, gT_omega});
    X_omega.set_formfactor(1, LamOmega);

    vector_exchange X_rho(kX, M_RHO, "#rho");
    X_rho.set_params({gX_rho, gV_rho, gT_rho});
    X_rho.set_formfactor(1, LamRho);

    std::vector<amplitude*> X_exchanges = {&X_omega, &X_rho};
    amplitude_sum X(kX, X_exchanges, "#it{X}(3872)");

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(&chi);
    amps.push_back(&X);

    int N = 100;
    bool PRINT = true;

    double  xmin = 4.;
    double  xmax = 7.;

    double  ymin = 2.E-3;
    double  ymax = 800.;

    std::string filename  = "X_FS.pdf";
    std::string ylabel    = "#it{#sigma(#gamma p #rightarrow X p)}  [nb]";
    std::string xlabel    = "#it{W_{#gammap}}  [GeV]";

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

        plotter->AddEntry(N, F, {xmin,xmax}, amps[n]->get_id(), PRINT);
    }

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(1);
    plotter->SetLegend(0.7, 0.2);

    // Output to file
    plotter->Plot(filename);

    return 0;
}
