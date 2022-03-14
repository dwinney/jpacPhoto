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
    linear_trajectory * alpha = new linear_trajectory(-1, 0.5, 0.9, "#rho - #omega");

    // ---------------------------------------------------------------------------
    // High-Energy Amplitudes
    // ---------------------------------------------------------------------------

    ////////////////////
    // Chi_c1(1P)
    vector_exchange Chi_omega(kChi, alpha, "#omega");
    Chi_omega.set_params({gChi_omega, gV_omega, gT_omega});
    Chi_omega.set_formfactor(true, LamOmega);

    vector_exchange Chi_rho(kChi, alpha, "#rho");
    Chi_rho.set_params({gChi_rho, gV_rho, gT_rho});
    Chi_rho.set_formfactor(true, LamRho);

    std::vector<amplitude*> chi_exchanges = {&Chi_omega, &Chi_rho};
    amplitude_sum chi(kChi, chi_exchanges, "#it{#chi_{c1}(1P)}");
    //////////////////

    //////////////////
    // X(3872)
    vector_exchange X_omega(kX, alpha, "#omega");
    X_omega.set_params({gX_omega, gV_omega, gT_omega});
    X_omega.set_formfactor(true, LamOmega);

    vector_exchange X_rho(kX, alpha, "#rho");
    X_rho.set_params({gX_rho, gV_rho, gT_rho});
    X_rho.set_formfactor(true, LamRho);

    std::vector<amplitude*> X_exchanges = {&X_omega, &X_rho};
    amplitude_sum X(kX, X_exchanges, "#it{X}(3872)");
    /////////////////

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(&chi);
    amps.push_back(&X);

    int N = 30;
    bool PRINT_TO_COMMANDLINE = true;

    double  xmin = 20.;
    double  xmax = 60.;

    double  ymin = 1.E-5;
    double  ymax = 1.;

    std::string filename  = "X_regge.pdf";
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
        std::cout << std::endl << "Printing amplitude: " << amps[n]->_identifier << "\n";

        auto F = [&](double x)
        {
            return amps[n]->integrated_xsection(x*x);
        };

        std::array<std::vector<double>, 2> x_fx, x_fx1;
        if (xmin < amps[n]->_kinematics->Wth())
        {
            x_fx = vec_fill(N, F, amps[n]->_kinematics->Wth() + EPS, xmax, PRINT_TO_COMMANDLINE);
            x_fx[0].insert(x_fx[0].begin(), amps[n]->_kinematics->Wth());
            x_fx[1].insert(x_fx[1].begin(), 0.);
        }
        else
        {
            x_fx = vec_fill(N, F, xmin, xmax, PRINT_TO_COMMANDLINE);
        }

        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->_identifier);
    }

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(1);
    plotter->SetLegend(0.73, 0.65);

    // Output to file
    plotter->Plot(filename);

    return 0;
}
