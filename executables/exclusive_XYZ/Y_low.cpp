// ---------------------------------------------------------------------------
// Prediction for Y(4260) and Psi(1S and 2S) based on effective vector
// pomeron exchange at low enegies.
//
// Reproduces left plot in FIG 5 of [1] 
// 
// USAGE:
// make Y_low && ./Y_low
//
// OUTPUT:
// Y_LE.pdf
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
#include "regge_trajectory.hpp"
#include "amplitudes/pomeron_exchange.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

    // ---------------------------------------------------------------------------
    // Preliminaries
    // ---------------------------------------------------------------------------

    // Low - energy trajectory and couplings
    linear_trajectory * alpha_LE = new linear_trajectory(1, 0.94, 0.36, "LE");
    double b_LE = 0.12;
    double A_LE = 0.38;

    // J/Psi
    reaction_kinematics * kJpsi = new reaction_kinematics(M_JPSI);
    kJpsi->set_JP(1, -1);
    double R_Jpsi = 1.;

    // Psi(2S)
    reaction_kinematics * kPsi2s = new reaction_kinematics(M_PSI2S);
    kPsi2s->set_JP(1, -1);
    double R_Psi2s = 0.55;

    // Y(4260)
    double mY = 4.220;
    reaction_kinematics * kY = new reaction_kinematics(M_Y4260);
    kY->set_JP(1, -1);
    double R_Y = 0.84;

    // ---------------------------------------------------------------------------
    // Low-Energy Amplitudes
    // ---------------------------------------------------------------------------

    pomeron_exchange Jpsi_LE(kJpsi, alpha_LE, false, "#it{J /#psi}");
    Jpsi_LE.set_params({A_LE * R_Jpsi, b_LE});

    pomeron_exchange Psi2s_LE(kPsi2s, alpha_LE, false, "#psi(2#it{S})");
    Psi2s_LE.set_params({A_LE * R_Psi2s, b_LE});

    pomeron_exchange Y_LE(kY, alpha_LE, false, "#it{Y}(4260)");
    Y_LE.set_params({A_LE * R_Y, b_LE});

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(&Jpsi_LE);
    amps.push_back(&Psi2s_LE);
    amps.push_back(&Y_LE);

    // Options
    int N = 25;
    double  xmin = 4.;
    double  xmax = 10.;

    double  ymin = 0.;
    double  ymax = 20.;

    std::string filename = "Y_LE.pdf";
    std::string ylabel  = ROOT_italics("#sigma(#gamma p #rightarrow Y p)") + "   [nb]";
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

        plotter->AddEntry(N, F, {xmin,xmax}, amps[n]->get_id(), PRINT);
    }
    
    plotter->SetXaxis(ROOT_italics("W_{#gammap}") + "  [GeV]", std::floor(xmin), xmax);

    plotter->SetYaxis(ylabel, ymin, ymax);
    
    plotter->SetLegend(0.2, 0.73);

    // Output to file
    plotter->Plot(filename);

    return 0;
};