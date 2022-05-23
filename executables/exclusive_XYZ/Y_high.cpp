// ---------------------------------------------------------------------------
// Prediction for Y(4260) and Psi(1S and 2S) based on helicity conserving 
// pomeron exchange at high energies.
//
// Reproduces right plot in FIG 5 of [1] 
// 
// USAGE:
// make Y_high && ./Y_high
//
// OUTPUT:
// Y_HE.pdf
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
    // Same but high-energy
    linear_trajectory * alpha_HE = new linear_trajectory(1, 1.15, 0.11, "HE");
    double b_HE = 1.01;
    double A_HE = 0.16;

    // J/Psi
    reaction_kinematics * kJpsi = new reaction_kinematics(M_JPSI);
    kJpsi->set_JP(1, -1);
    double R_Jpsi = 1.;

    // Psi(2S)
    reaction_kinematics * kPsi2s = new reaction_kinematics(M_PSI2S);
    kPsi2s->set_JP(1, -1);
    double R_Psi2s = 0.55;

    // Y(4260)
    reaction_kinematics * kY = new reaction_kinematics(M_Y4260);
    kY->set_JP(1, -1);
    double R_Y = 0.84;

    // ---------------------------------------------------------------------------
    // High-enegy Amplitudes
    // ---------------------------------------------------------------------------

    pomeron_exchange Jpsi_HE(kJpsi, alpha_HE, true, "#it{J /#psi}");
    Jpsi_HE.set_params({A_HE * R_Jpsi, b_HE});

    pomeron_exchange Psi2s_HE(kPsi2s, alpha_HE, true, "#psi(2#it{S})");
    Psi2s_HE.set_params({A_HE * R_Psi2s, b_HE});

    pomeron_exchange Y_HE(kY, alpha_HE, true, "#it{Y}(4260)");
    Y_HE.set_params({A_HE * R_Y, b_HE});

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(&Jpsi_HE);
    amps.push_back(&Psi2s_HE);
    amps.push_back(&Y_HE);

    // Options
    int N = 25;
    double  xmin = 30.;
    double  xmax = 100.;

    double  ymin = 0.;
    double  ymax = 100.;

    std::string ylabel  = ROOT_italics("#sigma(#gamma p #rightarrow Y p)") + "   [nb]";
    std::string filename = "Y_HE.pdf";
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