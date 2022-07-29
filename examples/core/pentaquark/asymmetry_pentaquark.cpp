// ---------------------------------------------------------------------------
// Predicted sensativity to LHCb pentaquarks in beam and parity asymmetries at
// GlueX at JLab
// 
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] 10.1103/PhysRevD.100.034019
// [2] 10.1103/PhysRevLett.115.072001
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "baryon_resonance.hpp"
#include "pomeron_exchange.hpp"
#include "amplitude_sum.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

void asymmetry_pentaquark()
{

    // ---------------------------------------------------------------------------
    // COMMAND LINE OPTIONS
    // ---------------------------------------------------------------------------

    // Default values
    double W = 4.45;
    int N = 200; // how many points to plot
    std::string filename = "5q_beam_asymmetry.pdf";

    // ---------------------------------------------------------------------------
    // AMPLITUDES
    // ---------------------------------------------------------------------------

    // Set up Kinematics for jpsi in final state
    reaction_kinematics kPsi (M_JPSI);
    kPsi.set_meson_JP(1, -1);

    // ---------------------------------------------------------------------------
    // T - CHANNEL 

    // Set up pomeron trajectory
    // best fit values from [1]
    linear_trajectory alpha(+1, 0.941, 0.364);

    // Create amplitude with kinematics and trajectory
    pomeron_exchange background(&kPsi, &alpha, false, "Background");

    // normalization and t-slope
    // best fit values from [1]
    background.set_params({0.379, 0.12});

    // ---------------------------------------------------------------------------
    // S - CHANNEL  // Two different pentaquarks 

    // masses and widths from 2015 LHCb paper [2]
    baryon_resonance P_c4450(&kPsi, 3, -1, 4.45, 0.040, "P_{c}(4450)");
    P_c4450.set_params({0.01, .7071});

    // 1% branching fraction and equal photocouplings for both
    baryon_resonance P_c4380(&kPsi, 5, +1, 4.38, 0.205, "P_{c}(4380)");
    P_c4380.set_params({0.01, .7071});

    // Incoherent sum of the s and t channels
    amplitude_sum sum(&kPsi, {&background, &P_c4450, &P_c4380}, "Sum");


    // ---------------------------------------------------------------------------
    // Choose which scenario to plot
    std::vector<amplitude*> amps = {&sum, &background, &P_c4450, &P_c4380};

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter objects
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // scan over theta
    for (int n = 0; n < amps.size(); n++)
    {

        auto F = [&](double theta)
        {
            double t = kPsi.t_man(W*W, theta * DEG2RAD);
            return amps[n]->beam_asymmetry_4pi(W*W, t);
        };

        plotter->AddEntry(N, F, {0., 90.}, amps[n]->get_id());
    }

    // Add a header to legend to specify the fixed energy
    std::ostringstream streamObj;
    streamObj << std::setprecision(4) << W;
    plotter->SetLegend(0.2, 0.7, "W = " + streamObj.str() + " GeV");

    plotter->SetXaxis("#theta", 0., 90.);

    // To change the range of the Y-axis or the position of the Legend change the arguments here
    plotter->SetYaxis("#Sigma");

    plotter->Plot(filename.c_str());

    delete plotter;  
};
