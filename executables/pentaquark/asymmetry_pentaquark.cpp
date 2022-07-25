// ---------------------------------------------------------------------------
// Predicted sensativity to LHCb pentaquarks in beam and parity asymmetries at
// GlueX at JLab
// 
// USAGE:
// make asymmetry_pentaquark && ./asymmetry_pentaquark -optional_flags flag_inputs
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] 10.1103/PhysRevD.100.034019
// [2] 10.1103/PhysRevLett.115.072001
// ---------------------------------------------------------------------------
// COMMAND LINE OPTIONS:
// -f string          # Select the filename of output pdf (default: "5q_beam_asymmetry.pdf")
// -e double          # Change CM fixed-energy (default: 4.45 Gev)
// -10q               # Plot 2 Pentaquark Scenario at fixed BR (default: false)
// -y "[y1:y2]"       # Custom y bounds in output plot
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/baryon_resonance.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

    // ---------------------------------------------------------------------------
    // COMMAND LINE OPTIONS
    // ---------------------------------------------------------------------------

    // Default values
    double y[2]; bool custom_y = false;
    double W = 4.45;
    int N = 200; // how many points to plot
    bool TENQ = false;
    std::string filename = "5q_beam_asymmetry.pdf";

    // Parse input string
    for (int i = 0; i < argc; i++)
    {
        if (std::strcmp(argv[i],"-e")==0) W = atof(argv[i+1]);
        if (std::strcmp(argv[i],"-f")==0) filename = argv[i+1];
        if (std::strcmp(argv[i],"-10q")==0) TENQ = true;
        if (std::strcmp(argv[i],"-y")==0)
        {
            custom_y = true;
            y_range(argv[i+1], y);
        }
    }

    // ---------------------------------------------------------------------------
    // AMPLITUDES
    // ---------------------------------------------------------------------------

    // Set up Kinematics for jpsi in final state
    reaction_kinematics kPsi (M_JPSI);
    kPsi.set_meson_JP(1, -1);

    // ---------------------------------------------------------------------------
    // T - CHANNEL // this is the same for all cases

    // Set up pomeron trajectory
    // best fit values from [1]
    linear_trajectory alpha(+1, 0.941, 0.364);

    // Create amplitude with kinematics and trajectory
    pomeron_exchange background(&kPsi, &alpha, false, "Background");

    // normalization and t-slope
    // best fit values from [1]
    background.set_params({0.379, 0.12});

    // ---------------------------------------------------------------------------
    // S - CHANNEL  // Two different pentaquarks (if -10q == true)

    // masses and widths from 2015 LHCb paper [2]
    baryon_resonance P_c4450(&kPsi, 3, -1, 4.45, 0.040, "P_{c}(4450)");
    P_c4450.set_params({0.01, .7071});

    // 1% branching fraction and equal photocouplings for both
    baryon_resonance P_c4380(&kPsi, 5, +1, 4.38, 0.205, "P_{c}(4380)");
    P_c4380.set_params({0.01, .7071});

    // Incoherent sum of the s and t channels
    amplitude_sum sum(&kPsi, {&background, &P_c4450, &P_c4380}, "Sum");

    // ---------------------------------------------------------------------------
    // S - CHANNEL  // 1 5q but different BR scenarios (if -10q == false)

    baryon_resonance P_c1(&kPsi, 3, -1, 4.45, 0.040, "1%");
    P_c1.set_params({0.01, .7071});

    baryon_resonance P_c05(&kPsi, 3, -1, 4.45, 0.040, "0.5%");
    P_c05.set_params({0.005, .7071});

    baryon_resonance P_c01(&kPsi, 3, -1, 4.45, 0.040, "0.1%");
    P_c01.set_params({0.001, .7071});

    // Add to the sum
    amplitude_sum sum1(&kPsi, {&background, &P_c1}, "1%");
    amplitude_sum sum2(&kPsi, {&background, &P_c05}, "0.5%");
    amplitude_sum sum3(&kPsi, {&background, &P_c01}, "0.1%");


    // ---------------------------------------------------------------------------
    // Choose which scenario to plot
    std::vector<amplitude*> amps;
    if (TENQ == true)
    {
        amps = {&sum, &background, &P_c4450, &P_c4380};
    }
    else
    {
        amps = {&background, &sum1, &sum2, &sum3};
    }

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
    (custom_y == true) ? (plotter->SetYaxis("#Sigma", y[0], y[1])) : (plotter->SetYaxis("#Sigma"));


    plotter->Plot(filename.c_str());

    delete plotter;
    
    return 0;
};
