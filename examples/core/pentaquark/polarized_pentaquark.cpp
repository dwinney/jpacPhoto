// ---------------------------------------------------------------------------
// Predicted sensativity to LHCb pentaquarks in double polarization Observables
// at Hall A at JLab.
//
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

void polarized_pentaquark()
{
    // ---------------------------------------------------------------------------
    // COMMAND LINE OPTIONS
    // ---------------------------------------------------------------------------

    // Default values
    double theta = 0.;
    int N = 100; // how many points to plot
    double max = 5.;
    std::string filename = "polarized_5q.pdf";

    // ---------------------------------------------------------------------------
    // AMPLITUDES
    // ---------------------------------------------------------------------------

    // Set up Kinematics
    reaction_kinematics kpsi (M_JPSI);
    kpsi.set_meson_JP(1, -1);

    // ---------------------------------------------------------------------------
    // S - CHANNEL

    // Two different pentaquarks
    // masses and widths from 2015 LHCb paper [2]
    baryon_resonance P_c4450(&kpsi, 1, 1, 4.45, 0.040, "P_{c}(4450)");
    P_c4450.set_params({0.01, .7071}); // 2% branching fraction and equal photocouplings

    baryon_resonance P_c4380(&kpsi, 5, +1, 4.38, 0.205, "P_{c}(4380)");
    P_c4380.set_params({0.01, .7071}); // 2% branching fraction and equal photocouplings

    // ---------------------------------------------------------------------------
    // T - CHANNEL

    // Set up pomeron trajectory
    // Best fit values from [1]
    linear_trajectory alpha(+1, 0.941, 0.364, "pomeron");

    // Create amplitude with kinematics and trajectory
    pomeron_exchange background(&kpsi, &alpha, false, "Background");

    // normalization and t-slope
    background.set_params({0.379, 0.12});

    // ---------------------------------------------------------------------------
    // SUM
    // ---------------------------------------------------------------------------
    // Incoherent sum of the s and t channels
    amplitude_sum sum5q(&kpsi, {&background, &P_c4450}, "5q Sum");
    amplitude_sum sum10q(&kpsi, {&background, &P_c4450, &P_c4380}, "10q Sum");

    std::vector<amplitude*> amps = {&background, &sum5q, &sum10q};

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter objects
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // scan over energy
    for (int n = 0; n < amps.size(); n++)
    {
        // find the desired observable
        auto F = [&](double W)
        {
            double t = kpsi.t_man(W*W, theta);
            return amps[n]->K_LL(W*W, t);
        };

        std::array<std::vector<double>, 2> x_fx = vec_fill(N, F, sqrt(kpsi.sth()) + 0.01, max);
        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->get_id());
    }

    // Add a header to legend to specify the fixed energy
    std::ostringstream streamObj;
    streamObj << std::setprecision(2) << theta;
    plotter->SetLegend(0.2, 0.7, "#theta = " + streamObj.str());

    // X axis
    plotter->SetXaxis("W  (GeV)", sqrt(kpsi.sth()) + 0.01, max);

    // To change the range of the Y-axis or the position of the Legend change the arguments here
    plotter->SetYaxis("K_{LL}");

    plotter->Plot(filename);

    delete plotter;  
};
