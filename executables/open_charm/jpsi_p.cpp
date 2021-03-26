
#include "constants.hpp"
#include "amplitudes/reaction_kinematics.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/dirac_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "photoPlotter.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // Form factor parameter
    double eta = 1.;

    // ---------------------------------------------------------------------------
    // D phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for Dbar LambdaC in final state
    auto kD = new reaction_kinematics(M_D, M_LAMBDAC, M_PROTON, M_JPSI);
    kD->set_JP(PSEUDO_SCALAR);

    auto d_dstarEx = new vector_exchange(kD, M_DSTAR, "FF");
    d_dstarEx->set_params({0.001, -13.2, 0.});
    d_dstarEx->set_formfactor(2, M_DSTAR + eta * 0.250);

    auto d_dstarExC = new vector_exchange(kD, M_DSTAR, "no FF");
    d_dstarExC->set_params({0.01, -13.2, 0.});
    d_dstarExC->set_formfactor(0, M_DSTAR + eta * 0.250);
    // d_dstarExC->set_debug(1);

   // ---------------------------------------------------------------------------
    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(d_dstarEx);
    amps.push_back(d_dstarExC);

    auto plotter = new photoPlotter(amps);

    plotter->N = 30;
    plotter->PRINT_TO_COMMANDLINE = true;
    plotter->LAB_ENERGY = false;

    plotter->xmin = 4.;
    plotter->xmax = 10.5;

    plotter->ymin = 0.;
    plotter->ymax = 540.;

    plotter->SHOW_LEGEND = true;
    plotter->xlegend = 0.2;
    plotter->ylegend = 0.75;
    plotter->SetLegendOffset(0.5, 0.1);

    plotter->filename  = "open_charm.pdf";
    plotter->ylabel    = "#it{#sigma(#psi p #rightarrow D #Lambda_{c}^{+})}   [nb]";
    plotter->xlabel    = "#it{W_{#gamma}}  [GeV]";

    plotter->Plot("integrated_xsection");

    return 1;
};