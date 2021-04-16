#include "constants.hpp"
#include "amplitudes/reaction_kinematics.hpp"
#include "amplitudes/pomeron_exchange.hpp"

#include "photoPlotter.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // ---------------------------------------------------------------------------
    // AMPLITUDES
    // ---------------------------------------------------------------------------

    // Set up Kinematics for jpsi in final state
    auto * ptr = new reaction_kinematics(M_JPSI);
    ptr->set_JP(1, -1);

    // ---------------------------------------------------------------------------
    // JPAC amplitude 2016

    auto alpha1 = new linear_trajectory(1, 1.15, 0.11, "HE");
    auto background1 = new pomeron_exchange(ptr, alpha1, 1, "JPAC 2016");

    // Normalization and t-slope
    background1->set_params({0.16, 1.01});

    // ---------------------------------------------------------------------------
    // JPAC amplitude 2019

    auto alpha0 = new linear_trajectory(+1, 0.941, 0.364);
    auto background0 = new pomeron_exchange(ptr, alpha0, 0, "JPAC 2019");

    // Normalization and t-slope
    background0->set_params({0.379, 0.12});

    // ---------------------------------------------------------------------------
    // Wang et al. amplitude
  
    auto alpha2 = new linear_trajectory(+1, 1. + 0.08, 0.25);
    auto background2 = new pomeron_exchange(ptr, alpha2, 2, "Wang et al 2019");

    // Pomeron-charm coupling and cutoff
    background2->set_params({sqrt(0.8), 1.2});

    // ---------------------------- -----------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(background1);
    amps.push_back(background0);
    amps.push_back(background2);

    auto plotter = new photoPlotter(amps);

    plotter->N = 30;
    plotter->PRINT_TO_COMMANDLINE = true;
    plotter->LAB_ENERGY = true;

    plotter->xmin = E_beam(ptr->Wth()) + EPS;
    plotter->xmax = 12.;

    plotter->ymin = 0.;
    plotter->ymax = 2.;

    plotter->SHOW_LEGEND = true;
    plotter->xlegend = 0.2;
    plotter->ylegend = 0.6;

    plotter->filename  = "jpsi_compare.pdf";
    plotter->ylabel    = "#it{#sigma(#gamma p #rightarrow J/#psi p)}  [nb]";
    plotter->xlabel    = "#it{E_{#gamma}}  [GeV]";

    plotter->Plot("integrated_xsection");

    delete ptr;
    delete alpha1;
    delete alpha2;
    delete background1;
    delete background2;
    delete plotter;
    return 1;
};
