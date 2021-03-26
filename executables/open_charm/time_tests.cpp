
#include "constants.hpp"
#include "amplitudes/reaction_kinematics.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/dirac_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "photoPlotter.hpp"

#include <cstring>
#include <ctime>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // Form factor parameter
    double eta = 1.;

    // ---------------------------------------------------------------------------
    // D* phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for D* LambdaC in final state
    auto kDstar = new reaction_kinematics(M_DSTAR, M_LAMBDAC, M_PROTON);
    kDstar->set_JP(VECTOR);

    auto dstar_dstarEx = new vector_exchange(kDstar, M_DSTAR, "D^{*} (Analytic)");
    dstar_dstarEx->set_params({0.641, -13.2, 0.});
    dstar_dstarEx->set_formfactor(2, M_DSTAR + eta * 0.250);
    dstar_dstarEx->set_debug(0);
    
    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(dstar_dstarEx);

    auto plotter = new photoPlotter(amps);

    plotter->N = 1000;
    plotter->PRINT_TO_COMMANDLINE = true;
    plotter->LAB_ENERGY = true;

    plotter->xmin = 8.5;
    plotter->xmax = 10.5;

    plotter->ymin = 0.;
    plotter->ymax = 250.;

    plotter->SHOW_LEGEND = true;
    plotter->xlegend = 0.2;
    plotter->ylegend = 0.6;

    plotter->filename  = "open_charm.pdf";
    plotter->ylabel    = "#it{#sigma(#gamma p #rightarrow #bar{D} #Lambda_{c}^{+})}  [nb]";
    plotter->xlabel    = "#it{E_{#gamma}}  [GeV]";

    clock_t begin = clock();
    plotter->Plot("integrated_xsection");
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Done in " << elapsed_secs << " seconds. \n";
    std::cout << "\n";

    return 1;
};