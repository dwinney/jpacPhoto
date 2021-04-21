
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
    double eta       = 1.,    lambdaQCD      = 0.25;  // Form factor parameters
    double gPsiDD    = 7.44;                          // Top couplings
    double gPsiDDstar = gPsiDD / sqrt(M_D * M_DSTAR);
    double gDNLam    = -4.3,  gDstarNLam     = -13.2; // Bottom couplings

    // ---------------------------------------------------------------------------
    // Dstar phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for Dbar LambdaC in final state
    auto kDstar = new reaction_kinematics(M_DSTAR, M_LAMBDAC, M_PROTON, M_JPSI);
    kDstar->set_JP(VECTOR);

    auto dstar_dstarEx = new vector_exchange(kDstar, M_DSTAR, "Analytic");
    dstar_dstarEx->set_params({gPsiDDstar, gDstarNLam, 0.});
    dstar_dstarEx->set_formfactor(2, M_DSTAR + eta * lambdaQCD);

    auto dstar_dstarExC = new vector_exchange(kDstar, M_DSTAR, "Covariant");
    dstar_dstarExC->set_params({gPsiDDstar, gDstarNLam, 0.});
    dstar_dstarExC->set_formfactor(2, M_DSTAR + eta * lambdaQCD);
    dstar_dstarExC->set_debug(1);

   // ---------------------------------------------------------------------------
    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(dstar_dstarEx);
    amps.push_back(dstar_dstarExC);

    auto plotter = new photoPlotter(amps);

    plotter->N = 10;
    plotter->PRINT_TO_COMMANDLINE = true;
    plotter->LAB_ENERGY = true;

    plotter->xmin = 4.;
    plotter->xmax = 10.5;

    plotter->ymin = 0.;
    plotter->ymax = 10000.;

    plotter->SHOW_LEGEND = true;
    plotter->xlegend = 0.2;
    plotter->ylegend = 0.75;
    plotter->SetLegendOffset(0.5, 0.1);

    plotter->filename  = "open_charm.pdf";
    plotter->ylabel    = "#it{#sigma(#psi p #rightarrow D #Lambda_{c}^{+})}   [nb]";
    plotter->xlabel    = "#it{W_{#gamma}}  [GeV]";

    plotter->Plot("integrated_xsection");

    return 0;
};