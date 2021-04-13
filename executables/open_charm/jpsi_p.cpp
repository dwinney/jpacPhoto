
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
    // D phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for Dbar LambdaC in final state
    auto kD = new reaction_kinematics(M_D, M_LAMBDAC, M_PROTON);
    kD->set_JP(VECTOR);

    auto d_dEx = new pseudoscalar_exchange(kD, M_D, "Analytic");
    d_dEx->set_params({gPsiDD, gDNLam});
    d_dEx->set_formfactor(2, M_D + eta * lambdaQCD);

    auto d_dExC = new pseudoscalar_exchange(kD, M_D, "Covariant");
    d_dExC->set_params({gPsiDD, gDNLam});
    d_dExC->set_formfactor(2, M_D + eta * lambdaQCD);
    d_dExC->set_debug(1);

    // auto d_dstarEx = new vector_exchange(kD, M_DSTAR, "D* exchange");
    // d_dstarEx->set_params({gPsiDDstar, gDstarNLam, 0.});
    // d_dstarEx->set_formfactor(2, M_DSTAR + eta * lambdaQCD);

   // ---------------------------------------------------------------------------
    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(d_dEx);
    amps.push_back(d_dExC);
    // amps.push_back(d_dstarEx);

    auto plotter = new photoPlotter(amps);

    plotter->N = 30;
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

    return 1;
};