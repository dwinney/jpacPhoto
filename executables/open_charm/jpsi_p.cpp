
#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/dirac_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "jpacGraph1D.hpp"

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
    auto kDstar = new reaction_kinematics(M_JPSI, M_PROTON, M_DSTAR, M_LAMBDAC);
    kDstar->set_meson_JP(1 ,-1);

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

    int N = 10;
    double PRINT = true;

    double xmin = 4.;
    double xmax = 10.5;

    double ymin = 0.;
    double ymax = 1.E4;

    std::string filename  = "open_charm.pdf";
    std::string ylabel    = "#it{#sigma(#psi p #rightarrow D #Lambda_{c}^{+})}   [nb]";
    std::string xlabel    = "E_{#gamma}  [GeV]";

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    jpacGraph1D * plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    for (int n = 0; n < amps.size(); n++)
    {
        auto F = [&](double e)
        {
            double w = W_cm(e);
            return amps[n]->integrated_xsection(w*w);
        };

        plotter->AddEntry(N, F, {xmin, xmax}, amps[n]->get_id(), PRINT);
    };

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetLegend(0.2, 0.75);
    plotter->SetLegendOffset(0.5, 0.1);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

    return 0;
};