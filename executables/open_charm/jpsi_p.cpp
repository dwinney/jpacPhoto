
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
    // amps.push_back(dstar_dstarExC);

    double theta = 45.;
    std::array<int,4> hel = {1, 1, 1, 1};

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    jpacGraph1D * plotter = new jpacGraph1D();

    double xmin = kDstar->Wth() + EPS;
    double xmax = 5.5;

    int N = 40;

    for (int n = 0; n < amps.size(); n++)
    {
        std::cout << "\nPrinting amplitude: " << amps[n]->_identifier << "\n";
        auto F = [&] (double W)
        {
            double s = pow(W, 2.);
            double t = amps[n]->_kinematics->t_man(s, theta * DEG2RAD);

            std::complex<double> amp = amps[n]->helicity_amplitude(hel, s, t);

            return std::abs(amp);
        };
        
        auto x_fx = vec_fill(N, F, xmin, xmax, true);
        plotter->AddEntry(std::get<0>(x_fx), std::get<1>(x_fx), amps[n]->_identifier);
    }

    plotter->SetLegend(0.2, 0.6);
    plotter->SetLegendOffset(0.2, 0.1);

    // plotter->SetYaxis("#sigma(#psi p #rightarrow D^{*} #Lambda)    [mb]", 0., 10.);
    plotter->SetXaxis("W    [GeV]", xmin, xmax);
    plotter->Plot("helamp.pdf");

    return 0;
};