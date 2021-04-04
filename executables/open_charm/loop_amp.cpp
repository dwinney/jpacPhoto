// ---------------------------------------------------------------------------
// Photoproduction cross-sections of the Lambda_c Dbar/D* final state
// by open charm exchanges 
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] arXiv:2009.08345v1
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "helicities.hpp"
#include "amplitudes/reaction_kinematics.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "one_loop/box_amplitude.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // Form factor parameter
    double eta  = 1.;
    double qmax = 1.;

    // ---------------------------------------------------------------------------
    // D phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for the overall process gamma p -> jpsi p
    auto kPsi = new reaction_kinematics(M_JPSI, M_PROTON);
    kPsi->set_JP(1, -1); // Vector production

    // Kinematics for the sub-processes (gamma / psi) p -> Lam D
    auto kgamD = new reaction_kinematics(M_D, M_LAMBDAC, M_PROTON);
    kgamD->set_JP(0, -1);  // Pseudo-scalar production

    auto kpsiD = new reaction_kinematics(M_D, M_LAMBDAC, M_PROTON, M_JPSI);
    kpsiD->set_JP(0, -1);  // Pseudo-scalar production

    // Couplings
    double lambdaQCD  = 0.250, e = sqrt(4. * PI * ALPHA);
    double gGamDDstar = 0.134, gGamDstarDstar = 0.641;
    double gDNLam     = -4.3,  gDstarNLam     = -13.2;
    double gPsiDD     = 7.44;
    double gPsiDDstar = gPsiDD / sqrt(M_D * M_DSTAR);

    // Vector exchange loop
    auto gamDDstarEx = new vector_exchange(kgamD, M_DSTAR); 
    gamDDstarEx->set_params({gGamDDstar, gDstarNLam, 0.});
    gamDDstarEx->set_formfactor(2, M_DSTAR + lambdaQCD * eta);

    auto psiDDstarEx = new vector_exchange(kpsiD, M_DSTAR);
    psiDDstarEx->set_params({gPsiDDstar, gDstarNLam, 0.});
    psiDDstarEx->set_formfactor(2, M_DSTAR + lambdaQCD * eta);

    // Combine sub-processes in a box_loop
    auto disc = new box_discontinuity(gamDDstarEx, psiDDstarEx);
    auto ddstar_box = new box_amplitude(kPsi, disc, "Charm Loop");

    double W_cut = sqrt(qmax*qmax + M2_LAMBDAC) + sqrt(qmax*qmax + M2_D);
    ddstar_box->set_cutoff(W_cut * W_cut);

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    auto alpha = new linear_trajectory(+1, 0.941, 0.364);
    auto pomeron = new pomeron_exchange(kPsi, alpha, 0, "Pomeron");
    pomeron->set_params({0.379, 0.12});

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(ddstar_box);
    amps.push_back(pomeron);

        // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------
    
    int N = 10;

    double xmin = 0.;
    double xmax = 90.;

    double theta = 0.;

   // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    std::cout << std::endl << "Printing amplitude: " << amps[n]->_identifier << "\n";

    auto F = [&](double x)
    {
        double W = W_cm(8.7);
        double s = W*W;
        double t = kPsi->t_man(s, x * DEG2RAD);
        
        return amps[n]->beam_asymmetry_4pi(s, t);
    };

    std::array<std::vector<double>, 2> x_fx;
    x_fx = vec_fill(N, F, xmin, xmax, true);

    return 1;
};