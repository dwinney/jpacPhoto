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
#include "amplitudes/dirac_exchange.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "one_loop/box_amplitude.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

    // ---------------------------------------------------------------------------
    // D loop
    // ---------------------------------------------------------------------------

    // Couplings
    double lambdaQCD  = 0.250, e = sqrt(4. * PI * ALPHA);
    double gGamDDstar = 0.134, gGamDstarDstar = 0.641;
    double gDNLam     = -4.3,  gDstarNLam     = -13.2;
    double gPsiDD     = 7.44;
    double gPsiDDstar = gPsiDD / sqrt(M_D * M_DSTAR);
    double gPsiLamLam = -1.4;

    // ---------------------------------------------------------------------------
    // Gamma amplitudes

    // Kinematics for the sub-processes (gamma) p -> Lam D
    auto kgamD = new reaction_kinematics(M_D, M_LAMBDAC, M_PROTON);
    kgamD->set_JP(0, -1);  // Pseudo-scalar production
    
    auto gamDDstarEx = new vector_exchange(kgamD, M_DSTAR); 
    gamDDstarEx->set_params({gGamDDstar, gDstarNLam, 0.});

    // ---------------------------------------------------------------------------
    // Psi amplitudes

    auto kpsiD = new reaction_kinematics(M_D, M_LAMBDAC, M_PROTON, M_JPSI);
    kpsiD->set_JP(0, -1);  // Pseudo-scalar production

    auto psiDDstarEx = new vector_exchange(kpsiD, M_DSTAR);
    psiDDstarEx->set_params({gPsiDDstar, gDstarNLam, 0.});

    // ---------------------------------------------------------------------------
    // Box Amplitude

    // Set up Kinematics for the overall process gamma p -> jpsi p
    auto kJPsi = new reaction_kinematics(M_JPSI, M_PROTON);
    kJPsi->set_JP(1, -1); // Vector production

    // Combine sub-processes in a box_loop
    auto disc = new box_discontinuity(gamDDstarEx, psiDDstarEx);
    auto dBox = new box_amplitude(kJPsi, disc, "Charm Loop");

    auto W_cut = [&](double qmax)
    { 
        return sqrt(qmax*qmax + M2_LAMBDAC) + sqrt(qmax*qmax + M2_D);
    };

    double qmax = 1.;
    dBox->set_cutoff(W_cut(qmax) * W_cut(qmax));

    double eta  = 1.;
    gamDDstarEx->set_formfactor(2, M_DSTAR + lambdaQCD * eta);
    psiDDstarEx->set_formfactor(2, M_DSTAR + lambdaQCD * eta);

    gamDDstarEx->set_debug(2);
    psiDDstarEx->set_debug(2);

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------
    
    int N = 10;

   // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    jpacGraph1D * plotter = new jpacGraph1D();


    double xmin = E_beam(kJPsi->Wth()) + EPS;
    double xmax = 9.2;

    auto F = [&](double egam)
    {
        double s = pow(W_cm(egam), 2.);
        return dBox->integrated_xsection(s);
    };
    
    auto x_fx = vec_fill(N, F, xmin, xmax, true);
    plotter->AddEntry(std::get<0>(x_fx), std::get<1>(x_fx), "Du et al.");

    plotter->SetLegend(0.2, 0.6);
    plotter->SetLegendOffset(0.2, 0.1);
    plotter->SetLegend(false);

    plotter->SetYaxis("#sigma(#gamma p #rightarrow J/#psi p)    [nb]", 0., 1.);
    plotter->SetXaxis("E_{#gamma}    [GeV]", xmin, xmax);
    plotter->Plot("helamp.pdf");

    delete plotter;
    delete kJPsi;
    delete kgamD;
    delete kpsiD;
    delete disc;
    delete dBox;
    delete psiDDstarEx;
    delete gamDDstarEx;

    return 1;
};