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
#include "amplitudes/reaction_kinematics.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/dirac_exchange.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "box/box_amplitude.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

    double eta  = 1.;
    double qmax = 1.;
    auto s_cut = [&](double qmax, double m)
    { 
        double W_cut = sqrt(qmax*qmax + M2_LAMBDAC) + sqrt(qmax*qmax + m*m);
        return W_cut*W_cut;
    };

    // Couplings
    double lambdaQCD  = 0.250, e = sqrt(4. * PI * ALPHA);
    double gGamDDstar = 0.134, gGamDstarDstar = 0.641;
    double gDNLam     = -4.3,  gDstarNLam     = -13.2;
    double gPsiLamLam = -1.4;
    double gPsiDD     = 7.44;
    double gPsiDDstar = gPsiDD / sqrt(M_D * M_DSTAR);
    double gPsiDstarDstar = gPsiDD * (M_DSTAR / M_D);

    // Set up Kinematics for the overall process gamma p -> jpsi p
    auto kJPsi = new reaction_kinematics(M_JPSI, M_PROTON);
    kJPsi->set_JP(1, -1); // Vector production

    // ---------------------------------------------------------------------------
    // D loop
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // Gamma amplitudes

    // Kinematics for the sub-processes (gamma) p -> Lam D
    auto kgamD = new reaction_kinematics(M_D, M_LAMBDAC, M_PROTON);
    kgamD->set_JP(0, -1);  // Pseudo-scalar production
    
    auto gamD_DstarEx = new vector_exchange(kgamD, M_DSTAR, "D* exchange"); 
    gamD_DstarEx->set_params({gGamDDstar, gDstarNLam, 0.});
    gamD_DstarEx->set_formfactor(2, M_DSTAR + lambdaQCD * eta);

    auto gamD_LamcEx = new dirac_exchange(kgamD, M_LAMBDAC, "#Lambda_{c} exchange");
    gamD_LamcEx->set_params({sqrt(4.* PI * ALPHA), -4.3});
    gamD_LamcEx->set_formfactor(2, M_LAMBDAC + eta * lambdaQCD);

    auto gamD_Sum = new amplitude_sum(kgamD, {gamD_DstarEx, gamD_LamcEx});

    // ---------------------------------------------------------------------------
    // Psi amplitudes

    auto kpsiD = new reaction_kinematics(M_D, M_LAMBDAC, M_PROTON, M_JPSI);
    kpsiD->set_JP(0, -1);  // Pseudo-scalar production

    auto psiD_DEx = new pseudoscalar_exchange(kpsiD, M_D, "D exchange");
    psiD_DEx->set_params({gPsiDD, gDNLam});
    psiD_DEx->set_formfactor(2, M_D + lambdaQCD * eta);

    auto psiD_DstarEx = new vector_exchange(kpsiD, M_DSTAR, "D* exchange");
    psiD_DstarEx->set_params({gPsiDDstar, gDstarNLam, 0.});
    psiD_DstarEx->set_formfactor(2, M_DSTAR + lambdaQCD * eta);

    auto psiD_LamcEx = new dirac_exchange(kpsiD, M_LAMBDAC, "#Lambda_{c} exchange");
    psiD_LamcEx->set_params({gPsiLamLam, -4.3});
    psiD_LamcEx->set_formfactor(2, M_LAMBDAC + eta * lambdaQCD);

    auto psiD_Sum = new amplitude_sum(kpsiD, {psiD_DEx, psiD_DstarEx, psiD_LamcEx});

    // ---------------------------------------------------------------------------
    // Box Amplitude

    // Combine sub-processes in a box_loop
    auto dDisc_Full = new box_discontinuity(gamD_DstarEx, psiD_DstarEx);
    auto dBox_Full  = new box_amplitude(kJPsi, dDisc_Full, "D Loop - Full");
    dBox_Full->set_cutoff(s_cut(qmax, M_D));

    auto dDisc_Swave = new swave_discontinuity(gamD_DstarEx, psiD_DstarEx);
    auto dBox_Swave  = new box_amplitude(kJPsi, dDisc_Swave, "D Loop - S-wave");
    dBox_Swave->set_cutoff(s_cut(qmax, M_D));

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------
    std::vector<box_amplitude*> amps;
    amps.push_back(dBox_Full);
    amps.push_back(dBox_Swave); 

    int N = 50;

   // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    jpacGraph1D * plotter = new jpacGraph1D();

    double xmin = E_beam(kJPsi->Wth()) + EPS;
    double xmax = 9.2;

    for (int n = 0; n < amps.size(); n++)
    {
        gamD_DstarEx->set_debug(n);
        psiD_DstarEx->set_debug(n);

        std::cout << "\nPrinting amplitude: " << amps[n]->_identifier << "\n";
        auto F = [&](double egam)
        {
            double s = pow(W_cm(egam), 2.);
            return amps[n]->integrated_xsection(s);
        };
        
        auto x_fx = vec_fill(N, F, xmin, xmax, true);
        vec_print(x_fx[0], x_fx[1], amps[n]->_identifier + ".dat");
        plotter->AddEntry(std::get<0>(x_fx), std::get<1>(x_fx), amps[n]->_identifier);
    }

    plotter->SetLegend(0.2, 0.6);
    plotter->SetLegendOffset(0.2, 0.1);

    plotter->SetYaxis("#sigma(#gamma p #rightarrow J/#psi  p)    [nb]", 0., 1.5);
    plotter->SetXaxis("E_{#gamma}    [GeV]", xmin, xmax);
    plotter->Plot("helamp.pdf");

    return 0;
};