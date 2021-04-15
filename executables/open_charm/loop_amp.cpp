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
#include "amplitudes/pseudoscalar_exchange.hpp"
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

    auto psiD_Sum = new amplitude_sum(kpsiD, {psiD_DEx, psiD_DstarEx});

    // ---------------------------------------------------------------------------
    // Box Amplitude

    // Combine sub-processes in a box_loop
    auto dDisc = new box_discontinuity(gamD_DstarEx, psiD_Sum);
    auto dBox  = new box_amplitude(kJPsi, dDisc, "D Loop");
    dBox->set_cutoff(s_cut(qmax, M_D));

    // ---------------------------------------------------------------------------
    // Dstar loop
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // Gamma amplitudes

    // Kinematics for the sub-processes gamma p -> Lam Dstar
    auto kgamDstar = new reaction_kinematics(M_DSTAR, M_LAMBDAC, M_PROTON);
    kgamDstar->set_JP(1, -1);  // Vector 
    
    auto gamDstar_DEx = new pseudoscalar_exchange(kgamDstar, M_D);
    gamDstar_DEx->set_params({gGamDDstar, gDstarNLam});
    gamDstar_DEx->set_formfactor(2, M_D + lambdaQCD * eta);

    auto gamDstar_DstarEx = new vector_exchange(kgamDstar, M_DSTAR); 
    gamDstar_DstarEx->set_params({gGamDstarDstar, gDstarNLam, 0.});
    gamDstar_DstarEx->set_formfactor(2, M_DSTAR + lambdaQCD * eta);

    auto gamDstar_Sum = new amplitude_sum(kgamDstar, {gamDstar_DEx, gamDstar_DstarEx});

    // ---------------------------------------------------------------------------
    // Psi amplitudes

    // Kinematics for the sub-processes psi p -> Lam Dstar
    auto kpsiDstar = new reaction_kinematics(M_DSTAR, M_LAMBDAC, M_PROTON, M_JPSI);
    kpsiDstar->set_JP(1, -1); 

    auto psiDstar_DEx = new pseudoscalar_exchange(kpsiDstar, M_D);
    psiDstar_DEx->set_params({gPsiDDstar, gDstarNLam});
    psiDstar_DEx->set_formfactor(2, M_D + lambdaQCD * eta);

    auto psiDstar_DstarEx = new vector_exchange(kpsiDstar, M_DSTAR); 
    psiDstar_DstarEx->set_params({gPsiDstarDstar, gDstarNLam, 0.});
    psiDstar_DstarEx->set_formfactor(2, M_DSTAR + lambdaQCD * eta);

    auto psiDstar_Sum = new amplitude_sum(kpsiDstar, {psiDstar_DEx, psiDstar_DstarEx});

    // ---------------------------------------------------------------------------
    // Box Amplitude

    // Combine sub-processes in a box_loop
    auto dstarDisc = new box_discontinuity(gamDstar_Sum, psiDstar_Sum);
    auto dstarBox  = new box_amplitude(kJPsi, dstarDisc, "D* Loop");
    dstarBox->set_cutoff(s_cut(qmax, M_DSTAR));

    // ---------------------------------------------------------------------------
    // Combine both box loops

    auto combined_box = new box_amplitude(kJPsi, {dDisc, dstarDisc}, "Sum");

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------
    std::vector<box_amplitude*> amps;
    amps.push_back(dBox);
    amps.push_back(dstarBox);
    amps.push_back(combined_box);

    int N = 10;

   // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    jpacGraph1D * plotter = new jpacGraph1D();

    double xmin = E_beam(kJPsi->Wth()) + EPS;
    double xmax = 10.;

    for (int n = 0; n < amps.size(); n++)
    {
        std::cout << "\nPrinting amplitude: " << amps[n]->_identifier << "\n";
        auto F = [&](double egam)
        {
            double s = pow(W_cm(egam), 2.);
            return amps[n]->integrated_xsection(s);
        };
        
        auto x_fx = vec_fill(N, F, xmin, xmax, true);
        plotter->AddEntry(std::get<0>(x_fx), std::get<1>(x_fx), amps[n]->_identifier);
    }

    plotter->SetLegend(0.2, 0.3);
    plotter->SetLegendOffset(0.2, 0.1);

    plotter->SetYaxis("#sigma(#gamma p #rightarrow J/#psi  p)    [nb]", 0., 5.);
    plotter->SetXaxis("E_{#gamma}    [GeV]", xmin, xmax);
    plotter->Plot("helamp.pdf");

    return 1;
};