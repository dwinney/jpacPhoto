#include "inclusive/triple_regge.hpp"
#include "inclusive/inclusive_kinematics.hpp"
#include "regge_trajectory.hpp"

#include "jpacUtils.hpp"
#include "jpacGraph1D.hpp"

#include <iostream>
#include <cstring>
#include <memory>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

    // ---------------------------------------------------------------------------
    // Amplitudes
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // Global quantities

    // Couplings
    double g_b1 = 0.24;
    double g_N = sqrt(2.) * sqrt(4. * PI * 13.81); // Nucleon coupling same for all

    // Delta p pi coupling
    double g_delta = 2.16;
    double g_pi = sqrt(4. * M_PI / 137.);
    double Lambda = .450;

    double LamPi = .9;  // 900 MeV cutoff for formfactor
    
    // ---------------------------------------------------------------------------
    // Fixed-spin pion amplitudes

    // Kinematics for b1 Δ++ final state
    reaction_kinematics * kDelta = new reaction_kinematics(M_B1, M_DELTA);
    kDelta->set_meson_JP( 1, +1); 
    kDelta->set_baryon_JP(3, +1); 

    // Exclusive amplitude for delta
    std::unique_ptr<pseudoscalar_exchange> excDelta( new pseudoscalar_exchange(kDelta, M_PION, "b_{1}^{-} #Delta^{++}") );
    excDelta->set_params({g_b1, g_delta});
    excDelta->set_formfactor(1, LamPi);

    // Kinematics b1 proton
    reaction_kinematics * kN = new reaction_kinematics(M_B1, M_PROTON);
    kN->set_meson_JP( 1, +1);
    kN->set_baryon_JP(1, +1); 

    // Exclusive amplitude for nucleon
    std::unique_ptr<pseudoscalar_exchange> excN( new pseudoscalar_exchange(kN, M_PION, "b1 production") );
    excN->set_params({g_b1, g_N});
    excN->set_formfactor(1, LamPi);

    // We now can pass this to an inclusive amplitude
    std::unique_ptr<triple_regge> incB1( new triple_regge(excN.get()));
    incB1->set_high_energy_approximation(true);
    incB1->set_sigma_total(JPAC_pipp_withResonances);
    
    // // ---------------------------------------------------------------------------
    // // Plotting options
    // // ---------------------------------------------------------------------------

    int N = 50;

    double xmin = 2.4;   
    double xmax = 4.;

    double ymin = 0.;
    double ymax = 6.;

    std::string filename = "integrated.pdf";
    std::string ylabel   = "#sigma  [#mub]";
    std::string xlabel   = "W  [GeV]";

    // ---------------------------------------------------------------------------  
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    std::unique_ptr<jpacGraph1D> plotter( new jpacGraph1D() );
    
    auto F = [&](double w)
    {
        return excDelta->integrated_xsection(w*w) * 1.E-3; // in mub!
    };
    plotter->AddEntry(2*N, F, {xmin, xmax}, "#gamma p #rightarrow b_{1}^{-} #Delta^{++}", 1);

    auto G = [&](double w)
    {
        return incB1->integrated_xsection(w*w); // in mub!
    };  
    plotter->AddEntry(N, G, {xmin, xmax},   "#gamma p #rightarrow b_{1}^{+} X", 1);

    // ---------------------------------------------------------------------------
    // Finally make the plot pretty

    // Axes options
    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);

    // LEgend options
    plotter->SetLegend(0.6, 0.65);
    plotter->SetLegendOffset(0.3, 0.1);

    // Output to file
    plotter->Plot(filename);

    return 0;
};