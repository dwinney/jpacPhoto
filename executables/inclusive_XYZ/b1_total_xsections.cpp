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
    double g_NN = sqrt(2.) * sqrt(4. * PI * 13.81); // Nucleon coupling same for all
    double LamPi = .9;  // 900 MeV cutoff for formfactor
    double bPi = 1. / (LamPi * LamPi);
    
    // Kinematics
    reaction_kinematics * kb1 = new reaction_kinematics(M_B1);
    kb1->set_JP(1, +1);

    // ---------------------------------------------------------------------------
    // Fixed-spin pion amplitude 

    // Exclusive amplitude
    std::unique_ptr<pseudoscalar_exchange> excB1f( new pseudoscalar_exchange(kb1, M_PION, "b1 production") );
    excB1f->set_params({g_b1, g_NN});
    excB1f->set_formfactor(true, bPi);

    // We now can pass this to an inclusive amplitude
    std::unique_ptr<triple_regge> incB1f( new triple_regge(excB1f.get()));
    incB1f->set_high_energy_approximation(true);
    
    // // ---------------------------------------------------------------------------
    // // Plotting options
    // // ---------------------------------------------------------------------------

    int N = 400;

    double xmin = 2.;   
    double xmax = 12. ;

    double ymin = 0.;
    double ymax = 5.3;

    std::string filename = "integrated.pdf";
    std::string ylabel   = "#sigma  [#mub]";
    std::string xlabel   = "W  [GeV]";

    // ---------------------------------------------------------------------------  
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    std::unique_ptr<jpacGraph1D> plotter( new jpacGraph1D() );

    bool addExc, excOnly;
    auto F = [&](double w)
    {
        return !excOnly * incB1f->integrated_xsection(w*w) + addExc * excB1f->integrated_xsection(w*w) * 1.E-3; // in mub!
    };

    addExc = false, excOnly = false; 
    incB1f->set_sigma_total(JPAC_pipp_withResonances);
    plotter->AddEntry(N, F, {xmin, xmax}, "total b_{1}^{-} production");

    addExc = true,  excOnly = false; 
    incB1f->set_sigma_total(JPAC_pimp_withResonances);
    plotter->AddEntry(N, F, {xmin, xmax}, "total b_{1}^{+} production");

    addExc = true,  excOnly = true; 
    incB1f->set_sigma_total(JPAC_pimp_withResonances);
    plotter->AddEntry(N, F, {xmin, xmax}, "exclusive b_{1}^{+} production");

    // ---------------------------------------------------------------------------
    // Finally make the plot pretty

    // Axes options
    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);

    // LEgend options
    plotter->SetLegend(0.5, 0.65);
    plotter->SetLegendOffset(0.3, 0.15);

    // Output to file
    plotter->Plot(filename);

    return 0;
};