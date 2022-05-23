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
    double g_NN = sqrt(2.) * sqrt(4. * PI * 13.81); // Nucleon coupling same for all
    double LamPi = .9;  // 900 MeV cutoff for formfactor
    double bPi = 1. / (LamPi * LamPi);

    // Masses
    double M_ZC  = 3.8884; // GeV
    double M_ZB  = 10.6072;
    double M_ZBP = 10.6522;

    // ---------------------------------------------------------------------------
    // Zc(3900)

    // Kinematics
    reaction_kinematics * kZc = new reaction_kinematics(M_ZC);
    kZc->set_JP(1, 1);

    // Coupling
    double gc_Psi = 1.91; // psi coupling before VMD scaling
    double gc_Gamma = E * F_JPSI * gc_Psi / M_JPSI;

    // Exclusive amplitude
    std::unique_ptr<pseudoscalar_exchange> excZc( new pseudoscalar_exchange(kZc, M_PION, "total Z_{c}(3900)^{-} production") );
    excZc->set_params({gc_Gamma, g_NN});
    excZc->set_formfactor(true, bPi);

    // We now can pass this to an inclusive amplitude
    std::unique_ptr<triple_regge> incZc( new triple_regge(excZc.get()));
    incZc->set_high_energy_approximation(true);

    // ---------------------------------------------------------------------------
    // Zb(10610)

    // Kinematics
    reaction_kinematics * kZb = new reaction_kinematics(M_ZB);
    kZb->set_JP(1, 1);

    // Coupling 
    double gb_Ups1 = 0.49, gb_Ups2 = 3.30, gb_Ups3 = 9.22;
    double gb_Gamma = E * ( F_UPSILON1S * gb_Ups1 / M_UPSILON1S 
                          + F_UPSILON2S * gb_Ups2 / M_UPSILON2S
                          + F_UPSILON3S * gb_Ups3 / M_UPSILON3S);  

    // Exclusive amplitude
    std::unique_ptr<pseudoscalar_exchange> excZb( new pseudoscalar_exchange(kZb, M_PION, "total Z_{b}(10610)^{-} production") );
    excZb->set_params({gb_Gamma, g_NN});
    excZb->set_formfactor(true, bPi);

    // We now can pass this to an inclusive amplitude
    std::unique_ptr<triple_regge> incZb( new triple_regge(excZb.get()));
    incZb->set_high_energy_approximation(true);
    
    // ---------------------------------------------------------------------------
    // Zb(10650)

    // Kinematics
    reaction_kinematics * kZbp = new reaction_kinematics(M_ZBP);
    kZbp->set_JP(1, 1);

    // Coupling 
    double gbp_Ups1 = 0.21, gbp_Ups2 = 1.47, gbp_Ups3 = 4.8;
    double gbp_Gamma = E * ( F_UPSILON1S * gbp_Ups1 / M_UPSILON1S 
                           + F_UPSILON2S * gbp_Ups2 / M_UPSILON2S
                           + F_UPSILON3S * gbp_Ups3 / M_UPSILON3S);  

    // Exclusive amplitude
    std::unique_ptr<pseudoscalar_exchange> excZbp( new pseudoscalar_exchange(kZbp, M_PION, "total Z_{b}(10650)^{-} production") );
    excZbp->set_params({gbp_Gamma, g_NN});
    excZbp->set_formfactor(true, bPi);

    // We now can pass this to an inclusive amplitude
    std::unique_ptr<triple_regge> incZbp( new triple_regge(excZbp.get()));
    incZbp->set_high_energy_approximation(true);


    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    int N = 100;

    double xmin = 4.8;   
    double xmax = 25.;

    double ymin = 1.E-2;
    double ymax = 70.;

    std::string filename = "incZminus_total.pdf";
    std::string ylabel   = "#sigma  [nb]";
    std::string xlabel   = "W  [GeV]";
    bool print_to_CMD    = true;

    // ---------------------------------------------------------------------------  
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------
    
    // Plotter object
    std::unique_ptr<jpacGraph1D> plotter( new jpacGraph1D() );

    int amp = 0;
    std::vector<triple_regge*> inc = {incZc.get(), incZb.get(), incZbp.get()};

    auto F = [&](double w)
    {
        return inc[amp]->integrated_xsection(w*w) * 1.E3; // in nb!
    };

    for (int i = 0; i < 3; i++)
    {
        amp = i;
        inc[amp]->set_sigma_total(JPAC_pipp_withResonances);
        plotter->AddEntry(N, F, {xmin, xmax}, inc[i]->get_id(), print_to_CMD);
    }

    // ---------------------------------------------------------------------------
    // Finally make the plot pretty

    // Axes options
    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(true);

    // LEgend options
    plotter->SetLegend(0.22, 0.7);
    plotter->SetLegendOffset(0.3, 0.15);

    // Output to file
    plotter->Plot(filename);

    return 0;
};