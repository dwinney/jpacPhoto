#include "triple_regge.hpp"
#include "inclusive_kinematics.hpp"
#include "regge_trajectory.hpp"

#include "jpacUtils.hpp"
#include "jpacGraph1D.hpp"

#include <iostream>
#include <cstring>

using namespace jpacPhoto;

void Zp_totals()
{
    // ---------------------------------------------------------------------------
    // Amplitudes
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // Global quantities

    // Couplings
    double g_NN = sqrt(2.) * sqrt(4. * PI * 13.81); // Nucleon coupling same for all
    double LamPi = .9;  // 900 MeV cutoff for formfactor

    // Masses
    double M_ZC  = 3.8884; // GeV
    double M_ZB  = 10.6072;
    double M_ZBP = 10.6522;

    // ---------------------------------------------------------------------------
    // Zc(3900)

    // Kinematics
    reaction_kinematics kZc (M_ZC);
    kZc.set_meson_JP(1, 1);

    // Coupling
    double gc_Psi = 1.91; // psi coupling before VMD scaling
    double gc_Gamma = E * F_JPSI * gc_Psi / M_JPSI;

    // Exclusive amplitude
    pseudoscalar_exchange excZc (&kZc, M_PION, "total #it{Z}_{#it{c}}(3900)^{#plus} production");
    excZc.set_params({gc_Gamma, g_NN});
    excZc.set_formfactor(true, LamPi);

    // We now can pass this to an inclusive amplitude
    triple_regge incZc (&excZc);
    incZc.set_high_energy_approximation(false);
    incZc.set_sigma_total(JPAC_pipp_withResonances);
    incZc.set_sigma_total_option(5);

    // ---------------------------------------------------------------------------
    // Zb(10610)

    // Kinematics
    reaction_kinematics kZb (M_ZB);
    kZb.set_meson_JP(1, 1);

    // Coupling 
    double gb_Ups1 = 0.49, gb_Ups2 = 3.30, gb_Ups3 = 9.22;
    double gb_Gamma = E * ( F_UPSILON1S * gb_Ups1 / M_UPSILON1S 
                          + F_UPSILON2S * gb_Ups2 / M_UPSILON2S
                          + F_UPSILON3S * gb_Ups3 / M_UPSILON3S);  

    // Exclusive amplitude
    pseudoscalar_exchange excZb (&kZb, M_PION, "total #it{Z}_{#it{b}}(10610)^{#plus} production");
    excZb.set_params({gb_Gamma, g_NN});
    excZb.set_formfactor(true, LamPi);

    // We now can pass this to an inclusive amplitude
    triple_regge incZb (&excZb);
    incZb.set_high_energy_approximation(false);
    incZb.set_sigma_total(JPAC_pipp_withResonances);
    
    // ---------------------------------------------------------------------------
    // Zb(10650)

    // Kinematics
    reaction_kinematics kZbp (M_ZBP);
    kZbp.set_meson_JP(1, 1);

    // Coupling 
    double gbp_Ups1 = 0.21, gbp_Ups2 = 1.47, gbp_Ups3 = 4.8;
    double gbp_Gamma = E * ( F_UPSILON1S * gbp_Ups1 / M_UPSILON1S 
                           + F_UPSILON2S * gbp_Ups2 / M_UPSILON2S
                           + F_UPSILON3S * gbp_Ups3 / M_UPSILON3S);  

    // Exclusive amplitude
    pseudoscalar_exchange excZbp (&kZbp, M_PION, "total #it{Z}_{#it{b}}(10650)^{#plus} production");
    excZbp.set_params({gbp_Gamma, g_NN});
    excZbp.set_formfactor(true, LamPi);

    // We now can pass this to an inclusive amplitude
    triple_regge incZbp (&excZbp);
    incZbp.set_high_energy_approximation(false);
    incZbp.set_sigma_total(JPAC_pipp_withResonances);

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    int N = 100;

    double xmin = 4.8;   
    double xmax = 20.;

    double ymin = 2.E-1;
    double ymax = 1.E2;

    std::string filename = "incZplus_total.pdf";
    std::string ylabel   = "#sigma  [nb]";
    std::string xlabel   = "#it{W}_{#gamma#it{p}}  [GeV]";
    bool print_to_CMD    = true;

    // ---------------------------------------------------------------------------  
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------
 
    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();

    int amp = 0;
    std::vector<triple_regge*> inc = {&incZc, &incZb, &incZbp};

    auto F = [&](double w)
    {
        return inc[amp]->integrated_xsection(w*w) * 1.E3 + inc[amp]->exclusive_xsection(w*w); // in nb!
    };

    auto G = [&](double w)
    {
        return inc[amp]->exclusive_xsection(w*w); // in nb!
    };

    for (int i = 0; i < 3; i++)
    {
        amp = i;
        plotter->AddEntry(N, F, {xmin, xmax}, inc[i]->get_id(), print_to_CMD);
        plotter->AddDashedEntry(N, G, {xmin, xmax}, print_to_CMD);
    }

    // ---------------------------------------------------------------------------
    // Finally make the plot pretty

    // Axes options
    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(true);

    // LEgend options
    plotter->SetLegend(0.53, 0.69);
    plotter->SetLegendOffset(0.3, 0.15);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
};