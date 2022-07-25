#include "inclusive/triple_regge.hpp"
#include "inclusive/inclusive_kinematics.hpp"
#include "regge_trajectory.hpp"

#include "jpacUtils.hpp"
#include "jpacGraph1D.hpp"

#include <iostream>
#include <cstring>

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
    double g_delta = 18.5; // Dleta p pi coupling
    double LamPi = .9;  // 900 MeV cutoff for formfactor

    // Masses
    double M_ZC  = 3.8884; // GeV
    double M_ZB  = 10.6072;
    double M_ZBP = 10.6522;

    // Pion trajectory 
    int signature = +1;
    double alpha_prime = 0.7; // GeV^-2
    double alpha_0 =  - alpha_prime * M2_PION;
    linear_trajectory * alpha = new linear_trajectory(signature, alpha_0, alpha_prime);

    // ---------------------------------------------------------------------------
    // B1

    // Kinematics
    reaction_kinematics * kb1 = new reaction_kinematics(M_B1);
    kb1->set_meson_JP(1, +1);

    // Couplings
    double g_b1 = 0.24;

    // Exclusive amplitude
    pseudoscalar_exchange * excB1f = new pseudoscalar_exchange(kb1, alpha, "b1 production");
    excB1f->set_params({g_b1, g_NN});
    excB1f->set_formfactor(true, LamPi);

    // We now can pass this to an inclusive amplitude
    triple_regge * incB1 = new triple_regge(excB1f);
    incB1->set_high_energy_approximation(true);
    incB1->set_sigma_total(JPAC_pipp_withResonances);
  
    // ---------------------------------------------------------------------------
    // Zc(3900)

    // Kinematics
    reaction_kinematics * kZc = new reaction_kinematics(M_ZC);
    kZc->set_meson_JP(1, +1);

    // Coupling
    double gc_Psi = 1.91; // psi coupling before VMD scaling
    double gc_Gamma = E * F_JPSI * gc_Psi / M_JPSI;

    // Exclusive amplitude
    pseudoscalar_exchange * excZc = new pseudoscalar_exchange(kZc, alpha, "Zc(3900)^{#minus}  production");
    excZc->set_params({gc_Gamma, g_NN});
    excZc->set_formfactor(true, LamPi);

    // We now can pass this to an inclusive amplitude
    triple_regge * incZc = new triple_regge(excZc);
    incZc->set_high_energy_approximation(true);
    incZc->set_sigma_total(JPAC_pipp_withResonances);

    // ---------------------------------------------------------------------------
    // Zb(10610)

    // Kinematics
    reaction_kinematics * kZb = new reaction_kinematics(M_ZB);
    kZb->set_meson_JP(1, +1);

    // Coupling 
    double gb_Ups1 = 0.49, gb_Ups2 = 3.30, gb_Ups3 = 9.22;
    double gb_Gamma = E * ( F_UPSILON1S * gb_Ups1 / M_UPSILON1S 
                          + F_UPSILON2S * gb_Ups2 / M_UPSILON2S
                          + F_UPSILON3S * gb_Ups3 / M_UPSILON3S);  

    // Exclusive amplitude
    pseudoscalar_exchange * excZb = new pseudoscalar_exchange(kZb, alpha, "Zb(10610)^{#minus}  production");
    excZb->set_params({gb_Gamma, g_NN});
    excZb->set_formfactor(true, LamPi);

    // We now can pass this to an inclusive amplitude
    triple_regge * incZb = new triple_regge(excZb);
    incZb->set_high_energy_approximation(true);
    incZb->set_sigma_total(JPAC_pipp_withResonances);
    
    // ---------------------------------------------------------------------------
    // Zb(10650)

    // Kinematics
    reaction_kinematics * kZbp = new reaction_kinematics(M_ZBP);
    kZbp->set_meson_JP(1, +1);

    // Coupling 
    double gbp_Ups1 = 0.21, gbp_Ups2 = 1.47, gbp_Ups3 = 4.8;
    double gbp_Gamma = E * ( F_UPSILON1S * gbp_Ups1 / M_UPSILON1S 
                           + F_UPSILON2S * gbp_Ups2 / M_UPSILON2S
                           + F_UPSILON3S * gbp_Ups3 / M_UPSILON3S);  

    // Exclusive amplitude
    pseudoscalar_exchange * excZbp = new pseudoscalar_exchange(kZbp, alpha, "Zb(10650)^{#minus}  production");
    excZbp->set_params({gbp_Gamma, g_NN});
    excZbp->set_formfactor(true, LamPi);

    // We now can pass this to an inclusive amplitude
    triple_regge * incZbp = new triple_regge(excZbp);
    incZbp->set_high_energy_approximation(true);
    incZbp->set_sigma_total(JPAC_pipp_withResonances);

    // ---------------------------------------------------------------------------
    // 
    std::vector<triple_regge*> inc = {incB1, incZc, incZb, incZbp};

    for (int i = 0; i < inc.size(); i++)
    {
        std::cout << inc[i]->get_id() << std::endl;
        double e30 = inc[i]->exclusive_xsection(30.*30.) * 1.E3;
        double e60 = inc[i]->exclusive_xsection(60.*60.) * 1.E3;
        double e90 = inc[i]->exclusive_xsection(90.*90.) * 1.E3;
        debug(e30, e60, e90);

        double i30 = inc[i]->integrated_xsection(30.*30.) * 1.E6;
        double i60 = inc[i]->integrated_xsection(60.*60.) * 1.E6;
        double i90 = inc[i]->integrated_xsection(90.*90.) * 1.E6;
        debug(i30 + e30, i60 + e60, i90 + e90);

        std::cout << std::endl;
    };

    return 0;
};