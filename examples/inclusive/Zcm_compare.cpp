#include "triple_regge.hpp"
#include "inclusive_kinematics.hpp"
#include "regge_trajectory.hpp"

#include "jpacUtils.hpp"
#include "jpacGraph1D.hpp"

#include <iostream>
#include <cstring>

using namespace jpacPhoto;

void Zcm_compare()
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
    pseudoscalar_exchange excZc (&kZc, M_PION, "total #it{Z}_{#it{c}}(3900)^{#minus} production");
    excZc.set_params({gc_Gamma, g_NN});
    excZc.set_formfactor(true, LamPi);

    // We now can pass this to an inclusive amplitude
    triple_regge incZc (&excZc);
    incZc.set_high_energy_approximation(false);

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    int N = 100;

    double xmin = 4.8;   
    double xmax = 10.;

    double ymin = 0;
    double ymax = 45.;

    std::string filename = "incZminus_compare.pdf";
    std::string ylabel   = "#sigma  [nb]";
    std::string xlabel   = "#it{W}_{#gamma#it{p}}  [GeV]";
    bool print_to_CMD    = true;

    // ---------------------------------------------------------------------------  
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------
 
    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();

    auto F = [&](double w)
    {
        return incZc.integrated_xsection(w*w) * 1.E3; // in nb!
    };

    auto G = [&](double w)
    {
        return incZc.exclusive_xsection(w*w); // in nb!
    };

    auto H = [&](double w)
    {
        return G(w) + F(w); // in nb!
    };

    // Exlclusive only
    plotter->AddEntry(N, G, {xmin, xmax}, "#it{Z}_{c}(3900)^{#minus} #it{n}", print_to_CMD);

    incZc.set_sigma_total(JPAC_pimp_onlyDelta);
    plotter->AddEntry(N, F, {xmin, xmax}, "#it{Z}_{c}(3900)^{#minus} (#Delta^{0} #rightarrow #pi^{#minus}#it{p})", print_to_CMD);

    incZc.set_sigma_total(JPAC_pimp_withResonances);
    plotter->AddEntry(N, H, {xmin, xmax}, "Total", print_to_CMD);

    incZc.set_sigma_total(JPAC_pimp_onlyDelta);
    plotter->AddDashedEntry(N, H, {xmin, xmax}, print_to_CMD);
    

    // ---------------------------------------------------------------------------
    // Finally make the plot pretty

    // Axes options
    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    // plotter->SetYlogscale(true);

    // LEgend options
    plotter->SetLegend(0.23, 0.75);
    plotter->SetLegendOffset(0.3, 0.15);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
};