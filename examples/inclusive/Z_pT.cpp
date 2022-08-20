#include "triple_regge.hpp"
#include "inclusive_kinematics.hpp"
#include "regge_trajectory.hpp"

#include "jpacUtils.hpp"
#include "jpacGraph1D.hpp"

#include <iostream>
#include <cstring>

using namespace jpacPhoto;
void Z_pT()
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

    // ---------------------------------------------------------------------------
    // Zc(3900)

    // Kinematics
    reaction_kinematics kZc (M_ZC);
    kZc.set_meson_JP(1, +1);

    // Coupling
    double gc_Psi = 1.91; // psi coupling before VMD scaling
    double gc_Gamma = E * F_JPSI * gc_Psi / M_JPSI;

    // Exclusive amplitude
    pseudoscalar_exchange excZc (&kZc, M_PION, "inclusive Z_{c}");
    excZc.set_params({gc_Gamma, g_NN});
    excZc.set_formfactor(true, LamPi);

    // We now can pass this to an inclusive amplitude
    triple_regge incZc (&excZc);
    incZc.set_high_energy_approximation(false);
    incZc.set_sigma_total(JPAC_pipp_withResonances);

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    int N = 100;
    double w = 6.;

    incZc._kinematics->_s = w*w;
    double pmax = incZc._kinematics->pMax();

    double xmin = 0.01;   
    double xmax = 1.2;

    double ymin = 2.E-5;
    double ymax = 50;

    std::string filename = "incZ_qT.pdf";
    std::string ylabel   = "d#sigma/d#it{q}_{T}     [nb / GeV]";
    std::string xlabel   = "#it{q}_{T}    [GeV]";
    bool print_to_CMD    = true;

    // ---------------------------------------------------------------------------  
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();
    
    // ---------------------------------------------------------------------------
    // Intermediate functions needed

    // Just the inclusive curve
    auto F = [&](double pT)
    {
        double y = pT / pmax;
        return incZc.dsigma_dy(w*w, y) * 1.E3; // in nb!
    };

    incZc.set_sigma_total(JPAC_pipp_withResonances);
    plotter->AddEntry(N, F, {xmin, xmax}, "#it{Z}_{c}(3900)^{#minus}", print_to_CMD);

    incZc.set_sigma_total(JPAC_pimp_withResonances);
    plotter->AddEntry(N, F, {xmin, xmax}, "#it{Z}_{c}(3900)^{#plus}", print_to_CMD);

    // ---------------------------------------------------------------------------
    // Finally make the plot pretty

    // Axes options
    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(true);

    std::ostringstream streamObj;
    streamObj << std::setprecision(4) << "#it{W}_{#gammap} = " << w << " GeV";
    std::string header = streamObj.str();

    // LEgend options
    plotter->SetLegend(0.35, 0.35, header);
    plotter->SetLegendOffset(0.3, 0.12);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
};