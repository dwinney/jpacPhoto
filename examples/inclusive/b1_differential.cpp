#include "triple_regge.hpp"
#include "inclusive_kinematics.hpp"
#include "regge_trajectory.hpp"

#include "jpacUtils.hpp"
#include "jpacGraph1D.hpp"

#include <iostream>
#include <cstring>

using namespace jpacPhoto;

void b1_differential()
{
    // ---------------------------------------------------------------------------
    // Data from OmegaPHOTON 
    // ---------------------------------------------------------------------------
    
    double s = 75.9421;
    std::vector<double> dat_x = {0.65, 0.75, 0.85, 0.95};
    std::vector<double> err_x = {0.05, 0.05, 0.05, 0.05};

    std::vector<double> dat_sigma = {1.80957, 2.15690, 1.3661, 0.65901};
    std::vector<double> err_sigma = {2.36188 - 1.80957, 2.47490 - 2.15690, 1.53345 - 1.36611, 0.76779 - 0.65901};

    // ---------------------------------------------------------------------------
    // Amplitudes
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // Global quantities

    // Couplings
    double g_b1 = 0.24;
    double g_NN = sqrt(2.) * sqrt(4. * PI * 13.81); // Nucleon coupling same for all
    double LamPi = .9;  // 900 MeV cutoff for formfactor
    
    // Kinematics
    reaction_kinematics kb1 (M_B1);
    kb1.set_meson_JP(1, +1);

    // ---------------------------------------------------------------------------
    // Reggeized pion amplitude 

    // Pion trajectory 
    int signature = +1;
    double alpha_prime = 0.7; // GeV^-2
    double alpha_0 =  - alpha_prime * M2_PION;
    linear_trajectory alpha (signature, alpha_0, alpha_prime);
    alpha.set_minJ(0);

    // Exclusive amplitude
    pseudoscalar_exchange excB1r (&kb1, &alpha, "b1 production");
    excB1r.set_params({g_b1, g_NN});
    excB1r.set_formfactor(true, LamPi);

    // We now can pass this to an inclusive amplitude
    triple_regge incB1r (&excB1r);
    incB1r.set_high_energy_approximation(true);

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    int N = 1000;

    double xmin = 0.7;   
    double xmax = 1.;

    double ymin = 0.;
    double ymax = 2.5;

    std::string filename = "dsigmadx.pdf";
    std::string ylabel   = "d#sigma/d#it{x}   [#mub]";
    std::string xlabel   = "#it{x}";

    // ---------------------------------------------------------------------------  
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    //Add entries for all the different differential cross-sections at fixed s

    auto F = [&](double x)
    {
        return incB1r.dsigma_dx(s, x) * 1.E-3; // in mub!
    };

    incB1r.set_sigma_total(JPAC_pimp_withResonances);
    plotter->AddEntry(N, F, {xmin, xmax}, "Inclusive #it{b}_{1}(1235)^{#plus}");

    // Remove resonances
    incB1r.set_sigma_total(PDG_pimp_onlyRegge);
    plotter->AddDashedEntry(N, F, {xmin, xmax});

    // ---------------------------------------------------------------------------
    // Add the data points 
    plotter->AddDataPoints(dat_x, dat_sigma, err_x, err_sigma, "Omega Photon");

    // ---------------------------------------------------------------------------
    // Finally make the plot pretty

    // Axes options
    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);

    // LEgend options
    plotter->SetLegend(0.25, 0.25);
    plotter->SetLegendOffset(0.3, 0.1);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
};