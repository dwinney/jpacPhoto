// Script for extracting dynamical quantities form the fit results of the S-wave 
// K-matrices in J/psi photo production
//
// OUTPUT: eta_th, R_VMD, a_psip, and r_psip printed in table on cmdline
//         plot of eta as a function of E saved to file (inelasticity.pdf)
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#include "constants.hpp"
#include "kinematics.hpp"
#include "plotter.hpp"
#include "analytic/K_matrix.hpp"
#include "gluex/data.hpp"
#include "gluex/plots.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

void psip_dynamics()
{
    using namespace jpacPhoto;
    using K_matrix         = analytic::K_matrix;

    // J/psi proton final
    kinematics kJpsi = new_kinematics(M_JPSI, M_PROTON);
    kJpsi->set_meson_JP( {1, -1} );
    
    // Threshold beam energy
    double wth = kJpsi->Wth();
    double eth = E_beam(wth);

    // Open charm channels
    std::array<double,2> lower  = {M_D,     M_LAMBDAC};
    std::array<double,2> higher = {M_DSTAR, M_LAMBDAC};
    std::array<std::array<double,2>,2> both = {lower, higher};

    // ---------------------------------------------------------------------------
    // One-Channel (Elastic J/psi p)
    // ---------------------------------------------------------------------------

    std::vector<double> pars_1c = {-0.063170013, -8.3206195, 6.3813772};

    double a1_1c = -2.6613031;

    // Create it as a pointer to generic amplitude
    amplitude x1c = new_amplitude<K_matrix>(kJpsi, 0, "1C Elastic");
    x1c->set_option(EffectiveRange);
    x1c->set_parameters(pars_1c);
    
    // ---------------------------------------------------------------------------
    // Two-Channel (J/psi p and D* Lambdac)
    // ---------------------------------------------------------------------------

    std::vector<double> pars_2c = { -0.099814509, -3.1814832,             // n's
                                    -4.2314169, 0.098642893, 0.94256618,  // a's
                                    -3.590315, -2.9318495 };              // b's

    double a1_2c = -0.87250104;

    // Create it as a pointer to generic amplitude
    amplitude x2c = new_amplitude<K_matrix>(kJpsi, 0, higher, "2C Elastic");
    x2c->set_option(EffectiveRange);
    x2c->set_parameters(pars_2c);

    // ---------------------------------------------------------------------------
    // Three-Channel (J/psi p and both open chamrs)
    // ---------------------------------------------------------------------------

    std::vector<double> pars_3c = { -0.10284626, 0.10671937, 0.099786669,
                                    -4.7926673, 3.1841738, -2.4125145,
                                    -3.0198596, 4.2151254, 1.5583444 };

    double a1_3c = -1.219443;

    // Create it as a pointer to generic amplitude
    amplitude x3c = new_amplitude<K_matrix>(kJpsi, 0, both, "3C Elastic");
    x3c->set_parameters(pars_3c);

    // ---------------------------------------------------------------------------
    // Results
    // ---------------------------------------------------------------------------

    std::array<amplitude,3> amps = {x1c, x2c, x3c};
    std::array<double,3> a_pwave = {a1_1c, a1_2c, a1_3c};

    plotter plotter;
    plot p1 = plotter.new_plot();
    p1.set_curve_points(1000);

    line(); 
    divider(6);
    print("", "eta_th", "R_VMD",  "a_psip [am]",  "r_psip  [am]", "(q r_psip)^2");
    divider(6);

    for (int i = 0; i < 3; i++)
    {
        // Downcast to K_matrix to access non-virtual functions]
        std::string si = std::to_string(i+1);
        auto amp = std::dynamic_pointer_cast<K_matrix>( amps[i] );
        std::string label = si + "-channel";

        // Add curves of elasticity to the plot
        auto eta = [&](double e)
        {
            return amp->inelasticity(e);
        };
        p1.add_curve( {8.2, 11.8}, eta, label+ " ("+si+ "C)");

        // Calcuate  the range parameter
        double r2 = std::abs(a_pwave[i] / (8*PI*wth*amp->scattering_length()*5.068E-3));

        // Print in a little table
        print(label, amp->inelasticity(eth), amp->VMD_test(), amp->scattering_length(), sqrt(r2)/5.068E-3, pow(kJpsi->final_momentum(W_cm(12)*W_cm(12)),2)*r2);
    };
    divider(6); line();

    p1.set_ranges( {8, 12}, {0, 1});
    p1.set_legend( 0.6, 0.65);
    p1.set_labels( "#it{E}_{#gamma}  [GeV]", "1 #minus #eta" );
    p1.save("inelasticity.pdf");
};