// Script to examine the breakdown of different L components to the total 
// cross-section from the 2-channel model compared to a pomeron exchange
//
// OUTPUT: jpsi007_results.pdf, gluex_results.pdf
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#include "constants.hpp"
#include "kinematics.hpp"
#include "partial_wave.hpp"
#include "plotter.hpp"

#include "analytic/K_matrix.hpp"
#include "analytic/pomeron_exchange.hpp"
#include "gluex/data.hpp"
#include "gluex/plots.hpp"
#include "jpsi007/data.hpp"
#include "jpsi007/plots.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

void pw_breakdown()
{
    using namespace jpacPhoto;
    using pomeron_exchange = analytic::pomeron_exchange;
    using K_matrix         = analytic::K_matrix;

    // ---------------------------------------------------------------------------
    // Amplitude setup
    // ---------------------------------------------------------------------------

    // J/psi proton final state this is shared by both amplitudes
    kinematics kJpsi = new_kinematics(M_JPSI, M_PROTON);
    kJpsi->set_meson_JP( {1, -1} );

    // K-MATRIX
    std::array<double,2> higher = {M_DSTAR, M_LAMBDAC};
    amplitude s = new_amplitude<K_matrix>(kJpsi, 0, higher, "S-wave");
    s->set_option(EffectiveRange);

    amplitude p = new_amplitude<K_matrix>(kJpsi, 1, "P-wave");
    amplitude d = new_amplitude<K_matrix>(kJpsi, 2, "D-wave");
    amplitude f = new_amplitude<K_matrix>(kJpsi, 3, "F-wave");

    // Amplitude being fit is the sum of them all
    amplitude sum_kmatrix = s + p + d + f;
    sum_kmatrix->set_id("Sum up to L_{max} = 3");

    std::vector<double> pars = {-0.099814509, -3.1814832, -4.2314169, 0.098642893, 0.94256618, -3.590315, -2.9318495,
                                -0.012730907, -1.1314573,
                                -0.0038812986, 0.075979432,
                                -0.00064385586, 0.22214291 };
    sum_kmatrix->set_parameters(pars);

    // POMERON AMPLITUDE
    amplitude pomeron = new_amplitude<pomeron_exchange>(kJpsi, "Full Pomeron");
    
    std::vector<double> pars_p = {0.064647133, 0.2907296, 0.97284385, 0.19430205};
    pomeron->set_parameters(pars_p);

    amplitude s_p = project(0, false, pomeron, "S-wave");
    amplitude p_p = project(1, false, pomeron, "P-wave");
    amplitude d_p = project(2, false, pomeron, "D-wave");
    amplitude f_p = project(3, false, pomeron, "F-wave");

    amplitude sum_pomeron = s_p + p_p + d_p + f_p;
    sum_pomeron->set_id("Pomeron (L_{max} = 3)");

    std::vector<amplitude> amps = {sum_kmatrix, s, p, d, f};
    std::vector<amplitude> amps_dashed = {sum_pomeron, s_p, p_p, d_p, f_p};

    // ---------------------------------------------------------------------------
    // Plot the results
    // ---------------------------------------------------------------------------

    plotter plotter;

    // All GlueX 2022 data
    std::vector<data_set> gluex   = gluex::all();
    std::vector<data_set> jpsi007 = jpsi007::all(); // And all Jpsi-007 data

    // -----------------------------------------
    // GlueX 

    std::vector<plot> gluex_plots;

    // Grab each pre-set plot but add the theory curve with add_curve
    plot pint = gluex::plot_integrated(plotter);

    for (int i = 0; i < 5; i++)
    {
        pint.add_curve(jpacPhoto::sigma_Egam, amps[i], {8, 11.8});
        pint.add_dashed(jpacPhoto::sigma_Egam, amps_dashed[i], {8, 11.8});
    }
    gluex_plots.push_back(pint);

    // Do the same with the differential sets
    for (int i = 0; i <= 2; i++)
    {
        double Eavg = gluex[i]._avg_w;
        double Wavg = W_cm(Eavg);
        double tmin = -kJpsi->t_min(Wavg*Wavg);
        double tmax = -kJpsi->t_max(Wavg*Wavg); 

        plot dif = gluex::plot_slice(plotter, i);

        for (int i = 0; i < 5; i++)
        {
            dif.add_curve(jpacPhoto::dsigmadt_Egam, amps[i], Eavg, {tmin, tmax});
            dif.add_dashed(jpacPhoto::dsigmadt_Egam, amps_dashed[i], Eavg, {tmin, tmax});
        }
        
        gluex_plots.push_back(dif);
    };
    
    // Print to file as a 2x2 grid
    plotter.combine({2,2}, gluex_plots, "gluex_results.pdf");
    
    // -----------------------------------------
    // J/psi-007
    
    std::vector<plot> jpsi007_plots;

    for (int i = 1; i <= 12; i++)
    {
        double Eavg = jpsi007[i-1]._avg_w;
        double Wavg = W_cm(Eavg);
        double tmin = -kJpsi->t_min(Wavg*Wavg);
        double tmax = -kJpsi->t_max(Wavg*Wavg);

        plot dif =  jpsi007::plot_slice(plotter, i);

        for (int i = 0; i < 5; i++)
        {
            auto dsig_tp = [&](double tp)
            {
                double t = tp + tmin;
                return amps[i]->differential_xsection(Wavg*Wavg, -t);
            };
            dif.add_curve({0, tmax-tmin}, dsig_tp, amps[i]->id());

            auto dsig_tpd = [&](double tp)
            {
                double t = tp + tmin;
                return amps_dashed[i]->differential_xsection(Wavg*Wavg, -t);
            };
            dif.add_dashed({0, tmax-tmin}, dsig_tpd);
        }
        jpsi007_plots.push_back(dif);
    };

    plotter.combine({3,4}, jpsi007_plots, "jpsi007_results.pdf");
};