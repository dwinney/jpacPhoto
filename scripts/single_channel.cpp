
#include "constants.hpp"
#include "kinematics.hpp"
#include "fitter.hpp"
#include "plotter.hpp"

#include "analytic/scattering_length.hpp"
#include "gluex/data.hpp"
#include "gluex/plots.hpp"
#include "jpsi007/data.hpp"
#include "jpsi007/plots.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

void single_channel()
{
    using namespace jpacPhoto;
    using scattering_length = analytic::scattering_length;

    // ---------------------------------------------------------------------------
    // Amplitude setup
    // ---------------------------------------------------------------------------

    // J/psi proton final
    kinematics kJpsi = new_kinematics(M_JPSI, M_PROTON);
    kJpsi->set_meson_JP( {1, -1} );

    // Any number of individual partial waves
    amplitude s = new_amplitude<scattering_length>(kJpsi, 0, "S-wave");
    amplitude p = new_amplitude<scattering_length>(kJpsi, 1, "P-wave");
    amplitude d = new_amplitude<scattering_length>(kJpsi, 2, "D-wave");
    amplitude f = new_amplitude<scattering_length>(kJpsi, 3, "F-wave");
    amplitude g = new_amplitude<scattering_length>(kJpsi, 4, "G-wave");

    // Amplitude being fit is the sum of them all
    amplitude to_fit = s + p + d /* + f  + g */;
    to_fit->set_id("Sum");

    // ---------------------------------------------------------------------------
    //  Fitter setup
    // ---------------------------------------------------------------------------

    fitter fitter(to_fit, "Migrad", 1.E-4);

    // -----------------------------------------
    // Choose which data we want to fit against

    // All GlueX 2022 data
    std::vector<data_set> gluex   = gluex::all();
    std::vector<data_set> jpsi007 = jpsi007::all(); // And all Jpsi-007 data

    // Load everything into the fitter
    fitter.add_data(gluex);
    fitter.add_data(jpsi007);
    
    // -----------------------------------------
    // Set up out initial guess

    // Fit N times with randomly sampled inital parameters
    fitter.set_parameter_limits("b[0]", {0., 5.});  // Force the first normalization to be positive
    fitter.do_fit(50);

    // ---------------------------------------------------------------------------
    // Plot the results
    // ---------------------------------------------------------------------------

    plotter plotter;

    // -----------------------------------------
    // GlueX 

    std::vector<plot> gluex_plots;

    // Grab each pre-set plot but add the theory curve with add_curve
    plot pint = gluex::plot_integrated(plotter);
    pint.add_curve(jpacPhoto::sigma_Egam, to_fit, {8., 11.8});
    gluex_plots.push_back(pint);

    // Do the same with the differential sets
    for (int i = 0; i <= 2; i++)
    {
        double Eavg = gluex[i]._avg_w;
        double Wavg = W_cm(Eavg);
        double tmin = -kJpsi->t_min(Wavg*Wavg);
        double tmax = -kJpsi->t_max(Wavg*Wavg); 

        plot dif = gluex::plot_slice(plotter, i);
        dif.add_curve(jpacPhoto::dsigmadt_Egam, to_fit, Eavg, {tmin, tmax});

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

        auto dsig_tp = [&](double tp)
        {
            double t = tp + tmin;
            return to_fit->differential_xsection(Wavg*Wavg, -t);
        };

        plot dif =  jpsi007::plot_slice(plotter, i);
        dif.add_curve({0, tmax-tmin}, dsig_tp, to_fit->id());
        jpsi007_plots.push_back(dif);
    };

    plotter.combine({4,3}, jpsi007_plots, "jpsi007_results.pdf");
};