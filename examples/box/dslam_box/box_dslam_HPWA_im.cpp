// ---------------------------------------------------------------------------
// Read in 2D interpolations of the PWA for b^J and c^J and assemble the 
// imaginary of the partial-wave of the full box amplitude
// 
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "reaction_kinematics.hpp"
#include "interpolated_discontinuity.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

void box_dslam_HPWA_im()
{
    int h;
    int J = 1;
    int N = 200;

    std::string filename = "boxDs_hpwa_im.pdf";

    auto Wcut = [&] (double qmax)
    {
        return sqrt(qmax*qmax + M_D*M_D) + sqrt(qmax*qmax + M_LAMBDAC*M_LAMBDAC);
    };

    double qmax = 1.;
    double eta  = 1.;
    int    jMax = 1;

    // Need the kinematics of the intermediate reactions to get phases and helicity combinations
    reaction_kinematics kBox (M_JPSI, M_PROTON);
    kBox.set_meson_JP(1, -1);

    interpolated_discontinuity disc(&kBox, jMax);
    disc.import_data("./grid_data/boxDs_");
    disc.set_intermediate_threshold(M_DSTAR + M_LAMBDAC + 1.E-3);
    disc.set_params( { Wcut(1.), eta } );

    // Plotter object
    jpacGraph1D* plotter = new jpacGraph1D();

    double smin = kBox.Wth();
    double smax = Wcut(1.);

    double ymin = -0.07;
    double ymax =  0.1;

    bool verbose = true;
    bool PRINT = true;

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    auto F = [&](double w)
    {
        return std::imag( disc.helicity_pwa(h, J, w*w) );
    };
    
    h = 7;
    plotter->AddEntry(N, F, {smin,smax}, "#{}{ + +, + + }", PRINT);

    h = 9;
    plotter->AddEntry(N, F, {smin,smax}, "#{}{ + +, 0 + }", PRINT);
    
    h = 8;
    plotter->AddEntry(N, F, {smin,smax}, "#{}{ + +, 0 #minus }", PRINT);

    h = 10;
    plotter->AddEntry(N, F, {smin,smax}, "#{}{ + +, #minus #minus }", PRINT);

    plotter->SetXaxis("#it{W}   [GeV]", smin, smax);
    plotter->SetYaxis("Im #it{a}_{#{}{#lambda}}^{1/2}(s)", ymin, ymax);
    plotter->SetLegend(0.2, 0.73);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
};