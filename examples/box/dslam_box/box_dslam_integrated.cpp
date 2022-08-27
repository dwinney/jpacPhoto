// ---------------------------------------------------------------------------
// Make sure that interpolated_discontinuity is importing amplitudes correctly
// 
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "reaction_kinematics.hpp"
#include "interpolated_discontinuity.hpp"
#include "box_amplitude.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

void box_dslam_integrated()
{
    reaction_kinematics kBox(M_JPSI, M_PROTON);
    kBox.set_meson_JP( {1, -1} );

    double Wcut = 4.61 ;
    double eta  = 1.;
    int    jMax = 1;

    interpolated_discontinuity disc(&kBox, jMax);
    disc.import_data("./grid_data/boxDs_");

    box_amplitude box(&kBox, &disc, "#bar{D}* #Lambda_{c} box");
    box.set_intermediate_threshold(M_DSTAR + M_LAMBDAC + 1.E-3);
    box.set_params( {Wcut, eta} );

    // ---------------------------------------------------------------------------
    // Plot settings
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(&box);

    // Options
    int N =  100;
    double  xmin = 8.0;
    double  xmax = 10.5;

    double  ymin = 0.;
    double  ymax = 1.427;

    std::string filename = "box.pdf";
    std::string ylabel  = "#sigma(#gamma#it{p} #rightarrow #it{J}/#psi #it{p})  [nb]";
    bool PRINT = true;

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    auto F = [&](double x)
    {
        double w = W_cm(x);
        return box.integrated_xsection(w*w);
    };  
    
    plotter->AddEntry(N, F, {xmin,xmax}, "#bar{D} #Lambda_{c}^{+} box", PRINT);

    plotter->SetXaxis("#it{E}_{#gamma}  [GeV]", xmin, xmax);

    plotter->SetYaxis(ylabel, ymin, ymax);
    
    plotter->SetLegend(0.2, 0.73);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

};