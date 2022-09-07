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


#include <chrono>
using namespace std::chrono;

using namespace jpacPhoto;

void box_dlam_integrated()
{
    reaction_kinematics kBox(M_JPSI, M_PROTON);
    kBox.set_meson_JP( {1, -1} );

    double Wcut = 4.61 ;
    double eta  = 1.;
    int    jMax = 1;

    interpolated_discontinuity disc1(&kBox, 1);
    disc1.import_data("./grid_data/boxD_");

    box_amplitude box1(&kBox, &disc1, "#it{J}_{max} = 1/2");
    box1.set_intermediate_threshold(M_D + M_LAMBDAC);
    box1.set_params( {Wcut, eta} );

    interpolated_discontinuity disc3(&kBox, 3);
    disc3.import_data("./grid_data/boxD_");

    box_amplitude box3(&kBox, &disc3, "#it{J}_{max} = 3/2");
    box3.set_intermediate_threshold(M_D + M_LAMBDAC);
    box3.set_params( {Wcut, eta} );

    // ---------------------------------------------------------------------------
    // Plot settings
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(&box1);
    amps.push_back(&box3);

    // Options
    int N =  50;
    double  xmin = 8.0;
    double  xmax = 10.5;

    double  ymin = 0.;
    double  ymax = 0.2;

    std::string filename = "sigma.pdf";
    std::string ylabel  = "#sigma(#gamma#it{p} #rightarrow #it{J}/#psi #it{p})  [nb]";
    bool PRINT = true;

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D* plotter = new jpacGraph1D();

    for (int i = 0; i < amps.size(); i++)
    {
        // ---------------------------------------------------------------------------
        // Print the desired observable for each amplitude
        auto F = [&](double x)
        {
            double w = W_cm(x);
            return amps[i]->integrated_xsection(w*w);
        };  
        
        plotter->AddEntry(N, F, {xmin,xmax}, amps[i]->get_id(), PRINT);
    }


    plotter->SetXaxis("#it{E}_{#gamma}  [GeV]", xmin, xmax);

    plotter->SetYaxis(ylabel, ymin, ymax);
    
    plotter->SetLegend(0.2, 0.73);
    plotter->SetLegendOffset(0.5, 0.1);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

};