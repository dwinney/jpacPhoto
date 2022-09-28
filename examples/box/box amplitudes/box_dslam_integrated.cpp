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

    auto Wcut = [&] (double qmax)
    {
        return sqrt(qmax*qmax + M_DSTAR*M_DSTAR) + sqrt(qmax*qmax + M_LAMBDAC*M_LAMBDAC);
    };

    double eta  = 1.;

    interpolated_discontinuity disc1(&kBox, 1);
    disc1.import_data("./grid_data/boxDs_");

    box_amplitude box1(&kBox, &disc1, "#it{J}_{max} = 1/2");
    box1.set_intermediate_threshold(M_DSTAR + M_LAMBDAC + 1.E-4);

    interpolated_discontinuity disc3(&kBox, 3);
    disc3.import_data("./grid_data/boxDs_");

    box_amplitude box3(&kBox, &disc3, "#it{J}_{max} = 3/2");
    box3.set_intermediate_threshold(M_DSTAR + M_LAMBDAC + 1.E-4);

    interpolated_discontinuity disc5(&kBox, 5);
    disc5.import_data("./grid_data/boxDs_");

    box_amplitude box5(&kBox, &disc5, "#it{J}_{max} = 5/2");
    box5.set_intermediate_threshold(M_DSTAR + M_LAMBDAC + 1.E-4);

    // ---------------------------------------------------------------------------
    // Plot settings
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;
    amps.push_back(&box1);
    amps.push_back(&box3);
    amps.push_back(&box5);


    // Options
    int N = 50;
    double  xmin = 8.0;
    double  xmax = 10.5;

    double  ymin = 0.;
    double  ymax = 0.8;

    std::string filename = "sigma.pdf";
    std::string ylabel  = "#sigma(#gamma#it{p} #rightarrow #it{J}/#psi #it{p})  [nb]";
    bool PRINT = true;

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    for (int n = 0; n < amps.size(); n++)
    {
        auto F = [&](double x)
        {
            double w = W_cm(x);
            return amps[n]->integrated_xsection(w*w);
        };  
        
        amps[n]->set_params( {Wcut(1.), eta} );
        plotter->AddEntry(N, F, {xmin,xmax}, amps[n]->get_id(), PRINT);
        amps[n]->set_params( {Wcut(1.2), eta} );
        plotter->AddDashedEntry(N, F, {xmin,xmax}, PRINT);
    }

    plotter->SetXaxis("#it{E}_{#gamma}  [GeV]", xmin, xmax);

    plotter->SetYaxis(ylabel, ymin, ymax);
    
    plotter->SetLegend(0.2, 0.73);
    plotter->SetLegendOffset(0.5, 0.15);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

};