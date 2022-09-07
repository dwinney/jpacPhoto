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

void box_dlam_differential()
{
    reaction_kinematics kBox(M_JPSI, M_PROTON);
    kBox.set_meson_JP( {1, -1} );

    double Wcut = 4.61 ;
    double eta  = 1.;

    interpolated_discontinuity disc1(&kBox, 1);
    disc1.import_data("./grid_data/boxD_");

    box_amplitude box1(&kBox, &disc1, "#it{J}_{max} = 1/2");
    box1.set_intermediate_threshold(M_D + M_LAMBDAC + 100.*EPS);
    box1.set_params( {Wcut, eta} );

    // ---------------------------------------------------------------------------
    // Plot settings
    // ---------------------------------------------------------------------------

    // Options
    int N =  100;

    double  ymin = 0.;
    double  ymax = 0.05;

    std::string filename = "box.pdf";
    std::string ylabel  = "d#sigma/d#it{t} (#gamma#it{p} #rightarrow #it{J}/#psi #it{p})  [nb / GeV^{2}]";
    bool PRINT = false;

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude

    double egam = 9.5;
    double W = W_cm(egam);
    double s = W*W;
    double  xmin = - kBox.t_man(s, 0.);
    double  xmax = - kBox.t_man(s, PI);

    auto F = [&](double mt)
    {
        return box1.differential_xsection(s, -mt);
    };  

    plotter->AddEntry(N, F, {xmin,xmax}, "#it{E}_{#gamma} = 9.5 GeV", PRINT);      
    plotter->SetXaxis("#minus#it{t}  [GeV^{2}]", xmin, xmax);

    egam = 9.0; W = W_cm(egam); s = W*W;  
    xmin = - kBox.t_man(s, 0.); xmax = - kBox.t_man(s, PI);
    plotter->AddEntry(N, F, {xmin,xmax}, "#it{E}_{#gamma} = 9.0 GeV", PRINT);   

    egam = 8.5; W = W_cm(egam); s = W*W;  
    xmin = - kBox.t_man(s, 0.); xmax = - kBox.t_man(s, PI);
    plotter->AddEntry(N, F, {xmin,xmax}, "#it{E}_{#gamma} = 8.5 GeV", PRINT);      

    plotter->SetYaxis(ylabel, ymin, ymax);
    
    plotter->SetLegend(0.35, 0.73);
    plotter->SetLegendOffset(0.5, 0.1);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

};