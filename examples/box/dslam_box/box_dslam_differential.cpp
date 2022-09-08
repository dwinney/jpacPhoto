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

void box_dslam_differential()
{
    reaction_kinematics kBox(M_JPSI, M_PROTON);
    kBox.set_meson_JP( {1, -1} );

    auto Wcut = [&] (double qmax)
    {
        return sqrt(qmax*qmax + M_DSTAR*M_DSTAR) + sqrt(qmax*qmax + M_LAMBDAC*M_LAMBDAC);
    };

    double eta  = 1.;
    interpolated_discontinuity disc3(&kBox, 3);
    disc3.import_data("./grid_data/boxDs_");

    box_amplitude box3(&kBox, &disc3, "#it{J}_{max} = 3/2");
    box3.set_intermediate_threshold(M_DSTAR + M_LAMBDAC + 1.E-4);
    box3.set_params( {Wcut(1.), eta} );

    // ---------------------------------------------------------------------------
    // Plot settings
    // ---------------------------------------------------------------------------

    // Options
    int N =  100;

    double  ymin = 0.;
    double  ymax = 0.2;

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

    double egam = 10.;
    double W = W_cm(egam);
    double s = W*W;
    double  xmin = - kBox.t_man(s, 0.);
    double  xmax = - kBox.t_man(s, PI);

    auto F = [&](double mt)
    {
        return box3.differential_xsection(s, -mt);
    };  

    plotter->AddEntry(N, F, {xmin,xmax}, "#it{E}_{#gamma} = 10 GeV", PRINT);      
    plotter->SetXaxis("#minus#it{t}  [GeV^{2}]", xmin, xmax);

    egam = 9.5; W = W_cm(egam); s = W*W;  
    xmin = - kBox.t_man(s, 0.); xmax = - kBox.t_man(s, PI);
    plotter->AddEntry(N, F, {xmin,xmax}, "#it{E}_{#gamma} = 9.5 GeV", PRINT);      

    egam = 9.0; W = W_cm(egam); s = W*W;  
    xmin = - kBox.t_man(s, 0.); xmax = - kBox.t_man(s, PI);
    plotter->AddEntry(N, F, {xmin,xmax}, "#it{E}_{#gamma} = 9.0 GeV", PRINT);   

    plotter->SetYaxis(ylabel, ymin, ymax);
    
    plotter->SetLegend(0.25, 0.77);
    plotter->SetLegendOffset(0.5, 0.1);

    // Output to file
    plotter->Plot(filename);

    delete plotter;

};