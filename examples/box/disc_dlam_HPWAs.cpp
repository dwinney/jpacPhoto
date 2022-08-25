// ---------------------------------------------------------------------------
// Read in 2D interpolations of the PWA for b^J and c^J and assemble the 
// imaginary of the partial-wave of the full box amplitude
// 
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "reaction_kinematics.hpp"
#include "interpolation_2D.hpp"

#include "jpacGraph1D.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

void disc_dlam_HPWAs()
{

    double smin = 4.;
    double smax = 5.;

    double ymin = -0.05;
    double ymax =  0.08;

    double eta =  1.;
    bool verbose = true;
    bool PRINT = true;

    int h;
    int J = 1;
    int N = 200;

    std::string filename = "hpwa.pdf";

    // ---------------------------------------------------------------------------
    // Preliminaries 

    // Need the kinematics of the intermediate reactions to get phases and helicity combinations
    reaction_kinematics kBox (M_JPSI, M_PROTON);
    kBox.set_meson_JP(1, -1);
    int nBox = kBox.num_amps();

    // Set up a vector to hold the helicity PWA interpoaltions
    std::vector<interpolation_2D*> ampA;
    for (int i = 0; i < nBox; i++) ampA.push_back( new interpolation_2D(verbose) );

    std::string path = "./grid_data/";
    
    // ---------------------------------------------------------------------------
    // Assemble the box PWAS

    // File name prefix
    std::string prefix = path + "boxD_J_" + std::to_string(J) + "_H_";

    // Grab the grids (remember we only saved half the amplitudes!)
    // first the ImA amps 
    for (int i=0; i < nBox/2; i++)
    {
        std::string filename = prefix + std::to_string(i) + ".dat";
        ampA[i]->import_grid(filename);
    };

    // Plotter object
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    auto F = [&](double w)
    {
        return ampA[h]->eval(w*w, eta);
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