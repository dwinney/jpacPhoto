// ---------------------------------------------------------------------------
// Predictions for the Jpsi SDME's near-threhsold 
// 
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// REFERENCES:
// [1] arXiv:2008.01001 [hep-ph]
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "pomeron_exchange.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

void jpsi_SDME()
{

    // ---------------------------------------------------------------------------
    // Preliminaries
    // ---------------------------------------------------------------------------

    // J/Psi kinematics, same for both models
    reaction_kinematics kJpsi (M_JPSI);
    kJpsi.set_meson_JP(1, -1);

    // Low - energy trajectory and couplings from [1]
    linear_trajectory alpha_LE (1, 0.94, 0.36, "LE");
    double A_LE = 0.38;
    double b_LE = 0.12;

    // High - energy trajectory and couplings from [2]
    linear_trajectory alpha_HE (1, 1.15, 0.11, "HE");
    double b_HE = 1.01;
    double A_HE = 0.16;

    // ---------------------------------------------------------------------------
    // Low-Energy Amplitudes
    // ---------------------------------------------------------------------------

    // The false means not helicity conserving version 
    pomeron_exchange jpsi_LE (&kJpsi, &alpha_LE, false, "#it{J}/#psi");
    jpsi_LE.set_params({A_LE, b_LE});

    // flip to true if you want the conserving one
    pomeron_exchange jpsi_HE (&kJpsi, &alpha_HE, true, "#it{J}/#psi");
    jpsi_HE.set_params({A_HE, b_HE});

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------
    
    // Beam energy
    double ebeam = 10.;
    double w = W_cm(ebeam);
    double s = w*w + EPS;

    int N = 50; // Number of points to plot
    double xmin = -kJpsi.t_man(s, 0.) + EPS;               // tmin
    double xmax = std::min(1., -kJpsi.t_man(s, M_PI));     // min(tmax, 1 GeV2)

    double  ymin = -0.6;
    double  ymax =  0.7;

    std::string filename = "jpsi_SDME.pdf";
    std::string ylabel  = "";
    std::string xlabel  = "#minus#it{t}  [GeV^{2}]";
    bool PRINT = true;

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    bool real;
    int alpha, lam, lamp;

    // Solid line entries call jpsi_LE
    auto F = [&](double mt)
    {
        std::complex<double> sdme = jpsi_LE.SDME_GJ(alpha, lam, lamp, s, -mt);
        return real * std::real(sdme) + !real * std::imag(sdme);
    };

    // Dashed line entries calls jpsi_HC
    auto G = [&](double mt)
    {
        std::complex<double> sdme = jpsi_HE.SDME_GJ(alpha, lam, lamp, s, -mt);
        return real * std::real(sdme) + !real * std::imag(sdme);
    };

    // ---------------------------------------------------------------------
    // Add desired sdmes

    real = true;
    alpha = 0; lam = 0; lamp = 0;
    plotter->AddEntry(N, F, {xmin,xmax}, "#rho^{0}_{00}", PRINT);
    plotter->AddDashedEntry(N, G, {xmin,xmax}, PRINT);

    // And others
    real = true;
    alpha = 1; lam = 1; lamp = -1;
    plotter->AddEntry(N, F, {xmin,xmax}, "#rho^{1}_{1-1}", PRINT);
    plotter->AddDashedEntry(N, G, {xmin,xmax}, PRINT);
    
    real = false;
    alpha = 2; lam = 1; lamp = -1;
    plotter->AddEntry(N, F, {xmin,xmax}, "Im #rho^{2}_{1-1}", PRINT);
    plotter->AddDashedEntry(N, G, {xmin,xmax}, PRINT);

    //--------------------------------------------------------------------
    // Make the plot pretty

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);

    // Add a header to legend to specify the fixed energy
    std::ostringstream streamObj;
    streamObj << std::setprecision(4) << "#it{E}_{#gamma} = " << ebeam << " GeV";
    std::string header = streamObj.str();
    plotter->SetLegend(0.22, 0.6, header); 
    plotter->SetLegendOffset(0.3, 0.18);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
};