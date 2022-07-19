// ---------------------------------------------------------------------------
// Prediction for Y(4260) and Psi(1S and 2S) based on effective vector
// pomeron exchange at low enegies.
//
// Reproduces left plot in FIG 5 of [1] 
// 
// USAGE:
// make Y_low && ./Y_low
//
// OUTPUT:
// Y_LE.pdf
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// REFERENCES:
// [1] arXiv:2008.01001 [hep-ph]
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "amplitudes/reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/pomeron_exchange.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

    // ---------------------------------------------------------------------------
    // Preliminaries
    // ---------------------------------------------------------------------------

    // Low - energy trajectory and couplings
    linear_trajectory * alpha_LE = new linear_trajectory(1, 0.94, 0.36, "LE");
    double A_LE = 0.38;
    double b_LE = 0.12;

    // J/Psi
    reaction_kinematics * kJpsi = new reaction_kinematics(M_JPSI);
    kJpsi->set_JP(1, -1);

    // ---------------------------------------------------------------------------
    // Low-Energy Amplitudes
    // ---------------------------------------------------------------------------

    // The false means not helicity conserving version 
    // flip to true if you want the conserving one
    pomeron_exchange * jpsi_LE = new pomeron_exchange(kJpsi, alpha_LE, false, "#it{J}/#psi");
    jpsi_LE->set_params({A_LE, b_LE});

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<amplitude*> amps;

    // Options

    int N = 50;
    
    // forward SDME

    double ebeam = 10.;
    double w = W_cm(ebeam);
    double s = w*w + EPS;

    // Plotting ranges
    double xmin = -kJpsi->t_man(s, 0.) + EPS;
    double xmax = std::min(1., -kJpsi->t_man(s, M_PI));

    double  ymin = -0.05;
    double  ymax =  0.1;

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
    auto F = [&](double mt)
    {
        std::complex<double> sdme = jpsi_LE->SDME(alpha, lam, lamp, s, -mt);
        return real * std::real(sdme) + !real * std::imag(sdme);
    };

    // Add upolarized
    real = true;
    alpha = 0; lam = 0; lamp = 0;
    plotter->AddEntry(N, F, {xmin,xmax}, "#rho^{0}_{00}", PRINT);

    // And others
    real = true;
    alpha = 1; lam = 1; lamp = 0;
    plotter->AddEntry(N, F, {xmin,xmax}, "#rho^{1}_{10}", PRINT);
    
    real = false;
    alpha = 1; lam = 1; lamp = -1;
    plotter->AddEntry(N, F, {xmin,xmax}, "#rho^{1}_{1-1}", PRINT);
    
    // Make plot pretty
    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);

    // Add a header to legend to specify the fixed energy
    std::ostringstream streamObj;
    streamObj << std::setprecision(4) << "#it{E}_{#gamma} = " << ebeam << " GeV";
    std::string header = streamObj.str();
    plotter->SetLegend(0.22, 0.65, header); 
    plotter->SetLegendOffset(0.3, 0.18);

    // Output to file
    plotter->Plot(filename);

    return 0;
};