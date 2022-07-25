#include "inclusive/triple_regge.hpp"
#include "inclusive/inclusive_kinematics.hpp"
#include "regge_trajectory.hpp"

#include "jpacUtils.hpp"
#include "jpacGraph1D.hpp"

#include <iostream>
#include <cstring>

using namespace jpacPhoto;

int main( int argc, char** argv )
{

    // ---------------------------------------------------------------------------
    // Amplitudes
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // Global quantities

    // Couplings
    double g_b1 = 0.24;
    double g_NN = sqrt(2.) * sqrt(4. * PI * 13.81); // Nucleon coupling same for all
    double LamPi = .9;  // 900 MeV cutoff for formfactor
    
    // Kinematics
    reaction_kinematics * kb1 = new reaction_kinematics(M_B1);
    kb1->set_meson_JP(1, +1);

    // ---------------------------------------------------------------------------
    // Fixed-spin pion amplitude 

    // Exclusive amplitude
    pseudoscalar_exchange * excB1f = new pseudoscalar_exchange(kb1, M_PION, "b1 production");
    excB1f->set_params({g_b1, g_NN});
    excB1f->set_formfactor(true, LamPi);

    // We now can pass this to an inclusive amplitude
    triple_regge * incB1f = new triple_regge(excB1f);
    incB1f->set_high_energy_approximation(false);
  
    // // ---------------------------------------------------------------------------
    // // Plotting options
    // // ---------------------------------------------------------------------------

    int N = 50;

    double xmin = 2.;   
    double xmax = 5. ;

    double ymin = 0.;
    double ymax = 8.;

    std::string filename = "integrated.pdf";
    std::string ylabel   = "#sigma  [#mub]";
    std::string xlabel   = "#it{W}_{#gamma#it{p}}  [GeV]";

    // ---------------------------------------------------------------------------  
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();

    bool addExc;
    auto F = [&](double w)
    {
        return incB1f->integrated_xsection(w*w) + excB1f->integrated_xsection(w*w) * 1.E-3; // in mub!
    };
    auto G = [&](double w)
    {
        return excB1f->integrated_xsection(w*w) * 1.E-3; // in mub!
    };
    auto H = [&](double w)
    {
        return incB1f->integrated_xsection(w*w); // in mub!
    };

    incB1f->set_sigma_total(JPAC_pimp_withResonances);
    std::array<std::vector<double>, 2> x_fx, x_fx2;
    x_fx  = vec_fill(2*N,   F, xmin, 3, 1);
    x_fx2 = vec_fill(N, F, 3.1, xmax, 1);
    for (int i = 0; i < x_fx2[0].size(); i++)
    {
        x_fx[0].push_back(x_fx2[0][i]);
        x_fx[1].push_back(x_fx2[1][i]);
    }
    plotter->AddEntry(x_fx[0], x_fx[1], "total #it{b}_{1}^{#plus} production");
    
    // Add the exclusive curve
    x_fx[0].clear(); x_fx[1].clear();
    x_fx2[0].clear(); x_fx2[1].clear();
    x_fx  = vec_fill(2*N,   G, xmin, 3, 1);
    x_fx2 = vec_fill(N, G, 3.1, xmax, 1);
    for (int i = 0; i < x_fx2[0].size(); i++)
    {
        x_fx[0].push_back(x_fx2[0][i]);
        x_fx[1].push_back(x_fx2[1][i]);
    }
    plotter->AddDashedEntry(x_fx[0], x_fx[1]);

    incB1f->set_sigma_total(JPAC_pipp_withResonances);
    plotter->AddEntry(N, H, {xmin, xmax}, "total #it{b}_{1}^{#minus} production", 1);

    // ---------------------------------------------------------------------------
    // Finally make the plot pretty

    // Axes options
    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);

    // LEgend options
    plotter->SetLegend(0.6, 0.65);
    plotter->SetLegendOffset(0.3, 0.1);

    // Output to file
    plotter->Plot(filename);

    return 0;
};