#include "inclusive/triple_regge.hpp"
#include "inclusive/inclusive_kinematics.hpp"
#include "regge_trajectory.hpp"

#include "jpacUtils.hpp"
#include "jpacGraph1D.hpp"

#include <iostream>
#include <cstring>
#include <memory>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // First we set up the exclusive amplitude

    // Couplings
    double g_b1 = 0.24;
    double g_NN = sqrt(2.) * sqrt(4. * PI * 13.81); // Nucleon coupling same for all
    double LamPi = .9;  // 900 MeV cutoff for formfactor
    double bPi = 1. / (LamPi * LamPi);
    
    // Kinematics
    reaction_kinematics * kb1 = new reaction_kinematics(M_B1);
    kb1->set_JP(1, +1);
    
    // Pion trajectory 
    int signature = +1;
    double alpha_prime = 0.7; // GeV^-2
    double alpha_0 =  - alpha_prime * M2_PION;
    std::unique_ptr<linear_trajectory> alpha( new linear_trajectory(signature, alpha_0, alpha_prime));

    // Amplitude object
    std::unique_ptr<pseudoscalar_exchange> excB1( new pseudoscalar_exchange(kb1, alpha.get(), "b1 production") );
    excB1->set_params({g_b1, g_NN});
    excB1->set_formfactor(true, bPi);

    // We now can pass this to an inclusive amplitude
    std::unique_ptr<triple_regge> incB1( new triple_regge(excB1.get()));
    
    // // Set up the kinematics
    // std::unique_ptr<inclusive_kinematics> kinematics( new inclusive_kinematics(M_B1));

    // // Now we set up the triple regge amplitude
    // std::unique_ptr<triple_regge> pi_exchange( new triple_regge(M_B1, alpha.get(), "pion exchange"));
    // // std::unique_ptr<triple_regge> pi_exchange( new triple_regge(M_B1, M_PION, "pion exchange"));

    // // Coupling function
    // double gb1 = 0.24;
    // auto coupling = [&] (double t)
    // {
    //     return (gb1 / M_B1) * (t - M2_B1);
    // };
    // pi_exchange->set_coupling(coupling);
    
    // // Total hadronic cross-ection
    // std::unique_ptr<sigma_tot> sigma_tot_pimp( new sigma_tot_PDG(M_PION, M_PROTON, {-1., 1., 9.56, 1.767, 18.75}, "rpp2020-pimp_total.dat"));
    // pi_exchange->set_sigma_tot(sigma_tot_pimp.get());

    // // Form factor 
    // pi_exchange->set_form_factor( 1./(LamPi*LamPi) );

    // // Specifiy that we want to use the high-energy approximation 
    // pi_exchange->set_high_energy_approximation(true);

    incB1->set_high_energy_approximation(true);
    
    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    int N = 200;
    
    double s = 75.9421;

    double xmin = 0.6;   
    double xmax = 0.99;

    std::string filename = "dsigmadx.pdf";
    std::string ylabel   = "d#sigma / dx  [#mub]";
    std::string xlabel   = "x";

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    std::unique_ptr<jpacGraph1D> plotter( new jpacGraph1D() );

    // ---------------------------------------------------------------------------
    // Print the phase-space for each kinematics

    auto F = [&](double x)
    {
        return incB1->dsigma_dx(s, x) * 1.E3; // in mub!
    };

    std::array<std::vector<double>, 2> x_fx; 
    x_fx = vec_fill(N, F, xmin, xmax, true);

    plotter->AddEntry(x_fx[0], x_fx[1], "");

    plotter->SetXaxis(xlabel, xmin, 1.);
    plotter->SetYaxis(ylabel);

    plotter->SetLegend(false);

    // Output to file
    plotter->Plot(filename);
    return 1.;
};