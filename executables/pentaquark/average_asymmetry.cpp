#include "constants.hpp"
#include "amplitudes/reaction_kinematics.hpp"
#include "amplitudes/baryon_resonance.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // ---------------------------------------------------------------------------
    // AMPLITUDES
    // ---------------------------------------------------------------------------

    // Set up Kinematics for jpsi in final state
    auto ptr = new reaction_kinematics(M_JPSI);
    ptr->set_JP(1, -1);

    // ---------------------------------------------------------------------------
    // T - CHANNEL // this is the same for all cases

    // Set up pomeron trajectory
    // best fit values from [1]
    auto alpha = new linear_trajectory(+1, 0.941, 0.364);

    // Create amplitude with kinematics and trajectory
    auto background = new pomeron_exchange(ptr, alpha, false, "Background");

    // normalization and t-slope
    // best fit values from [1]
    std::vector<double> back_params = {0.379, 0.12};
    background->set_params(back_params);

    // ---------------------------------------------------------------------------
    // S - CHANNEL 

    // Pc4312
    auto Pc4312_1m = new baryon_resonance(ptr, 1, -1, 4.3119, 9.8E-3, "P_{c}(4312)");
    Pc4312_1m->set_params({0.01, .7071});

    auto Pc4312_3m = new baryon_resonance(ptr, 3, -1, 4.3119, 9.8E-3, "P_{c}(4312)");
    Pc4312_3m->set_params({0.01, .7071});

    // Pc4440
    auto Pc4440_1m = new baryon_resonance(ptr, 1, -1, 4.4403, 20.6E-3, "P_{c}(4440)");
    Pc4440_1m->set_params({0.01, .7071});

    auto Pc4440_3m = new baryon_resonance(ptr, 3, -1, 4.4403, 20.6E-3, "P_{c}(4440)");
    Pc4440_3m->set_params({0.01, .7071});

    auto Pc4440_3p = new baryon_resonance(ptr, 3, +1, 4.4403, 20.6E-3, "P_{c}(4440)");
    Pc4440_3p->set_params({0.01, .7071});

    // Pc4457
    auto Pc4457_1m = new baryon_resonance(ptr, 1, -1, 4.4573, 6.4E-3, "P_{c}(4457)");
    Pc4457_1m->set_params({0.01, .7071});

    auto Pc4457_3m = new baryon_resonance(ptr, 3, -1, 4.4573, 6.4E-3, "P_{c}(4457)");
    Pc4457_3m->set_params({0.01, .7071});

    auto Pc4457_5m = new baryon_resonance(ptr, 5, 1, 4.4573, 6.4E-3, "P_{c}(4457)");
    Pc4457_5m->set_params({0.01, .7071});


    // Add to the sum
    auto sumA = new amplitude_sum (ptr, {background, Pc4312_1m, Pc4440_3m, Pc4457_1m}, "A");
    auto sumB = new amplitude_sum (ptr, {background, Pc4312_3m, Pc4440_1m, Pc4457_3m}, "B");
    auto sumC = new amplitude_sum (ptr, {background, Pc4312_3m, Pc4440_3p, Pc4457_5m}, "C");

    // ---------------------------------------------------------------------------
    // Choose which scenario to plot
    std::vector<amplitude*> amps;
    amps.push_back(background);
    amps.push_back(sumA);
    amps.push_back(sumB);
    amps.push_back(sumC);

    int N = 100;
    bool PRINT_TO_COMMANDLINE = true;

    double ymax =  0.2;
    double ymin = -0.05;

    double Emin = E_beam(ptr->Wth()) + EPS;
    double Emax = 12.;

    std::string filename = "sigma_integrated.pdf";

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter objects
    jpacGraph1D* plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // scan over theta
    for (int n = 0; n < amps.size(); n++)
    {
        std::cout << "\nPrinting amplitude " << amps[n]->_identifier << ".\n";

        auto F = [&](double theta)
        {

            auto f = [&](double E)
            {
                double W = W_cm(E);
                double t = ptr->t_man(W*W, theta * DEG2RAD);
                return amps[n]->beam_asymmetry_4pi(W*W, t);
            };

            ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
            ROOT::Math::Functor1D wF(f);
            ig.SetFunction(wF);
            
            return ig.Integral(Emin, Emax) / (Emax - Emin);
        };

        std::array<std::vector<double>, 2> x_fx = vec_fill(N, F, 0., 90., PRINT_TO_COMMANDLINE);
        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->_identifier);
    }

    // Add a header to legend to specify the fixed energy
    plotter->SetLegend(0.2, 0.7);

    plotter->SetXaxis("#theta", 0., 90.);

    // To change the range of the Y-axis or the position of the Legend change the arguments here
    plotter->SetYaxis("#Sigma averaged over E_{#gamma}", ymin, ymax);

    std::string file1 = "theta_" + filename;
    plotter->Plot(file1.c_str());

    // clear
    plotter->ClearData();
    
    // ---------------------------------------------------------------------------
    // scan over energy
    for (int n = 0; n < amps.size(); n++)
    {
        std::cout << "\nPrinting amplitude " << amps[n]->_identifier << ".\n";

        auto F = [&](double egam)
        {
            double W = W_cm(egam);
            auto f = [&](double t)
            {
                return amps[n]->beam_asymmetry_4pi(W*W, t);
            };

            ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS61);
            ROOT::Math::Functor1D wF(f);
            ig.SetFunction(wF);
            
            double t_min = amps[n]->_kinematics->t_man(W*W, 0.);
            double t_max = amps[n]->_kinematics->t_man(W*W, PI);
            return ig.Integral(t_max, t_min) / (t_min - t_max);
        };

        std::array<std::vector<double>, 2> x_fx = vec_fill(N, F, Emin, Emax, PRINT_TO_COMMANDLINE);
        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->_identifier);
    }

    // Add a header to legend to specify the fixed energy
    plotter->SetLegend(0.2, 0.7);

    plotter->SetXaxis("E_{#gamma}", Emin, Emax);

    // To change the range of the Y-axis or the position of the Legend change the arguments here
    plotter->SetYaxis("#Sigma averaged over t", ymin, ymax);

    std::string file2 = "egam_" + filename;
    plotter->Plot(file2.c_str());

    return 0;
};
