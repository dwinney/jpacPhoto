// ---------------------------------------------------------------------------
// Script to generate 2D grids of helicity partial waves of gamma p -> D Lam_c
//
// This generates 4 files, one for each independent helicity combination, per J.
// By default this considers up to 2J = 5, so 12 .dat files total
// 
// the grid is 50 points in s in range [4^2, 6^2]
// and 20 points in eta in range [0, 1.5]
// 
// filenames follow convention: gamD_J_%1_H_%2.dat with %1 and %2 the values 2J and helicity index.
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        daniel.winney@gmail.edu
// ---------------------------------------------------------------------------
// References:
// [1] arXiv:2009.08345v1
// ---------------------------------------------------------------------------

#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "pseudoscalar_exchange.hpp"
#include "vector_exchange.hpp"
#include "dirac_exchange.hpp"
#include "amplitude_sum.hpp"

#include "helicity_PWA.hpp"
#include "interpolation_2D.hpp"
#include "projected_amplitude.hpp"

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

void gamp_dslam_grid()
{
    // Form factor parameter
    double lambdaQCD = 0.25;

    double e     = sqrt(4.* PI * ALPHA);
    double gDsDs = 0.641;
    double gDDs  = 0.134;
    double gNDsL = -4.3;
    double gNDL  = -13.2;

    // ---------------------------------------------------------------------------
    // D phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for Dbar LambdaC in final state
    reaction_kinematics kDs (M_DSTAR, M_LAMBDAC);
    kDs.set_meson_JP(1, -1);

    dirac_exchange ds_lamcEx (&kDs, M_LAMBDAC, "#Lambda_{c} exchange");
    ds_lamcEx.set_params({e, gNDsL});
    ds_lamcEx.force_covariant(true);

    vector_exchange ds_dstarEx (&kDs, M_DSTAR, "D* exchange");
    ds_dstarEx.set_params({gDsDs, gNDsL, 0.});
    ds_dstarEx.force_covariant(true);

    pseudoscalar_exchange ds_dEx (&kDs, M_D, "D exchange");
    ds_dEx.set_params({gDDs, gNDL});
    ds_dEx.force_covariant(true);
    
    amplitude_sum ds_sum (&kDs,  {&ds_dEx, &ds_dstarEx, &ds_lamcEx}, "Sum");
    // ---------------------------------------------------------------------------
    // PWA projection 
    // ---------------------------------------------------------------------------

    // We want partial waves up to 5/2
    int Jmax = 1;
    
    // Pass the full amplitude to the helicity pwa amplitude
    helicity_PWA hpwa(&ds_sum);

    // ---------------------------------------------------------------------------
    // Inteprolation
    // ---------------------------------------------------------------------------

    // Prefix identifier for the name of files
    // These are the pw of the amplitude "B"
    std::string prefix = "./grid_data/gamDs";

    //Grid size parameters
    double Wmin = kDs.Wth() + 0.001, Wmax = 6.;
    double etamin = 0., etamax = 1.5;
    int nS = 200, nEta = 20;

    // Interpolation object that actually generates the grid
    interpolation_2D interpolator;
    interpolator.set_verbose(1);
    interpolator.set_limits({Wmin*Wmin, Wmax*Wmax}, {etamin, etamax});
    interpolator.set_grid_size(nS, nEta); // 50 points in s and 20 in eta

    // Loop over only the first half of helicity combinations
    int nAmps = kDs.num_amps();
    for (int j = 0; (2*j+1) <= Jmax; j++)
    {
        int J = 2*j+1;
        hpwa.set_J(J);

        for (int i = 0; i < nAmps/2; i++)
        {
            // Create a helicity partial wave amplitude per helicity combination
            std::array<int,4> ith_helicities = kDs.helicities(i);

            // Pass projection parameters
            hpwa.set_helicities(ith_helicities);

            // Function which will be binned
            auto f = [&](double s, double eta)
            {
                // Pass eta to the form-factors of the constituent amplitudes
                ds_dEx.set_formfactor(    2, M_D       + eta * lambdaQCD);
                ds_dstarEx.set_formfactor(2, M_DSTAR   + eta * lambdaQCD);
                ds_lamcEx.set_formfactor( 2, M_LAMBDAC + eta * lambdaQCD);
                
                return hpwa.real_part(s);
            };
            
            // Finish the filename with J and H index values
            std::string filename = prefix + "_J_" + std::to_string(J) + "_H_" + std::to_string(i) + ".dat";
            
            std::string hset = "{" + std::to_string(ith_helicities[0]) + ", " + std::to_string(ith_helicities[1]) + ", " + std::to_string(ith_helicities[2]) + ", " + std::to_string(ith_helicities[3]) + "} ";
            std::cout << "Generating grid for J = " << std::to_string(J) <<" and H = " + hset << std::endl;

            // Make grid 
            // (int 1 is to skip interpolation step since we only want to print to file not eval)
            interpolator.generate_grid(f, 1);

            // Save it to file
            interpolator.export_grid(filename);
            std::cout << std::endl;
        };
    }

};