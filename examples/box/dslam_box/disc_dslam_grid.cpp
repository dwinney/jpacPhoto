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

#include <cstring>
#include <iostream>
#include <iomanip>

using namespace jpacPhoto;

void disc_dslam_grid()
{
    // ---------------------------------------------------------------------------
    // Preliminaries 

    // How many partial waves to process
    int  Jmax = 1;
    bool verbose = true;

    // Need the kinematics of the intermediate reactions to get phases and helicity combinations
    reaction_kinematics kgamDs (M_DSTAR, M_LAMBDAC);
    kgamDs.set_meson_JP(1, -1);
    int nGamDs = kgamDs.num_amps();

    // Set up a vector to hold the nGamD helicity PWA interpoaltions
    std::vector<interpolation_2D*> ampB;
    for (int i = 0; i < nGamDs; i++) ampB.push_back( new interpolation_2D(verbose) );

    reaction_kinematics kpsiDs (M_JPSI, M_PROTON, M_DSTAR, M_LAMBDAC);
    kpsiDs.set_meson_JP(1, -1);
    int nPsiDs = kpsiDs.num_amps();

    // Also set up containers for the C PWAs
    std::vector<interpolation_2D*> ampC;
    for (int i = 0; i < nPsiDs; i++) ampC.push_back( new interpolation_2D(verbose) );
    
    // path to where the grid files are
    std::string path = "./grid_data/";
    
    // ---------------------------------------------------------------------------
    // Set up kinematics for the overall process
    reaction_kinematics kbox (M_JPSI, M_PROTON);
    kbox.set_meson_JP(1, -1);
    int nBoxAmps = kbox.num_amps();

    // Need two-body phase-space function
    auto rho = [&] (double s)
    {
        if ( s <= kgamDs.sth() ) return 0.;
        double k = kgamDs.final_momentum(s);
        return k / (8. * PI * sqrt(s));
    };

    //Grid size parameters (same as those used in the consitutent grids but this is not necessary)
    double Wmin = sqrt(18.4704) + EPS, Wmax = 6.;
    double etamin = 0., etamax = 1.5;
    int nS = 200, nEta = 20;

    // Interpolation object that actually generates the grid
    interpolation_2D output;
    output.set_verbose(verbose);
    output.set_limits({Wmin*Wmin, Wmax*Wmax}, {etamin, etamax});
    output.set_grid_size(nS, nEta); // 50 points in s and 20 in eta

    // ---------------------------------------------------------------------------
    // Assemble the box PWAS

    // Loop over all desired J values
    for (int j=0; (2*j+1) <= Jmax; j++)
    {   
        int J = 2*j+1;

        // File name prefix
        std::string gamPrefix = path + "gamDs_J_" + std::to_string(J) + "_H_";
        std::string psiPrefix = path + "psiDs_J_" + std::to_string(J) + "_H_";

        // Grab the grids (remember we only saved half the amplitudes!)
        // first the B amps 
        for (int i=0; i < nGamDs/2; i++)
        {
            std::string filename = gamPrefix + std::to_string(i) + ".dat";
            ampB[i]->import_grid(filename);
        };

        // then the C amps
        for (int i=0; i < nPsiDs/2; i++)
        {
            std::string filename = psiPrefix + std::to_string(i) + ".dat";
            ampC[i]->import_grid(filename);
        };

        for (int i = 0; i < nBoxAmps/2; i++)
        {
            std::array<int,4> ith_helicities = kbox.helicities(i);

            // Tricky bit is making sure the correct product of b and c is done 
            // while summing over intermediate helicities
            auto f = [&] (double s, double eta)
            {
                double a = 0.;
                double b = 0., c = 0.;

                for (int n = 0; n < 6; n++)
                {   

                    int b_ind = (i>=6) * 6;
                    b = ampB[b_ind + n]->eval(s, eta);

                    int c_ind = (i - (i>=6)*6) * 6;
                    if (c_ind <= 12)
                    {
                        c = ampC[c_ind + n]->eval(s, eta);
                    } 
                    else
                    {
                        // with phase the kth amplitude gets replaced by N-1-kth amplitude where N is TOTAL number of amplitudes
                        c = kpsiDs.parity_phase(c_ind + n, S)*ampC[nPsiDs -1 - c_ind - n]->eval(s,eta);
                    };

                    a += rho(s) * b * c;
                }

                return a;
            };

            // Finish the filename with J and H index values
            std::string filename = path + "boxDs_J_" + std::to_string(J) + "_H_" + std::to_string(i) + ".dat";
            
            if (verbose)
            {
                std::string hset = "{" + std::to_string(ith_helicities[0]) + ", " + std::to_string(ith_helicities[1]) + ", " + std::to_string(ith_helicities[2]) + ", " + std::to_string(ith_helicities[3]) + "} ";
                std::cout << "Generating grid for J = " << std::to_string(J) <<" and H = " + hset << std::endl;
            }

            // Make grid 
            // (int 1 is to skip interpolation step since we only want to print to file not eval)
            output.generate_grid(f, 1);

            // Save it to file
            output.export_grid(filename);
            if (verbose) std::cout << std::endl;
        };
    };

    // Clean up out pointers
    for (int i = 0; i < ampB.size(); i++) delete ampB[i];
    for (int i = 0; i < ampC.size(); i++) delete ampC[i];
};