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

void disc_dlam_grid()
{
    // ---------------------------------------------------------------------------
    // Preliminaries 

    // How many partial waves to process
    int  Jmax = 3;
    bool verbose = true;

    // Need the kinematics of the intermediate reactions to get phases and helicity combinations
    reaction_kinematics kgamD (M_D, M_LAMBDAC);
    kgamD.set_meson_JP(0, -1);
    int nGamD = kgamD.num_amps();

    // Set up a vector to hold the nGamD helicity PWA interpoaltions
    std::vector<interpolation_2D*> ampB;
    for (int i = 0; i < nGamD; i++) ampB.push_back( new interpolation_2D(verbose) );

    reaction_kinematics kpsiD (M_JPSI, M_PROTON, M_D, M_LAMBDAC);
    kpsiD.set_meson_JP(0, -1);
    int nPsiD = kpsiD.num_amps();

    // Also set up containers for the C PWAs
    std::vector<interpolation_2D*> ampC;
    for (int i = 0; i < nPsiD; i++) ampC.push_back( new interpolation_2D(verbose) );
    
    // path to where the grid files are
    std::string inpath = "./grid_data/";
    // std::string outpath = "../../combined/grid_data/";
    std::string outpath = inpath;
    
    
    // ---------------------------------------------------------------------------
    // Set up kinematics for the overall process
    reaction_kinematics kbox (M_JPSI, M_PROTON);
    kbox.set_meson_JP(1, -1);
    int nBoxAmps = kbox.num_amps();

    // Need two-body phase-space function
    auto rho = [&] (double s)
    {
        if ( s <= kgamD.sth() ) return 0.;
        double k = kgamD.final_momentum(s);
        return k / (8. * PI * sqrt(s));
    };

    //Grid size parameters (same as those used in the consitutent grids but this is not necessary)
    double Wmin = kgamD.Wth() + 1.1E-4, Wmax = 6.;
    double etamin = 0.95, etamax = 1.05;
    int nS = 100, nEta = 3;

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
        std::string gamPrefix = inpath + "gamD_J_" + std::to_string(J) + "_H_";
        std::string psiPrefix = inpath + "psiD_J_" + std::to_string(J) + "_H_";

        // Grab the grids (remember we only saved half the amplitudes!)
        // first the B amps 
        for (int i=0; i < nGamD/2; i++)
        {
            std::string filename = gamPrefix + std::to_string(i) + ".dat";
            ampB[i]->import_grid(filename);
        };

        // then the C amps
        for (int i=0; i < nPsiD/2; i++)
        {
            std::string filename = psiPrefix + std::to_string(i) + ".dat";
            ampC[i]->import_grid(filename);
        };

        // Now for each helicity combination in kbox 
        // we grab the grids of the 4 & 6 combinations of kgamD and kpsiD
        for (int i = 0; i < nBoxAmps/2; i++)
        {
            std::array<int,4> ith_helicities = kbox.helicities(i);

            // Tricky bit is making sure the correct product of b and c is done 
            // while summing over intermediate helicities
            auto f = [&] (double s, double eta)
            {
                double a = 0.;

                for (int n = 0; n < 2; n++)
                {
                    double b, c;

                    // if i >= 6 hel[1] flips sign and we add 2 
                    int b_ind = (i>=6) * 2;
                    b  = ampB[b_ind+n]->eval(s, eta);

                    // the psi case is symmetric from 0-5 and 6-11
                    int c_ind = ( (i>=6)*(i-6) + !(i>=6)*i ) * 2;

                    // However we need to incoporate parity phases since we dont have all 12 amplitudes stored
                    if (c_ind <= 4)
                    {
                        c = ampC[c_ind+n]->eval(s, eta);
                    } 
                    else
                    {
                        // with phase the kth amplitude gets replaced by N-1-kth amplitude where N is TOTAL number of amplitudes
                        c = kpsiD.intrinsic_parity(S)*ampC[nPsiD-1-(c_ind+n)]->eval(s,eta);
                    };                    

                    a += rho(s) * b * c;
                };

                return a;
            };

            // Finish the filename with J and H index values
            std::string filename = outpath + "boxD_J_" + std::to_string(J) + "_H_" + std::to_string(i) + ".dat";
            
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