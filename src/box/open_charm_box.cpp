// Amplitude class specifically for box diagrams through gam p -> DLambda_c -> psi p
// 
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "open_charm_box.hpp"

// ---------------------------------------------------------------------------
// These first methods are of the generic open_charm box

// Generate the 2D interpolation grid for the discontinuity of the box partial waves
// This first generates 2D interpolation grids for the b and c helicity partial waves (gamma/psi p -> open charm).
void jpacPhoto::open_charm_box::generate_grid(double Wmax, std::array<double,2> eta_bounds, std::array<int,2> ns, std::string path)
{
    // Save the grid parameters because these will be used to generate sub-amplitude grids as well
    _wbounds = {_m + M_LAMBDAC + _thOffset, Wmax}; _etabounds = eta_bounds;
    _ns = ns; _path = path;

    // Generate the gamma p and psi p grids
    generate_gamp_grid();
    generate_psip_grid();

    // Set up containers for importing back in the newly generated grids

    // vector full of the B amplitudes (photon)
    std::vector<interpolation_2D*> ampB;
    for (int i = 0; i < _nGam; i++) ampB.push_back( new interpolation_2D(_verbose) );

    // Also set up containers for the C amplitude (jpsi)
    std::vector<interpolation_2D*> ampC;
    for (int i = 0; i < _nPsi; i++) ampC.push_back( new interpolation_2D(_verbose) );
    
    // Interpolation object that actually generates the grid
    interpolation_2D output;
    output.set_verbose(_verbose);
    output.set_limits({_wbounds[0]*_wbounds[0], _wbounds[1]*_wbounds[1]}, {_etabounds[0], _etabounds[1]});
    output.set_grid_size(_ns[0], _ns[1]);

    // Loop over all desired spin's (J <= Jmax)
    for (int j=0 ; (2*j+1) <= _Jmax; j++)
    {   
        // Spin of the PWA
        int J = 2*j+1;

        // File name prefix
        std::string gamPrefix = _path + _file_id + "_gamp_J_" + std::to_string(J) + "_H_";
        std::string psiPrefix = _path + _file_id + "_psip_J_" + std::to_string(J) + "_H_";

        // Import the appropriate amplitudes we made
        // first the B amps 
        for (int i=0; i < _nGam; i++)
        {
            ampB[i]->import_grid(gamPrefix + std::to_string(i) + ".dat");
        };

        // then the C amps
        for (int i=0; i < _nPsi; i++)
        {
            ampC[i]->import_grid(psiPrefix + std::to_string(i) + ".dat");
        };

        // Now we sum over intermediate helicities for each set of the external helicities 
        for (int i = 0; i < _nBox; i++)
        {
            // External (gamma p -> jpsi p) helicities
            std::array<int,4> ith_helicities = _kinematics->helicities(i);

            // Function to compute the discontinutiy at fixed s and eta 
            // This combines the imported grids B and C together
            auto f = [&] (double s, double eta)
            {
                // Output
                double Ima = 0.;

                // sum over all possible intermediate state helicity combos
                int _nIntermediate = (2*_kGam->get_meson_JP()[0]+1)*2; // Charmed meson * Lambda spin combinations
                for (int n = 0; n < _nIntermediate; n++)
                {
                    double b, c;

                    // if i >= 6 hel[1] flips sign and we add 2 
                    int b_ind = (i>=_nBox/2) * _nIntermediate;
                    b  = ampB[b_ind+n]->eval(s, eta);

                    // the psi case is symmetric from 0-5 and 6-11
                    int c_ind = ( (i>=_nBox/2)*(i-_nBox/2) + !(i>=_nBox/2)*i )*_nIntermediate;

                    // However we need to incoporate parity phases since we dont have all 12 amplitudes stored
                    if (c_ind <= 2*_nIntermediate)
                    {
                        c = ampC[c_ind+n]->eval(s, eta);
                    } 
                    else
                    {
                        // with phase the kth amplitude gets replaced by N-1-kth amplitude where N is TOTAL number of amplitudes
                        c = _kPsi->intrinsic_parity(S)*ampC[2*_nPsi-1-(c_ind+n)]->eval(s,eta);
                    };                    

                    Ima += phase_space(s) * b * c;
                };

                return Ima;
            };

            // Finish the filename with J and H index values
            std::string filename = path +  _file_id + "_disc_J_" + std::to_string(J) + "_H_" + std::to_string(i) + ".dat";
            
            if (_verbose)
            {
                std::cout << "Generating grid for J = " << std::to_string(J) <<" and H = " + print_helicities(ith_helicities) << std::endl;
            }

            // Make grid 
            // (int 1 is to skip interpolation step since we only want to print to file not eval)
            output.generate_grid(f, 1);

            // Save it to file
            output.export_grid(filename);
            if (_verbose) std::cout << std::endl;
        };
    };

    // Clean up the pointers we made
    for (int i = 0; i < ampB.size(); i++) delete ampB[i];
    for (int i = 0; i < ampC.size(); i++) delete ampC[i];

    // Finally, import the disc grids we just made
    import_grid(path);
};

// ---------------------------------------------------------------------------
// Methods to generate 2D interpolation grids for the sub-process amplitudes 

void jpacPhoto::open_charm_box::generate_gamp_grid()
{
    // Interpolation object to generate the output files
    interpolation_2D interpolator;
    interpolator.set_verbose(_verbose);
    interpolator.set_limits({_wbounds[0]*_wbounds[0], _wbounds[1]*_wbounds[1]}, {_etabounds[0], _etabounds[1]});
    interpolator.set_grid_size(_ns[0], _ns[1]);

    // We will be taking partial-wave projections so pass the photon amplitude to a projector object
    helicity_PWA hpwa(_gamp_amp);

    // sum over J values untim max is reached
     for (int j = 0; (2*j+1) <= _Jmax; j++)
    {
        int J = 2*j+1; hpwa.set_J(J);

        // Sum over the helicity combinations of gamma p -> D(*) Lambda
        for (int i = 0; i < _nGam; i++)
        {
            // Create a helicity partial wave amplitude per helicity combination
            std::array<int,4> ith_helicities = _kGam->helicities(i);

            // Pass projection parameters
            hpwa.set_helicities(ith_helicities);

            // Function which will be binned
            auto f = [&](double s, double eta)
            {
                update_gamp(eta);
                if (_charmedMeson == D) return hpwa.real_part(s);
                else                    return hpwa.imag_part(s);                
            };
            
            // Finish the filename with J and H index values
            std::string filename = _path + _file_id + "_gamp_J_" + std::to_string(J) + "_H_" + std::to_string(i) + ".dat";
            
            std::cout << "Generating grid for J = " << std::to_string(J) <<" and H = " + print_helicities(ith_helicities) << std::endl;

            // Make grid 
            // (int 1 is to skip interpolation step since we only want to print to file not eval)
            interpolator.generate_grid(f, 1);

            // Save it to file
            interpolator.export_grid(filename);
            std::cout << std::endl;
        };
    };
};

void jpacPhoto::open_charm_box::generate_psip_grid()
{
    // Interpolation object to generate the output files
    interpolation_2D interpolator;
    interpolator.set_verbose(_verbose);
    interpolator.set_limits({_wbounds[0]*_wbounds[0], _wbounds[1]*_wbounds[1]}, {_etabounds[0], _etabounds[1]});
    interpolator.set_grid_size(_ns[0], _ns[1]);

    // We will be taking partial-wave projections so pass the photon amplitude to a projector object
    helicity_PWA hpwa(_psip_amp);

    // sum over J values untim max is reached
     for (int j = 0; (2*j+1) <= _Jmax; j++)
    {
        int J = 2*j+1; hpwa.set_J(J);

        // Sum over the helicity combinations of psi p -> D(*) Lambda
        for (int i = 0; i < _nPsi; i++)
        {
            // Create a helicity partial wave amplitude per helicity combination
            std::array<int,4> ith_helicities = _kPsi->helicities(i);

            // Pass projection parameters
            hpwa.set_helicities(ith_helicities);

            // Function which will be binned
            auto f = [&](double s, double eta)
            {
                update_psip(eta);
                if (_charmedMeson == D) return hpwa.real_part(s);
                else                    return hpwa.imag_part(s);                
            };
            
            // Finish the filename with J and H index values
            std::string filename = _path +  _file_id + "_psip_J_" + std::to_string(J) + "_H_" + std::to_string(i) + ".dat";
            
            std::cout << "Generating grid for J = " << std::to_string(J) <<" and H = " + print_helicities(ith_helicities) << std::endl;

            // Make grid 
            // (int 1 is to skip interpolation step since we only want to print to file not eval)
            interpolator.generate_grid(f, 1);

            // Save it to file
            interpolator.export_grid(filename);
            std::cout << std::endl;
        };
    };
};
