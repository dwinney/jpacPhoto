// Class that calculates the projection of a given exchange amplitude
// onto s-channel total spin-J
//
// This means calculating and storing all helicity combinations for a given J
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "projected_amplitude.hpp"

// ---------------------------------------------------------------------------
// These functions are for if you want to treat the projection as a full amplitude

// Evaluate the helicity amplitude 
std::complex<double> jpacPhoto::projected_amplitude::helicity_amplitude(std::array<int, 4> helicities, double s, double t)
{
    int n_amps = _kinematics->num_amps();

    if ( _projections.size() != n_amps ) initialize_amplitudes();

    // Find where the desired helicities is in the vectors
    double mJ     = _kinematics->get_meson_JP()[0];
    double bJ     = _kinematics->get_baryon_JP()[0];
    bool massless = _kinematics->is_photon();

    int hel_id = find_helicity(helicities, mJ, bJ, massless);

    std::complex<double> pwa;
    if ( hel_id < n_amps / 2 )
    {
        pwa = _projections[hel_id]->helicity_partial_wave(s);
    }
    else 
    {
        pwa = parity_phase(helicities) * _projections[hel_id - n_amps/2]->helicity_partial_wave(s);
    }

    // Net helicities
    int lam  = 2 * helicities[0] - helicities[1]; // Photon - Target
    int lamp = 2 * helicities[2] - helicities[3]; // Meson  - Recoil

    double theta = _kinematics->theta_s(s, t);

    return double(_J + 1) * wigner_d_half(_J, lam, lamp, theta) * pwa;
};

// ---------------------------------------------------------------------------
// Delete all the saved projections and re-do them in case parameters have changed
void jpacPhoto::projected_amplitude::initialize_amplitudes()
{
    if (_projections.size() > 0)
    {
        for (int i = 0; i < _projections.size(); i++)
        {
            delete _projections[i];
        }
        _projections.clear();
    }

    for (int i = 0; i < _kinematics->num_amps() / 2; i++)
    {
        std::array<int,4> helicities = _kinematics->helicities(i);
        _projections.push_back( new helicity_pwa(_amplitude, _J, helicities) );
    };
};