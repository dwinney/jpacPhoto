// Class to sum up any number of generic amplitudes and build observables.
// Amplitudes are loaded up in a vector and summed incoherently
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitudes/amplitude_sum.hpp"

// Make sure all amplitudes being summed are compatible kineamtically
bool jpacPhoto::amplitude_sum::check_compatibility(amplitude * amp)
{
    bool COMPATIBLE = true;
    if (amp->_kinematics != _kinematics)
    {
        std::cout << "Warning! Amplitude (" << amp->get_id() 
            << ") not compatible kinematically in sum (" << get_id() << "). Skipping..." 
            << std::endl; 
        COMPATIBLE = false;
    };

    return COMPATIBLE;
};

// Evaluate the sum for given set of helicites, and mandelstam invariant s and t
std::complex<double> jpacPhoto::amplitude_sum::helicity_amplitude(std::array<int, 4> helicities, double s, double t)
{
    // Get the index in the cache corresponding to these helicities
    int index = find_helicity(helicities, _kinematics->_jp[0], _kinematics->get_beam_mass());
    
    // Sum over the caches of each constituate ampltiude
    std::complex<double> result = 0.;
    for (int i = 0; i < _amps.size(); i++)
    {
        _amps[i]->update_cache(s, t);
        result += _amps[i]->get_cached_helicity_amplitude(index);
    }

    return result;
};
