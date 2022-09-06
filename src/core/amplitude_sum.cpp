// Class to sum up any number of generic amplitudes and build observables.
// Amplitudes are loaded up in a vector and summed incoherently
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitude_sum.hpp"

// Make sure all amplitudes being summed are compatible kineamtically
bool jpacPhoto::amplitude_sum::check_compatibility(amplitude * amp)
{
    bool COMPATIBLE = true;
    std::string reason;

    // Check that all amplitudes in the sum share the same kinematics pointer
    if (amp->_kinematics != _kinematics)
    {
        COMPATIBLE = false;
        reason = "Kinematics pointers do not match up.";
    };

    // Also check that all amplitudes have helicities defined with respect to the same
    // center-of-mass frame
    // TODO: add wigner rotations such that this is not a requirement
    if (_amps.size() > 0)
    {
        COMPATIBLE = (amp->helicity_CM_frame() == _amps[0]->helicity_CM_frame());
        reason = "Helicity frames do not match up.";
    };

    if (!COMPATIBLE)
    {
        std::cout << "Warning! Amplitude (" << amp->get_id() 
            << ") not compatible in sum (" << get_id() << ")." 
            << std::endl; 
        std::cout << reason +  "Skipping..." << std::endl;
    };
            
    return COMPATIBLE;
};

// Evaluate the sum for given set of helicites, and mandelstam invariant s and t
std::complex<double> jpacPhoto::amplitude_sum::helicity_amplitude(std::array<int, 4> helicities, double s, double t)
{
    // Get the index in the cache corresponding to these helicities
    int index = find_helicity(helicities, _kinematics->get_meson_JP()[0],  _kinematics->get_baryon_JP()[0], _kinematics->is_photon());
    
    // Sum over the caches of each constituate ampltiude
    std::complex<double> result = 0.;
    for (int i = 0; i < _amps.size(); i++)
    {
        if (_debug == 1)
        {
            result += _amps[i]->helicity_amplitude(helicities, s, t);
        }
        else
        {
            _amps[i]->update_cache(s, t);
            result += _amps[i]->get_cached_helicity_amplitude(index);
        }
    }

    return result;
};
