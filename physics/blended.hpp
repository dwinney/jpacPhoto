// If a process is described by to different models at low and high energy,
// this function allows your to blend the two together to a simple "all energies"
// interpolation
//
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
//               University of Bonn (HISKP)
// Email:        daniel.winney@gmail.com  
// ---------------------------------------------------------------------------

#ifndef BLENDED_HPP
#define BLENDED_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto
{
    class blended : public raw_amplitude
    {
        public:

        blended(amplitude_key key, kinematics xkinem, std::array<amplitude,2> amps, std::array<double,2> bounds, std::string id = "blended")
        : raw_amplitude(key, xkinem, id), _low(amps[0]), _high(amps[1]), _bounds(bounds)
        {
            initialize(0);
            if (bounds[0] > bounds[1]) warning("jpacPhoto::blended", "Lower limit bigger than upper limit!");
            if ((amps[0]->_kinematics != amps[1]->_kinematics) || (xkinem != amps[0]->_kinematics))
            {
                warning("jpacPhoto::blended", "Kinematics instances not compatible with each other!");
            };
        };

        // Calculate each helicity amplitude 
        inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
        {
            int index = find_helicity(helicities, _kinematics->get_meson_JP()[0], _kinematics->get_baryon_JP()[0], _kinematics->is_photon());
            if (s < _bounds[0]) return _low->helicity_amplitude(helicities, s, t);
            if (s > _bounds[1]) return _high->helicity_amplitude(helicities, s, t);

            double blend_range  = _bounds[1] - _bounds[0];
            double fraction_low = 1 - (s - _bounds[0])/blend_range;

            complex low_contrib  = _low->helicity_amplitude(helicities,  s, t) * fraction_low;
            complex high_contrib = _high->helicity_amplitude(helicities, s, t) * (1 - fraction_low);

            return low_contrib + high_contrib;
        };

        // Get all info from _low since both amplitudes SHOULD be compatible
        inline helicity_frame native_helicity_frame(){  return _low->native_helicity_frame(); };
        inline std::vector<quantum_numbers> allowed_mesons() { return _low->allowed_mesons(); };
        inline std::vector<quantum_numbers> allowed_baryons(){ return _low->allowed_baryons(); };

        private:

        amplitude _low, _high;

        // Bound up to which to use _low and from which to use _high
        // in between we linear interpolate
        std::array<double, 2> _bounds;
    };
};

#endif