// Class to sum up any number of generic amplitudes and build observables.
// Amplitudes are loaded up in a vector and summed incoherently
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _SUM_
#define _SUM_

#include "amplitude.hpp"

// ---------------------------------------------------------------------------
// The amplitude_sum class can take a vector of the above amplitude objects
// and add them together to get observables!
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
    class amplitude_sum : public amplitude
    {
        private:
        // Store a vector of all the amplitudes you want to sum incoherently
        std::vector<amplitude*> _amps;

        // Make sure all amplitudes being summed are compatible kineamtically
        bool check_compatibility(amplitude* amp);

        public:
        // Empty constructor
        amplitude_sum(reaction_kinematics * xkinem, std::string identifer = "amplitude_sum")
        : amplitude(xkinem, "amplitude_sum", identifer)
        {
            _isSum = true;
        };

        // Constructor with a vector already set up
        amplitude_sum(reaction_kinematics * xkinem, std::vector<amplitude*> vec, std::string identifer = "amplitude_sum")
        : amplitude(xkinem, "amplitude_sum", identifer)
        {
            _isSum = true;
            for (int i = 0; i < vec.size(); i++)
            {
                if (check_compatibility(vec[i]))
                {
                    _amps.push_back(vec[i]);
                };
            }
        };

        // Add a new amplitude to the vector
        void add_amplitude(amplitude * new_amp)
        {
            if (check_compatibility(new_amp))
            {
                _amps.push_back(new_amp);
            };
        };

        // Add all the members of an existing sum to a new sum
        void add_amplitude(amplitude_sum * new_sum)
        {
            for (int i = 0; i < new_sum->_amps.size(); i++)
            {
                if (check_compatibility(new_sum->_amps[i]))
                {
                    _amps.push_back(new_sum->_amps[i]);
                };
            }
        };

        // empty allowedJP, leave the checks to the individual amps instead
        inline std::vector<std::array<int,2>> allowed_meson_JP() { return {}; };
        inline std::vector<std::array<int,2>> allowed_baryon_JP(){ return {}; };

        // TODO: Add a set_params which timesi in one vector and allocates approriaten number of
        // params to each sub amplitude

        // Evaluate the sum for given set of helicites, energy, and cos
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t);

        // Covariant quantities define lambda in S-channel
        // Analytic ones in the T-channel
        helicity_channel helicity_CM_frame()
        {
            return _amps[0]->helicity_CM_frame();
        };
    };
};

#endif
