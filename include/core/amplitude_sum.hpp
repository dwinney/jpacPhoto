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
        // -------------------------------------------------------------------
        public:

        // Empty constructor
        amplitude_sum(reaction_kinematics * xkinem, std::string identifer = "amplitude_sum")
        : amplitude(xkinem, "amplitude_sum", identifer)
        {
            set_nParams(0);
        };

        // Constructor with a vector already set up
        amplitude_sum(reaction_kinematics * xkinem, std::vector<amplitude*> vec, std::string identifer = "amplitude_sum")
        : amplitude(xkinem, "amplitude_sum", identifer)
        {
            for (int i = 0; i < vec.size(); i++)
            {
                if (check_compatibility(vec[i]))
                {
                    // Save the amplitude 
                    _amps.push_back(vec[i]);

                    // But also add its parameter details to running total
                    set_nParams(get_nParams() + vec[i]->get_nParams());
                };
            }
        };

        // Add a new amplitude to the vector
        void add_amplitude(amplitude * new_amp)
        {
            if (check_compatibility(new_amp))
            {
                _amps.push_back(new_amp);
                set_nParams(get_nParams() + new_amp->get_nParams());
            };
        };

        // Add all the members of an existing sum to a new sum
        void add_amplitude(amplitude_sum * new_sum)
        {
            for (int i = 0; i < new_sum->_amps.size(); i++)
            {
                // Get each subamplitude in the sum
                amplitude * new_amp = new_sum->_amps[i];

                // Add to existing set
                add_amplitude(new_amp);
            }
        };

        // empty allowedJP, leave the checks to the individual amps instead
        inline std::vector<std::array<int,2>> allowed_meson_JP() { return {}; };
        inline std::vector<std::array<int,2>> allowed_baryon_JP(){ return {}; };

        // Allocate an aggragated vector of parameters to individual amplitudes
        inline void set_params(std::vector<double> x)
        {
            check_nParams(x);
        
            int N = 0;
            // Allocate parameters to the individual amplitudes
            for (int i = 0; i < _amps.size(); i++)
            {
                // Extract the subvector corresponding to the i-th amplitude
                auto start = x.begin() + N;
                auto end   = x.begin() + N + _amps[i]->get_nParams();
                std::vector<double> pars(start, end);

                _amps[i]->set_params(pars);
                N += _amps[i]->get_nParams();
            };

            // At the end check that the number of params allocated is the same as expected
            // However is this happpens the damage is done so we just send out a warning...
            if (N != get_nParams())
           {
                std::cout << "Warning! amplitude_sum::set_params() : Number of parameters allocated doesnt match those expected..." << std::endl;
           };
        };

        // Evaluate the sum for given set of helicites, energy, and cos
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t);

        // Covariant quantities define lambda in S-channel
        // Analytic ones in the T-channel
        helicity_channel helicity_CM_frame()
        {
            // The compatibility check assumes the helicity are in the same frame 
            // So we return the first amplitude

            // CAUTION though, methods like amplitude::force_covariant() can change this after initialization
            return _amps[0]->helicity_CM_frame();
        };

        // -------------------------------------------------------------------
        private:

        // Store a vector of all the amplitudes you want to sum incoherently
        std::vector<amplitude*> _amps;

        // Make sure all amplitudes being summed are compatible kineamtically
        bool check_compatibility(amplitude* amp);
    };
};

#endif
