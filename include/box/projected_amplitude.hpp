// Class that calculates the projection of a given exchange amplitude
// onto s-channel total spin-J
//
// This means calculating and storing all helicity combinations for a given J
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef PARTIAL_WAVE
#define PARTIAL_WAVE

#include "helicity_PWA.hpp"

namespace jpacPhoto
{
    class projected_amplitude : public amplitude
    {
        // ---------------------------------------------------------------------------
        public:

        // Constructor 
        // just needs an exchange amplitude and total spin J
        projected_amplitude(amplitude * amp, int J, std::string id = "")
        : amplitude(amp->_kinematics, "projected_amplitude", id), 
         _amplitude(amp), _J(J)
        {
            check_amplitude(amp);
        };

        ~projected_amplitude()
        {
            for (int i = 0; i < _projections.size(); i++) 
            {
                delete _projections[i];
            };
        };

        // Evaluate the helicity amplitude 
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t);

        // Allowed quantum numbers are the same as the amplitude being projected
        inline std::vector<std::array<int,2>> allowed_meson_JP() { return _amplitude->allowed_meson_JP();  };
        inline std::vector<std::array<int,2>> allowed_baryon_JP(){ return _amplitude->allowed_baryon_JP(); };

        // We assume all helicities are in s channel
        inline helicity_channel helicity_CM_frame(){ return S; };

        // ---------------------------------------------------------------------------
        private:

        // Spin we are projecting onto
        int _J; 

        // The total amplitude being projected
        amplitude * _amplitude;
        inline void check_amplitude(amplitude * amp)
        {
            if (amp->helicity_CM_frame() != S)
            {
                std::cout << "Warning! projected_amplitude initialized with an input amplitude whose helicities are not defined in the s-channel. Results may vary... \n";
                std::cout << "Try using force_covariant()." << std::endl;
            }
        };

        // Update the saved interpolations
        void initialize_amplitudes();

        // Vector containing all the independent helicity projection combinations
        std::vector<helicity_PWA*> _projections;
    };
};

#endif