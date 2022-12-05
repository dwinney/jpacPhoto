// Class for the polarization vector of vector particles
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _POLVEC_
#define _POLVEC_

#include <iostream>
#include "constants.hpp"
#include "two_body_state.hpp"

// ---------------------------------------------------------------------------
// Polarization vectors for vector particles
// in the s-channel center of mass frame
// Vector particles are always particle 1
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
    class polarization_vector
    {
        public:

        // Constructor
        polarization_vector(two_body_state * xstate)
            : _state(xstate)
        {};

        // Destructor
        ~polarization_vector(){};

        // Components
        std::complex<double> component(int i, int lambda, double s, double theta);
        std::complex<double> conjugate_component(int i, int lambda, double s, double theta);

        inline std::complex<double> field_tensor(int i, int j, int lambda, double s, double theta)
        {
            std::complex<double> result;
            result  = _state->q(i, s, theta) * component(j, lambda, s, theta);
            result -= _state->q(j, s, theta) * component(i, lambda, s, theta);

            return result;
        };

        private:
        
        two_body_state * _state;
    };
};

#endif
