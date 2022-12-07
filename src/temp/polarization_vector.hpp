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
#include "four_momenta.hpp"

// ---------------------------------------------------------------------------
// Polarization vectors for spin-1 particles
// in the s-channel center of mass frame
// Mesons are always "particle 1"
// ---------------------------------------------------------------------------

namespace jpacPhoto
{

    class polarization_vector
    {
        public:

        // Constructor
        polarization_vector(four_momenta xstate, bool OUTGOING = false)
        : _state(xstate), _star(OUTGOING)
        {};

        // Individual Lorentz components
        complex component(lorentz_index i, int lambda, double s, double theta);

        private:
        
        // This saves the momenta and energy info
        four_momenta _state;

        // Whether this is a eps* (outgoing particle) or simply eps
        bool _star = false;
    };
};

#endif
