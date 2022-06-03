// Class for the spinor of a massive spin-3/2 particle
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef THREEHALVES
#define THREEHALVES

#include <iostream>
#include "constants.hpp"
#include "two_body_state.hpp"
#include "polarization_vector.hpp"
#include "dirac_spinor.hpp"

// ---------------------------------------------------------------------------
// Rarita spinor combines a spin-1 polarization vector and spin-1/2 dirac spinor
// coupled together with appropriate Clebsh-Gordon coefficients
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
    class rarita_spinor
    {
        public:

        // Constructor
        rarita_spinor(two_body_state * xstate)
        {
            initialize(xstate);
        }; 

        ~rarita_spinor()
        {
            delete _particle_one;
            delete _spin_one;
            delete _spin_half;
        };

        // Access individual components
        std::complex<double> component        (int i, int mu, int lambda, double s, double theta);
        std::complex<double> adjoint_component(int i, int mu, int lambda, double s, double theta);

        private:
        
        // Need two two_body_state objects because we need to be able to 
        // define a polarization_vector with the baryon mass
        two_body_state * _particle_one;
        two_body_state * _particle_two;

        // Make sure we get all the masses allocated correctly.
        void initialize(two_body_state * instate)
        {
            // Assume the inpue is a "normal" kinematics with meson is particle 1 and baryon in particle 2
            _particle_two = instate;

            double meson_mass  = instate->get_meson_mass();
            double baryon_mass = instate->get_baryon_mass();

            // Create a new two_body_state with the masses reversed
            _particle_one = new two_body_state(baryon_mass, meson_mass);

            // The spinor component comes from particle 2
            _spin_half = new dirac_spinor(_particle_two);

            // The "vector" component comes from particle 1
            _spin_one = new polarization_vector(_particle_one);
        };

        // Constituate objects
        polarization_vector * _spin_one;
        dirac_spinor        * _spin_half;
    };
};

#endif
