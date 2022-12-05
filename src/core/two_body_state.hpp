// Energy and momentum of a two-particle state in the center of mass scattering frame
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _TWO_BODY_
#define _TWO_BODY_

#include <string>
#include <complex>
#include <iostream>

#include "reaction_kinematics.hpp"

// ---------------------------------------------------------------------------
// The two_body_state is the base object for defining a reaction in the
// s-channel center of mass scatering frame.
//
// Two particles of mass mV and mB are defined with momenta opposite along the
// same axis such that the energy and momenta of both particles is entirely
// determined by the center-of-mass energy, s, and the cosing of the angle
// from the z-axis (define to be at theta = 0).
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
  class two_body_state
  {
        private:

        bool _photon = false;
        double _mM2; // Meson mass
        double _mB2; // Baryon mass (allowed to be float for N* or Δ)

        public:

        // Constructor
        two_body_state(double m, double b)
        : _mM2(m*m), _mB2(b*b)
        {
            if (m < 1.E-3) _photon = true;
        };

        // Whether or not our meson is massless photon
        inline bool if_photon(){ return _photon; };

        // return mass
        inline double get_meson_mass() 
        { 
            if (_mM2 >= 0.) 
            {
                return sqrt(_mM2);
            }
            else
            {
                return sqrt(-_mM2);
            }
        };

        inline double get_baryon_mass() { return sqrt(_mB2); };

        // set masses independently
        inline void set_meson_mass(double m)
        {
            _mM2 = m*m;
        };

        inline void set_baryon_mass(double m)
        {
            _mB2 = m*m;
        };

        // Momenta
        // meson is always particle 1 in + z direction, 
        inline std::complex<double> momentum(double s)
        {
            return sqrt( Kallen(XR * s, XR *_mM2, XR * _mB2)) / (2. * sqrt(XR * s));
        };

        // Energies
        inline std::complex<double> meson_energy(double s)
        {
            return (s + _mM2 - _mB2) / (2. * sqrt(XR * s));
        };

        inline std::complex<double> baryon_energy(double s)
        {
            return (s - _mM2 + _mB2) / (2. * sqrt(XR * s));
        };

        // Full 4-momenta 
        std::complex<double> q(int mu, double s, double theta); // 4vector of vector, particle 1
        std::complex<double> p(int mu, double s, double theta); // 4vector of baryon, particle 2
    };
};

#endif
