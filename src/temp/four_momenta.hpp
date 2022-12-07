// This provides implementation of the 4-momenta of a two-body system in their CM frame
//
// Here masses the two masses of the particles are all thats needed to specify the system
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef FOUR_MOMENTA_HPP
#define FOUR_MOMENTA_HPP

#include <string>
#include <complex>
#include <iostream>

#include "kinematics.hpp"

namespace jpacPhoto
{
    class raw_four_momenta;
    using four_momenta = std::shared_ptr<jpacPhoto::raw_four_momenta>;

    class raw_four_momenta
    {
        // -----------------------------------------------------------------------
        public:

        // Constructor
        raw_four_momenta(double m, double b)
        : _mM2(m*m), _mB2(b*b)
        {
            if (is_zero(m)) _photon = true;
        };
        
        // "Constructor" method to create a shared_ptr
        inline four_momenta create(double m, double b)
        {
            return std::make_shared<raw_four_momenta>(m, d);
        };

        // Whether or not our meson is massless photon
        inline bool if_photon(){ return _photon; };

        // Update masses that are saved
        inline void set_masses( std::array<double, 2> m)
        {
            _mM2 = m[0]*m[0]; _mB2 = m[1]*m[1]; 
            
            // If masses chance make sure to update the cache
            recalculate(_cached_s);
        };

        // Grab cached energies and momenta from outside
        inline double get_momenta(double s)
        {
            if (updated(s)) recalculate(s);
            return _momentum;
        };

        inline std::array<double,2> get_energies(double s)
        {
            if (updated(s)) recalculate(s);
            return {_meson_energy, _baryon_energy};
        };

        // 4-vector of vector, particle 1
        complex q(lorentz_index mu, double s, double theta); 

        // 4-vector of baryon, particle 2
        complex p(lorentz_index mu, double s, double theta); 

        // -----------------------------------------------------------------------
        private:

        bool _photon = false;
        double _mM2; // Meson mass
        double _mB2; // Baryon mass

        // Cache the energy dependence
        // Often we will be integrating over theta at fixes s, so this prevents recalculating
        // energy and momenta at each step
        bool updated(double s){ return is_equal(_cached_s, s); };
        void recalculate(double s);

        double  _cached_s = 0.;
        complex _momentum = 0.;
        complex _meson_energy = 0., _baryon_energy = 0.;
    };
};

#endif
