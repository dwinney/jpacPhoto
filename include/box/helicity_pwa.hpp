// Container class that calculates a specific helicity partial-wave projection
// of a given amplitude
// So all helicities and total spin J are assumed fixed.
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef HWA_AMP
#define HWA_AMP

#include "constants.hpp"
#include "amplitude.hpp"
#include "reaction_kinematics.hpp"

#include <Math/Interpolator.h>

#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace jpacPhoto
{
    class helicity_pwa
    {
        // ---------------------------------------------------------------------------
        public: 

        // Constructor
        // Needs pointer to an amplitude object which we will take projections of
        // Also the spins and helicities
        helicity_pwa(amplitude * amp, int J, std::array<int,4> helicities)
        : _amplitude(amp), _kinematics(amp->_kinematics), 
          _J(J), _helicities(helicities)
        {};

        // Destructor should clean up interpolation if that was initiated
        ~helicity_pwa()
        {
            if (_interpSaved) delete _interp;
        };

        // Output the helicity parial wave projection onto the s-channel with total spin j/2
        double helicity_partial_wave(double s);

        // Update and output the helicity PWA from saved interpolation
        void   update_interpolation();
        double interpolation(double s);

        // Access private members
        inline int get_J(){ return _J; };
        inline std::array<int, 4> get_helicities(){ return _helicities; };

        // Set the interval of interpolation
        inline void set_interpolation_cutoff(double s){ _smax = s; update_interpolation(); };

        // ---------------------------------------------------------------------------
        private: 

        // Saved amplitude which we are projecting
        reaction_kinematics * _kinematics;
        amplitude * _amplitude;

        // Partial-wave projection spin and helicities
        int _J;
        std::array<int,4> _helicities;

        // Save an interpolation to make subsequent evaluation faster 
        bool   _interpSaved = false;
        int    _Ninterp = 50;  // Number of interpolation points
        double _smin = (_kinematics->sth() + EPS);
        double _smax = -1;       
        std::vector<double> _x, _fx;
        ROOT::Math::Interpolator * _interp;
    };
};

#endif