// Discontinuity across unitarity cut for vector production via a one-loop box diagram
// Used as a container class in integration processes
// 
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _BOX_DISC_
#define _BOX_DISC_

#include "constants.hpp"
#include "amplitudes/amplitude.hpp"
#include "reaction_kinematics.hpp"

#include "Math/IntegratorMultiDim.h"

namespace jpacPhoto
{
    class box_discontinuity
    {
        public: 
        box_discontinuity(double threshold)
        : _threshold(threshold)
        {};

        box_discontinuity(amplitude * left, amplitude * right)
        : _initialAmp(left), _finalAmp(right)
        {
            _threshold = left->_kinematics->sth();

            // Make sure the left and right amplitudes match!
            // Check the spins of the intermediate state
            _jp_left  = left->_kinematics->get_meson_JP();
            _jp_right = right->_kinematics->get_meson_JP();

            if (_jp_left != _jp_right) 
            {
                std::cout << std::left;
                std::cout << std::setw(40) << "box_amplitude: Intermediate state between sub-amplitudes dont match!" << std::endl;
                std::cout << std::setw(20) << _initialAmp->get_id() << ": \t (" << _jp_left[0]  << ", " << _jp_left[1]  << ")\n";
                std::cout << std::setw(20) << _finalAmp->get_id()   << ": \t (" << _jp_right[0] << ", " << _jp_right[1] << ")\n";
                std::cout << "Returning 0!" << std::endl;

                _matchError = true;
            };

            // But also masses
            if ((std::abs(left->_kinematics->get_meson_mass() - right->_kinematics->get_meson_mass()) > 1.E-4) 
             || (std::abs(left->_kinematics->get_recoil_mass() - right->_kinematics->get_recoil_mass()) > 1.E-4))
            {
                _matchError = true;
            };
            
            // IF they match, get the spin and therefor helicities of the intermediate meson
            _intermediate_helicities = get_helicities(_jp_left[0], _initialAmp->_kinematics->get_baryon_JP()[0]);
        };

        // Evaluate the discontinuity integrated over intermediate phase space
        virtual double eval(double s);

        // Pass values that come from the external gamma p -> V p reaction
        inline void set_externals(std::array<int,4> helicities, double theta)
        {
            _external_theta = theta; 
            _external_helicities = helicities;
        };

        double _threshold; 

        private:
        double _external_theta;
        std::array<int,4> _external_helicities;

        // Individual tree amplitudes
        amplitude * _initialAmp;
        amplitude * _finalAmp;

        std::array<int,2> _jp_left, _jp_right;
        std::vector< std::array<int,4> > _intermediate_helicities;
        bool _matchError = false;
    };

    class test_disc : public box_discontinuity
    {
        public: 

        test_disc(double w)
        : box_discontinuity(w*w)
        {};

        inline double eval(double s)
        {
            if (s < _threshold) return 0.;
            return sqrt(s - _threshold);
        };
    };
};

#endif