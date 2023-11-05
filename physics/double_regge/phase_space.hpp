//  A two-meson photoproduction "amplitude" that just returns 1
//
// --------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
//               University of Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// --------------------------------------------------------------------------------

#ifndef DR_PHASESPACE_HPP
#define DR_PHASESPACE_HPP

#include "constants.hpp"
#include "kinematics2.hpp"
#include "amplitude2.hpp"

namespace jpacPhoto
{
    class phase_space : public two_meson::raw_amplitude
    {
        public: 

        // Basic constructor
        phase_space(two_meson::amplitude_key key, two_meson::kinematics xkinem, std::string id = "phase_space")
        : two_meson::raw_amplitude(key, xkinem, id)
        {
            initialize(0);
        };

        inline complex helicity_amplitude(std::array<int,3> helicities, double s, double t, double s12, double thetaGJ, double phiGJ)
        {
            return 1.;
        };
    };
};

#endif