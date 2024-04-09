// Reproduces a semi-incusive distribution following phasespace with constant amplitude
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#ifndef PHASESPACE_HPP       
#define PHASESPACE_HPP

#include "constants.hpp"
#include "semi_inclusive.hpp"
#include "kinematics.hpp"

namespace jpacPhoto
{
    namespace inclusive
    {
        class phase_space : public raw_semi_inclusive
        {
            public: 

            phase_space(key k, kinematics mX, double mMin, std::string id = "")
            : raw_semi_inclusive(k, mX, id), _mMin2(mMin*mMin)
            {
                set_N_pars(0);
            };

            // Minimum mass is the proton 
            double minimum_M2(){ return _mMin2; };

            // Amplitude is a constant
            double invariant_xsection(double s, double t, double mm){ return 1; };

            private:

            // Set the minimum mass
            double _mMin2 = (M_PROTON + M_PION)*(M_PROTON + M_PION);
        };
    };
};

#endif