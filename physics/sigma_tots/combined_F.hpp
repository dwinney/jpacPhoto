// Implementation of unpolarized structure functions at high energies
// combines the low energy resonance region from [1]
// and the asymptotic regge region in [2]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] - https://arxiv.org/abs/hep-ph/0402081
// [2] - https://arxiv.org/abs/0712.3731
// ------------------------------------------------------------------------------

#ifndef COMBINED_F_HPP
#define COMBINED_F_HPP

#include "constants.hpp"
#include "data_set.hpp"
#include "inclusive_function.hpp"
#include "sigma_tots/CB_F.hpp"
#include "sigma_tots/DL_F.hpp"

#include <Math/Interpolator.h>

namespace jpacPhoto
{
    class combined_F : public raw_inclusive_function
    {
        public: 

        combined_F(unsigned x, int p_or_n)
        : raw_inclusive_function({0, M_PROTON}),
          _low( new_inclusive_function<CB_F>(x, p_or_n)),
          _high(new_inclusive_function<DL_F>(x))
        {
            if (x != 1 && x != 2)
            {
                error("combined_F", "Integer argument must be 1 or 2 for F_1 and F_2 respectively! Defaulting to F_2...");
                return;
            }

            _mode = x;
        };

        inline double evaluate(double s, double q2)
        {
            if (sqrt(s) <= _low_W)
            {
                return _low->evaluate(s, q2);
            }
            if (sqrt(s) >=_high_W) 
            {
                return _high->evaluate(s, q2);
            }   
            
            double low_F  = _low->evaluate( _low_W*_low_W,   q2);
            double high_F = _high->evaluate(_high_W*_high_W, q2);

            return low_F + (s - _low_W*_low_W)/(_high_W*_high_W - _low_W*_low_W) * (high_F - low_F);
        };

        // Flags to denote which cross-sections we're looking for
        static const int kProton = 0, kNeutron = 1;

        private:

        int _mode = 2;

        inclusive_function _low, _high;

        double _low_W  = 2.5;
        double _high_W = 5;
        double _Lam    = 0.5;

        // Saved external variables
        double _s, _w, _Q2;
    };
};

#endif