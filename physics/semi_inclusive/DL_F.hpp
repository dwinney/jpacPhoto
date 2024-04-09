// Implementation of unpolarized structure functions at high energies
// Uses the Regge pole parameteriztion from [1]
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
// ------------------------------------------------------------------------------

#ifndef DL_F_HPP
#define DL_F_HPP

#include "constants.hpp"
#include "data_set.hpp"
#include "inclusive_function.hpp"
#include <Math/Interpolator.h>

namespace jpacPhoto
{
    class DL_F : public raw_inclusive_function
    {
        public: 

        DL_F(unsigned x)
        : raw_inclusive_function({0, M_PROTON}) // Massless photon
        {
            if (x != 1 && x != 2)
            {
                error("DL_F", "Integer argument must be 1 or 2 for F_1 and F_2 respectively! Defaulting to F_2...");
                return;
            }

            // If F_1 requested, import R values
            if (x == 1) _RBS.SetData(_interpQs, _interpRs);
            _mode = x;
        };

        inline double evaluate(double s, double q2)
        {
            // Save external variables
            _s = s; _w = sqrt(s); _Q2 = -q2;
            _nu = nu(); _x = x();

            // Calculate residues for each exchange
            for (int i = 0; i < 3; i++)
            {
                _f[i] = _A[i] * pow(_Q2 / (1.+ _Q2/_Q20[i]), 1.+_eps[i]);    
                if (i == 0) _f[i] *= pow(1.+_Q2/_Q20[i], _eps[i]/2.);        
            };

            // Sum terms with the appropriate Regge propagator
            double F = 0;
            for (int i = 0; i < 3; i++)
            {
                int del = (i == 2) ? 1 : 5;
                double exponent = (_mode == 2) ? -_eps[i] : -(1+_eps[i]);
                F += _f[i] * pow(_x, exponent) * pow(1.-_x, del);
            };

            if (_mode == 1) F /= 2*(1+R());

            return F;
        };

        inline double R(double Q2)
        {
            return (Q2 < 1.5) ? _RBS.Eval(Q2) : _R0;
        };
        inline double R(){ return R(_Q2); };

        private:

        int _mode = 2;

        // Saved external variables
        double _s, _w, _Q2;

        // Saved internal variables
        double _nu, _x;

        // Residues
        std::array<double,3> _f;

        // R = sigma_L / sigma_T parameters
        double _R0 = 0.23;
        // Instead of calcualting R from Christy & Bolsted,
        // Digitize and interpolate Fig. 15 of [1]
        ROOT::Math::Interpolator _RBS;

        // parameters
        std::array<double,3> _A   = {0.00151, 0.658,   1.01};  // normalizations
        std::array<double,3> _eps = {0.452,   0.0667, -0.476}; // intercepts
        std::array<double,3> _Q20 = {7.85,    0.6,     0.214}; // scale parameters

        // Energy variables
        inline double nu(){ return (_s + _Q2 - M2_PROTON) / (2*M_PROTON); };
        inline double  x(){ return _Q2 / (_s + _Q2 - M2_PROTON); };

        // Interpolate the plot for R(Q2) instead of evaluating the ratio explicitly
        std::vector<double> _interpQs = {
            0.0,
            0.01699641458,
            0.08440235397,
            0.1363918851,
            0.2574187985,
            0.4021573047,
            0.5311381068,
            0.6737513805,
            0.7905496135,
            0.943401764,
            1.149486392,
            1.276574333,
            1.453479372,
            1.525638039
        };
        
        std::vector<double> _interpRs = {
            0.0,
            0.005070422535,
            0.06422535211,
            0.1312676056,
            0.243943662,
            0.3025352113,
            0.316056338,
            0.3115492958,
            0.3014084507,
            0.2856338028,
            0.2630985915,
            0.2495774648,
            0.2326760563,
            0.2292957746
        };
    };
};

#endif