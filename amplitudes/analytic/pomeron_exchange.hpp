// Implementation of a simple pomeron exchange assuming perfect helicity conservation
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef ANALYTIC_POMERON_HPP
#define ANALYTIC_POMERON_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto
{
    namespace analytic
    {
        class pomeron_exchange : public raw_amplitude
        {
            public: 

            pomeron_exchange(amplitude_key key, kinematics xkinem, std::string id = "pomeron_exchange")
            : raw_amplitude(key, xkinem, id)
            {
                initialize(4);
            };

            // -----------------------------------------------------------------------
            // Virtuals 

            // We assume spinless particles so we have no helicity structure
            inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
            {
                // Save inputes
                store(helicities, s, t);

                // Intermediate quantities
                double t_prime = t - _kinematics->t_min(s);
                double alpha   = _alpha0 + t * _alphaP;

                // Helicity conservation
                bool helicity_conserving = (_lamB == _lamX) && (_lamT == _lamR);

                return _A * exp(_b0 * t_prime) * pow((s - _kinematics->sth()) / _s0, alpha) * helicity_conserving;
            };

            // Even though its "analytic" we require s-channel helicity conservation
            inline helicity_channel native_helicity_frame(){ return helicity_channel::S_CHANNEL; };

            // Vector mesons and half plus only
            inline std::vector<particle> allowed_mesons() { return { vector }; };
            inline std::vector<particle> allowed_baryons(){ return { halfplus }; };

            // Parameter names are a[J] and b[J] for scattering length and normalization respectively
            inline std::vector<std::string> parameter_labels()
            {
                return { "A", "b0", "alpha_0", "alpha_prime" }; 
            };
            
            protected:

            inline void allocate_parameters(std::vector<double> pars)
            {
                _A      = pars[0];
                _b0     = pars[1];
                _alpha0 = pars[2];
                _alphaP = pars[3];
            };

            // Free parameters
            double _A      = 1; // Overall noramalization
            double _b0     = 0; // constant t-slope parameter [GeV-2]
            double _alpha0 = 0; // Trajectory intercept
            double _alphaP = 0; // Trajectory slope [GeV-2]

            // Fixed scale parameter
            double _s0     = 1; // GeV2
        };
    };
};

#endif