// Implementation of a simple pomeron exchange assuming spin-less external particles
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
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
            : raw_amplitude(key, xkinem, "pomeron_exchange", id)
            {
                set_N_pars(4);
                check_QNs(xkinem);
            };

            // -----------------------------------------------------------------------
            // Virtuals 

            // We assume spinless particles so we have no helicity structure
            inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
            {
                // Arbitrarily pick one of the helicities to evaluate
                if (helicities != _kinematics->helicities(0)) return 0;

                // Save inputes
                store(helicities, s, t);

                // Intermediate quantities
                double t_prime = t - _kinematics->t_min(s);
                double alpha   = _alpha0 + t * _alphaP;

                // Normalization here to get rid of helicity dependence in amplitude::probability_distribution
                // First a 2 removes the factor 1/4 when averaging over initial helicities
                // then a 1/sqrt(2) removes the factor of 2 from the parity relation in amplitude::update_cache
                return sqrt(2) * _A * exp(_b0 * t_prime) * pow((s - _kinematics->sth()) / _s0, alpha);
            };

            // Explicitly require t-channel helicities
            inline helicity_channel native_helicity_frame(){ return helicity_channel::S_CHANNEL; };

            // We can have any quantum numbers
            // but for now explicitly put only pseudo-scalar, vector, and axial-vector
            // and either parity spin-1/2
            inline std::vector<std::array<int,2>> allowed_meson_JP() { return { PSUEDOSCALAR, VECTOR, AXIALVECTOR }; };
            inline std::vector<std::array<int,2>> allowed_baryon_JP(){ return { HALFPLUS, HALFMINUS }; };

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

            // -----------------------------------------------------------------------
            private:

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