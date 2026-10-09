// Implementation of a simple pomeron exchange 
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef POMERON_HPP
#define POMERON_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto
{
    namespace regge 
    {
        class pomeron_exchange : public raw_amplitude
        {
            public: 

            pomeron_exchange(key k, kinematics xkinem)
            : raw_amplitude(k, xkinem, "pomeron_exchange")
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

                // Helicity structure comes from contractin the top and bottom vertices
                // We remove one factor of s to account for the asymptotic scaling of the dirac_spinors
                complex helicity_structure;
                if (_option == kVecPom)
                {
                    _covariants->update(helicities, s,t);
                    helicity_structure = contract(top_vertex(), bottom_vertex()) / _s;
                }
                else helicity_structure = (_lamB==_lamX)*(_lamT==_lamR);

                // Nothing here depends on helicities
                return _A*exp(_b0*t_prime)*pow((s - _kinematics->sth())/_s0,alpha)*helicity_structure;
            };

            // Even though its "analytic" we require s-channel helicity conservation
            inline helicity_frame native_helicity_frame(){ return helicity_frame::S_CHANNEL; };

            // Vector mesons and half plus only
            inline std::vector<quantum_numbers> allowed_mesons() { return { VECTOR }; };
            inline std::vector<quantum_numbers> allowed_baryons(){ return { HALFPLUS }; };

            static const int kVecPom = 0;
            static const int kHelCon = 1;
            inline void set_option( int opt ){ _option = opt; };

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

            // Photon -- Pomeron -- Vector meson vertex
            inline lorentz_tensor<complex,1> top_vertex()
            {
                // Beam polarization and momentum
                auto eps   = _covariants->eps();
                auto q     = _covariants->q();

                // Vector polarization
                auto eps_p = _covariants->eps_prime();

                return contract(eps, eps_p) * q - contract(q, eps_p) * eps;
            };  

            // Nucleon -- Pomeron -- Nucleon vertex
            inline lorentz_tensor<complex,1> bottom_vertex()
            {
                // Spinors
                auto u    = _covariants->u();    // Target
                auto ubar = _covariants->ubar(); // Recoil

                return bilinear(ubar, gamma_vector(), u);;
            };
        };
    };
};

#endif