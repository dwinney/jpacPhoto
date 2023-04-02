// Amplitude defining production via a spin-0 exchange in the t-channel,
// this is evaluated via contraction of covariant quantities in the s-channel
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef COVARIANT_PSEUDOSCALAR_EXCHANGE_HPP
#define COVARIANT_PSEUDOSCALAR_EXCHANGE_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "form_factor.hpp"
#include "amplitude_options.hpp"
#include "amplitude.hpp"

namespace jpacPhoto
{
    namespace covariant
    {
        class pseudoscalar_exchange : public raw_amplitude
        {
            public:

            // Constructor we specify the exchange particle mass
            pseudoscalar_exchange(amplitude_key key, kinematics xkinem, double exchange_mass, std::string id = "pseudoscalar_exchange")
            : raw_amplitude(key, xkinem, id), _mEx(exchange_mass)
            {
                initialize(3);
            }
            
            // ---------------------------------------------------------------------------
            // VIRTUALS

            inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
            {
                // Save inputs
                store(helicities, s, t);
                _covariants->update(helicities, s, t);
                
                // Parse which argument should go into the form-factor
                // The exponential takes t' = t - tmin while monopole takes just t
                complex FF = (_option == amplitude_option::ExpFF) ? _FF->eval(_t - _kinematics->t_min(s))
                                                                  : _FF->eval(_t);
                
                // Multiply couplings with propagator
                return FF * top_coupling() * propagator() * bottom_coupling();
            }

            // Covariants are s-channel amplitudes
            inline helicity_frame native_helicity_frame(){return S_CHANNEL; };

            // We can have pseudo-scalar, vector, and axial-vector
            inline std::vector<particle> allowed_mesons()
            {
                return { PSEUDOSCALAR, VECTOR, AXIALVECTOR };
            };

            inline std::vector<std::array<int,2>> allowed_baryons()
            {
                return { HALFPLUS, THREEPLUS };
            };
            
            // The options here are the type of form_factor used
            // Default assumed exponential
            inline void set_option( amplitude_option opt )
            {
                switch (opt)
                {
                    case (Default):
                    case (ExpFF): 
                    {
                        _FF = new_FF<exponential>();
                        _option = ExpFF;
                        break;
                    };
                    case (MonopoleFF):
                    {
                        _FF = new_FF<monopole>(_mEx);
                        _option = MonopoleFF;
                        break;
                    }
                    default: 
                    {
                        option_error();
                        return;
                    };
                };

                // Before leaving make sure the new FF gets the last saved cutoff
                _FF->set_cutoff(_ffCutoff);
            };
            
            // Parameter names
            inline std::vector<std::string> parameter_labels()
            {
                return { "gPhoton", "gNucleon", "Cutoff" };
            };
            
            // -----------------------------------------------------------------------
            // Internal data members 

            protected:

            // Have three free parameters: 
            // [0] Top (beam-exchange-meson) coupling
            // [1] Bottom (target-exchange-recoil) coupling
            // [2] Form-factor cutoff
            inline void allocate_parameters(std::vector<double> x)
            {
                _gTop     = x[0];
                _gBot     = x[1];
                _ffCutoff = x[2];

                // Pass cutoff to FF as well
                _FF->set_cutoff(_ffCutoff);
                return;
            };

            // Exchange mass
            double _mEx;

            // Free parameters: two couplings and a form-factor cutoff
            double _gTop = 0, _gBot = 0, _ffCutoff = 0;

            // We include a t-channel form-factor for the propagator
            // By default this is the exponential one
            form_factor _FF = new_FF<exponential>();

            // COVARIANT PIECES

            // Scalar propagator in the t-channel
            inline complex propagator()
            {
                return - I / (_t - _mEx*_mEx);
            };

            inline complex top_coupling()
            {
                // Beam
                auto q     = _covariants->q();
                auto eps   = _covariants->eps();
                
                // Outgoing meson
                auto q_p   = _covariants->q_prime();
                auto eps_p = _covariants->eps_prime();

                complex result = 0;

                switch (_kinematics->get_meson())
                {
                    case (AXIALVECTOR): { result = contract(eps,   eps_p) * contract(q,   q_p) 
                                                 - contract(eps, q_p)     * contract(eps_p, q);
                                          result /= _mX; break;};                                          
                    
                    case (VECTOR): { result = 4*levi_civita(eps_p, eps, q, q_p);
                                     if (!_kinematics->is_photon()) result /= 4;  break;}

                    case (PSEUDOSCALAR): { result = - contract(eps, q - 2 * q_p); break;}

                    default: break;
                };

                return _gTop * result; 
            };

            inline complex bottom_coupling()
            {
                auto u    = _covariants->u();    // Target
                auto ubar = _covariants->ubar(); // Recoil

                complex result = contract(ubar, gamma_5() * u);
            
                return _gBot * result;
            };

        };
    };
};

#endif
