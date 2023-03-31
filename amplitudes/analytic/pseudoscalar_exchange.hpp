// Amplitude defining production via a spin-0 exchange in the t-channel
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef ANALYTIC_PSEUDOSCALAR_EXCHANGE_HPP
#define ANALYTIC_PSEUDOSCALAR_EXCHANGE_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "form_factor.hpp"
#include "amplitude_options.hpp"
#include "amplitude.hpp"

namespace jpacPhoto
{
    namespace analytic
    {
        class pseudoscalar_exchange : public raw_amplitude
        {
            public:

            // Constructor
            pseudoscalar_exchange(amplitude_key key, kinematics xkinem, double exchange_mass, std::string id)
            : raw_amplitude(key, xkinem, id),
              _mEx(exchange_mass)
            {
                initialize(3);
            };

            // ---------------------------------------------------------------------------
            // Defining the virtual functions required of an amplitude
            
            inline complex helicity_amplitude(std::array<int, 4> helicities, double s, double t)
            {
                // Save inputs
                store(helicities, s, t);
                _qi = _kinematics->initial_momentum_tframe(_t);
                _qf = _kinematics->final_momentum_tframe(_t);

                // An on-shell spin-0 exchange cannot flip helicities
                bool helicity_conserving = (_lamB == _lamX) && (_lamR == _lamT);
                if (!helicity_conserving) return 0;
                
                // Parse which argument should go into the form-factor
                // The exponential takes t' = t - tmin while monopole takes just t
                complex FF = (_option == amplitude_option::ExpFF) ? _FF->eval(_t - _kinematics->t_min(s))
                                                                  : _FF->eval(_t);

                // Multiply couplings with propagator
                return FF * top_coupling() * propagator() * bottom_coupling();
            };

            // Explicitly require t-channel helicities
            inline helicity_channel native_helicity_frame(){ return helicity_channel::T_CHANNEL; };

            // We can have pseudo-scalar, vector, and axial-vector
            inline std::vector<std::array<int,2>> allowed_meson_JP()
            {
                return { PSUEDOSCALAR, VECTOR, AXIALVECTOR };
            };

            // But only either parity spin-1/2
            inline std::vector<std::array<int,2>> allowed_baryon_JP()
            {
                return { HALFPLUS, HALFMINUS };
            };

            // The options here are the type of form_factor used
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

            // Initial & final momentum in the t-channel frame
            complex _qi = 0, _qf = 0;

            // We include a t-channel form-factor for the propagator
            // By default this is the exponential one
            form_factor _FF = new_FF<exponential>();

            // Top coupling refers to the beam-exchange-meson interaction
            inline complex top_coupling()
            {
                // Coupling function depends on
                // the quantum numbers of the produced meson
                auto JP  = _kinematics->get_meson_JP();
                int  qns =  10 * JP[0] + (JP[1] == 1);

                complex result = 0;
                switch (qns)
                {
                    // Axial-vector
                    case (11): result = 1/_mX; break;
                    
                    // Vector
                    case (10): result = (_kinematics->is_photon()) ? -4 : -1; break;

                    // Pseudo-scalar
                    case ( 0): result = (_kinematics->is_photon()) ? 0 : 2*I/_mB; break;

                    default: break;
                };

                return _gTop*_qi*csqrt(_t)*result;
            };

            // Bottom coupling refers to the target-exchange-recoil interation
            inline complex bottom_coupling()
            {
                // This is also dependent on baryon quantum numbers
                auto JP = _kinematics->get_baryon_JP();
                int qns = 10 * JP[0] + (JP[1] == 1);

                complex result = 0;
                switch (qns)
                {
                    // 1/2+
                    case (11) : result = csqrt( _t - (_mT-_mR)*(_mT-_mR)); break;
                    
                    // 3/2+
                    case (31) : result = csqrt(2*_t/3)*(_qf/_mR)*csqrt( _t - (_mT+_mR)*(_mT+_mR))/2; break;

                    default: break;
                };

                return _gBot * result;
            };

            // Scalar propagator in the t-channel
            inline complex propagator()
            {
                return - I / (_t - _mEx*_mEx);
            };
            
        };
    };
};

#endif