// Reggeized exchange of an unnatural parity particle
// Adapted from the model considered in [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------
// References:
// [1] - https://arxiv.org/abs/1710.09394
// ------------------------------------------------------------------------------

#ifndef UNNATURAL_EXCHANGE_HPP
#define UNNATURAL_EXCHANGE_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto
{
    namespace regge
    {
        class unnatural_exchange : public raw_amplitude
        {
            public:

            // Constructor
            unnatural_exchange(key k, kinematics xkinem, int signature, std::string id)
            : raw_amplitude(k, xkinem, id), _signature(signature)
            {  initialize(6); };

            // ---------------------------------------------------------------------------
            // Defining the virtual functions required of an amplitude
            
            inline complex helicity_amplitude(std::array<int, 4> helicities, double s, double t)
            {
                // Save inputs
                store(helicities, s, t);
                return half_angle()*form_factor()*top()*bottom()*propagator();
            };

            // Explicitly require t-channel helicities
            inline helicity_frame native_helicity_frame()        { return   S_CHANNEL;      };
            inline std::vector<quantum_numbers> allowed_mesons() { return { PSEUDOSCALAR }; };
            inline std::vector<quantum_numbers> allowed_baryons(){ return { THREEPLUS    }; };

            // -----------------------------------------------------------------------
            // Internal data members 

            protected:

            double _s0 = 1;
            // Set parameters
            inline void allocate_parameters(std::vector<double> x)
            {
                _alpha0 = x[0]; _alphaP = x[1];
                _gT     = x[2];
                _gB     = x[3];
                _b      = x[4], _c      = x[5];
                return;
            };

            // Regge trajectory
            int _signature = +1;
            double _alpha0 = 0, _alphaP = 0;
            inline double trajectory(){ return _alpha0 + _alphaP * _t; };
            inline complex propagator() // includes ghost killing factor
            {
                double alpha = trajectory();
                complex sig_factor  = (_signature + exp(-I*PI*alpha));
                double ghost_factor = (_signature == +1) ? (alpha+2)/2. : (alpha+1);
                return (PI*_alphaP)/2.* sig_factor*ghost_factor/sin(PI*alpha) * pow(_s/_s0, alpha);
            };

            // Top coupling
            double _gT = 0.;
            double top(){ return (_signature == +1) ? _gT : sqrt(-_t)*_gT; };

            // Bottom coupling
            double _gB = 0.;
            double bottom()
            {
                int phase = 1, sign = 1;
                if (_lamT < 0)
                {
                    phase = pow(-1, (_lamT - _lamR)/ 2); sign = -1;
                }

                int n = std::abs(-_lamB - (_lamR - _lamT)/2); 
                int x = std::abs((_lamT - _lamR)/2) + abs(_lamB) - n;
                double tfactor = phase * pow(sqrt(-_t), n) * pow(-M_PION*M_PION, x/2);

                switch (sign*_lamR)
                {
                    case  3: return   tfactor*_gB*(_mT + _mR)/(sqrt(2.)*_mR);
                    case -1: return - tfactor*_gB*(-_mT*_mT + _mT*_mR + 2.*_mR*_mR + _t)/(sqrt(6)*_mR*_mR);
                    case  1: return - tfactor*_gB*(-_mT*_mT*_mT - _mT*_mT*_mR + _mR*_mR*_mR + 2.*_mR*_t + _mT*(_mR*_mR + _t))/(sqrt(6)*_mR*_mR);
                    case -3: return - tfactor*_gB/(sqrt(2.)*_mR);
                    default:  return std::nan("");
                }
                return std::nan("");
            };

            // Suppression 
            double _b = 0., _c = 0.;
            double form_factor(){ return _c*exp(_b*_t); };

            // Half ange factor
            complex half_angle()
            {
                double mui = double(_lamB) - _lamT/2., muf = double(_lamX) - _lamR/2.;
                double z = cos(_theta);
                return pow((_s/-_t)*(1-z)/2, std::abs(mui-muf)/2) * pow( (1+z)/2, std::abs(mui + muf)/2);
            };

        };
    };
};

#endif