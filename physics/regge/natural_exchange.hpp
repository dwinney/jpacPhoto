// Reggeized exchange of an natural parity particle
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

#ifndef NATURAL_EXCHANGE_HPP
#define NATURAL_EXCHANGE_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto
{
    namespace regge
    {
        class natural_exchange : public raw_amplitude
        {
            public:

            // Constructor
            natural_exchange(key k, kinematics xkinem, int signature, std::string id)
            : raw_amplitude(k, xkinem, id), _signature(signature)
            {  initialize(7); };

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

            // Options
            static const int kPole = 0;
            static const int kCut  = 1;

            // -----------------------------------------------------------------------
            // Internal data members 

            protected:
            
            // Fixed parameters
            int _signature = +1;
            double _s0 = 1;
            double _alpha0Pom = 1.08, _alphaPPom = 0.25; // Pomeron trajectory for cut model

            // Free parameters
            double _alpha0 = 0, _alphaP = 0;
            double _gT = 0., _gB1 = 0., _gB2 = 0, _gB3 = 0;
            double _b = 0.;

            // Set parameters
            inline void allocate_parameters(std::vector<double> x)
            {
                _alpha0 = x[0]; _alphaP = x[1];
                _gT     = x[2];
                _gB1    = x[3], _gB2    = x[4], _gB3 = x[5];
                _b      = x[6];
                return;
            };

            // Regge trajectory
            inline double slope()     { return (_option == kCut) ? _alphaPPom * _alphaP / (_alphaPPom + _alphaP) : _alphaP; };
            inline double trajectory(){ return (_option == kCut) ? _alpha0Pom - _alphaP*M_RHO*M_RHO + slope()*_t : _alpha0 + slope()*_t; };
            inline complex propagator() // includes ghost killing factor
            {
                double alpha = trajectory();
                complex sig_factor  = (_signature + exp(-I*PI*alpha));
                double ghost_factor = (_signature == +1) ? alpha*(alpha+2)/3 : (alpha+1)/2;
                if (_option == kCut) ghost_factor /= log(_s);

                return (PI*slope())/2 * sig_factor*ghost_factor/sin(PI*alpha) * pow(_s/_s0, alpha);
            };

            // Top coupling
            double top(){ return sqrt(-_t)*_gT; };

            // Bottom coupling
            double bottom()
            {
                int phase = 1, sign = 1;
                if (_lamT < 0)
                {
                    phase = pow(-1, 1 + (_lamT - _lamR)/ 2); sign = -1;
                }

                double tfactor = phase * pow(sqrt(-_t), abs(_lamT - _lamR)/2);

                switch (sign*_lamR)
                {
                    case  3: return - tfactor*(2*_mR    *_gB1 + _gB2*(_mT - _mR))                                     /(2.*_mR*_mR);
                    case -1: return - tfactor*(2*_mT*_mR*_gB1 + _gB2*(-_mT*_mR + _mR*_mR + 2*_t) + 2*_t        *_gB3) /(2*sqrt(3.)*_mR*_mR*_mR);
                    case  1: return - tfactor*(2.*_mR   *_gB1 + _gB2*(2.*_mT-3.*_mR)             + 2.*(_mT-_mR)*_gB3) /(2*sqrt(3.)*_mR*_mR*_mR)*(-_t);
                    case -3: return - tfactor*_gB2/(2.*_mR*_mR);
                    default:  return std::nan("");
                }
                return std::nan("");
            };

            // Suppression factor
            double form_factor(){ return exp(_b*_t); };

            // Half ange factor
            complex half_angle()
            {
                double mui = double(_lamB) - _lamT/2., muf = double(_lamX) - _lamR/2.;
                double z = cos(_theta);
                return pow((_s/-_t)*(1-z)/2, std::abs(mui-muf)/2.) * pow( (1+z)/2, std::abs(mui + muf)/2.);
            };

        };
    };
};

#endif