// Generic reggeon exchange model with simple vertices
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

#ifndef REGGEON_EXCHANGE_HPP
#define REGGEON_EXCHANGE_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto
{
    namespace regge
    {
        class reggeon_exchange : public raw_amplitude
        {
            public:

            // Constructor
            reggeon_exchange(key k, kinematics xkinem, int naturality, std::string id)
            : raw_amplitude(k, xkinem, id), _naturality(naturality)
            {  
                if (abs(naturality) > 1) warning("reggeon_exchange", "Invalid naturality passed to constructor!");
                if (naturality == -1) initialize(4); 
                else                  initialize(7);
            };

            // ---------------------------------------------------------------------------
            // Defining the virtual functions required of an amplitude
            
            inline complex helicity_amplitude(std::array<int, 4> helicities, double s, double t)
            {
                // Save inputs
                store(helicities, s, t);
                return _beta*exp(_b*t)*top()*propagator()*bottom();
            };

            // Explicitly require t-channel helicities
            inline helicity_frame native_helicity_frame()        { return   S_CHANNEL;  };
            inline std::vector<quantum_numbers> allowed_mesons() { return {  VECTOR };  };
            inline std::vector<quantum_numbers> allowed_baryons(){ return { HALFPLUS }; };

            // -----------------------------------------------------------------------
            // Internal data members 

            protected:

            // Set parameters
            inline void allocate_parameters(std::vector<double> x)
            {
                _alpha0 = x[0];
                _alphaP = x[1];
                _beta   = x[2]; 
                _b      = x[3];
                if (_naturality == -1)
                {
                    _gT[0] = 0.;   _gT[1] = 0.;   _kappa = 0.;
                }
                else
                {
                    _gT[0] = x[4]; _gT[1] = x[5]; _kappa = x[6];         
                };
            };

            //----------------------------------------------------
            int _signature  = +1; // Paper explicitly fixes all trajectories to + signature
            int _naturality = +1;
            double _s0 = 1;

            // Couplings
            double _beta = 0, _b = 0., _kappa = 0;
            std::array<double,2> _gT;

            // Linear trajectory
            double _alpha0 = 0, _alphaP = 0;
            //----------------------------------------------------

            inline double trajectory(){ return _alpha0 + _alphaP * _t; };
            
            // Regge propagator 
            inline complex propagator()
            {
                double alpha = trajectory();
                double gf = (_naturality == -1) ? 1. : alpha / _alpha0; // includes ghost killing factor for natural
                return _alphaP*M_PI*(_signature + exp(-I*PI*alpha))/2.*gf/sin(PI*alpha)*pow(_s/_s0, alpha);
            };

            double top()
            { 
                double sqt = sqrt(-_t)/_mX;
                int g = _lamB, V = _lamX;

                if (_naturality == -1) return g*(g == V) - sqrt(2.)*sqt*(V == 0) + g*sqt*sqt*(g == -V);

                return ((g == V) + _gT[0]*sqt*g/sqrt(2)*(V == 0) + _gT[1]*sqt*sqt*(g == -V));
            };

            double bottom()
            {
                double sqt = sqrt(-_t)/(_mT+_mR);
                int l = _lamT, lp = _lamR;

                if (_naturality == -1) return - sqt*(l == -lp);
                return ((l == lp) + 2.*_kappa*sqt*(l == -lp));
            };
        };
    };
};

#endif