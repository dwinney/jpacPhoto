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
#include "cgamma.hpp"

namespace jpacPhoto
{
    namespace regge
    {
        class reggeon_exchange : public raw_amplitude
        {
            public:

            // Constructor
            reggeon_exchange(key k, kinematics xkinem, int naturality, std::string id)
            : raw_amplitude(k, xkinem, id), _naturality(naturality), _signature(+1)
            {  
                if (abs(naturality) > 1) warning("reggeon_exchange", "Invalid naturality passed to constructor!");

                _J = xkinem->get_meson_JP()[0];
                initialize(6 + 2*_J);
            };

            reggeon_exchange(key k, kinematics xkinem, int naturality, int signature, std::string id)
            : raw_amplitude(k, xkinem, id), _naturality(naturality), _signature(signature)
            {  
                if (abs(naturality) > 1) warning("reggeon_exchange", "Invalid naturality passed to constructor!");

                _J = xkinem->get_meson_JP()[0];
                initialize(6 + 2*_J);
            };

            // ---------------------------------------------------------------------------
            // Defining the virtual functions required of an amplitude
            
            inline complex helicity_amplitude(std::array<int, 4> helicities, double s, double t)
            {
                // Save inputs
                store(helicities, s, t);
                auto result = exp(_b*t)*top()*propagator()*bottom();
                return result;
            };

            // Explicitly require t-channel helicities
            inline helicity_frame native_helicity_frame()        { return   S_CHANNEL;  };
            inline std::vector<quantum_numbers> allowed_mesons() { return {  ANY };  };
            inline std::vector<quantum_numbers> allowed_baryons(){ return { HALFPLUS, HALFMINUS }; };

            static const int k2018_Model   = 0;
            static const int k2020_Minimal = 1;
            static const int k2020_TMD     = 2;
            static const int kRemoveZero   = 3;
            static const int kDefault = k2018_Model;

            inline void set_option(int x)
            {
                if (x > 3) option_error();
                if (x == kRemoveZero) { _remove_zero = true; return; };
                _option = x;
            };

            // -----------------------------------------------------------------------
            // Internal data members 

            protected:
            
            // Whether or not to remove the zero in the regge propagator at alpha = 0
            bool _remove_zero = false;

            // Set parameters
            inline void allocate_parameters(std::vector<double> x)
            {
                _alpha0 = x[0];
                _alphaP = x[1];
                _b      = x[2];

                _gT.clear();
                for (int i = 0; i < 1+2*_J; i++) _gT.push_back(x[3+i]);

                _gB[0]  = x[4+2*_J]; _gB[1] = x[5+2*_J];   
            };

            //----------------------------------------------------
            int _signature  = +1; // Paper explicitly fixes all trajectories to + signature
            int _naturality = +1;
            double _s0 = 1;
            int _J = 0;

            // Couplings
            double _b = 0.;           // Form factor cutoff
            std::vector<double> _gT;  // non-flip, single flip, double flip
            std::array<double,2> _gB; // non-flip and flip

            // Linear trajectory
            double _alpha0 = 0, _alphaP = 0;

            //----------------------------------------------------

            inline double trajectory(){ return _alpha0 + _alphaP * _t; };
            
            // Regge propagator 
            inline complex propagator()
            {
                double alpha = _alpha0 + _alphaP * _t;
                complex sigf = (_signature + exp(-I*PI*alpha))/2.;

                if (_option == k2018_Model)
                {
                    // Ghost killing factor for natural exchanges
                    double gf = (_naturality == +1) ? alpha / _alpha0 : 1.; 
                    if (is_zero(gf)) gf = 1/_alpha0/PI; // Treat the pole at alpha(t) = 0 special
                    else gf /= sin(PI*alpha);
                    return -_alphaP*sigf*gf*PI*pow(_s/_s0, alpha);
                };

                bool if_remove_zero = (_remove_zero && _lamT == -_lamR);
                int ell = (_naturality < 0 || if_remove_zero) ?  0 : 1;
                return pow(-1, 1+ell)*_naturality*sigf*cgamma(ell-alpha)*pow(_alphaP*_s, alpha);               
            };

            inline double top()
            { 
                int li    = abs(_lamB - _lamX);
                int index = _J - _lamX; // _gT's are indexed from [0, 1, .., 2J-1, 2J] = [J, J-1, ..., -J+1, -J]
                int phase = (_option == k2018_Model) ? pow(-_lamB, li + (_naturality < 0)) : pow(_lamB, li + (_naturality < 0));;
                double t  = (_option == k2018_Model) ? _t : _t - _kinematics->t_min(_s);
                return phase*_gT[index]*pow( sqrt(-t)/_mX, li)*beta_T(index);
            };

            inline double beta_T(int i)
            {
                if (_option == k2020_TMD && ( i == 1 || i == 2)) return -_t/2/(_mX*_mX);
                return 1.;
            };

            inline double bottom()
            {
                int lf    = abs(_lamT - _lamR)/2;
                int phase = (_option == k2018_Model) ? pow(-_lamT, lf + (_naturality < 0)) : pow(_lamT, lf + (_naturality < 0));
                double t  = (_option == k2018_Model) ? _t : _t - _kinematics->t_min(_s);
                return phase*_gB[lf] * pow( sqrt(-t) / (_mT + _mR), lf);
            };
        };
    };
};

#endif