// Generic Reggeon exchange amplitude with a unnatural parity exchange
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#ifndef UNNATURAL_EXCHANGE_HPP
#define UNNATURAL_EXCHANGE_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "cgamma.hpp"

namespace jpacPhoto
{
    class unnatural_exchange : public raw_amplitude
    {
        public: 

        unnatural_exchange(amplitude_key key, kinematics xkinem, std::string id = "unnatural_exchange")
        : raw_amplitude(key, xkinem, id), _J(0), _m(M_PION)
        {
            initialize(3);
        };

        // -----------------------------------------------------------------------
        // REQUIRED VIRTUAL FUNCTIONS

        inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
        {
            store(helicities, s, t);

            // An on-shell spin-0 exchange cannot flip helicities
            bool helicity_conserving = (_lamB == _lamX) && (_lamR == _lamT);
            if (!helicity_conserving) return 0;

            return top_coupling() * bottom_coupling() * propagator();
        };

        inline helicity_frame native_helicity_frame(){ return T_CHANNEL; };

        // Specify which final state particles amplitude can acommodate
        inline std::vector<particle> allowed_mesons() { return {AXIALVECTOR}; };
        inline std::vector<particle> allowed_baryons(){ return {THREEPLUS};   };

        inline void allocate_parameters(std::vector<double> pars)
        {
            _gT = pars[0];
            _gB = pars[1];
            _aP = pars[2];
        };

        protected:

        // Coupling constants 
        // gT = gamma - b1    - J 
        // gB = gamma - Delta - J
        double _gT = 1, _gB = 1;

        // Exchange trajectory parameters
        double _aP; // Slope
        double _m;  // Mass of the lowest resonance on trajectory
        int    _J;  // Spin of lowest exchange
        double _s0 = 1; // Fixed scale

        // Photon - B1- J^P coupling
        inline complex top_coupling()
        {
            complex q = _kinematics->initial_momentum_tframe(_s);
            return _gT * csqrt(_t) * q;
        };

        // N - Delta - J^P coupling
        inline complex bottom_coupling()
        {
            // Momentum
            complex p  = _kinematics->final_momentum_tframe(_s);
            return _gB*sqrt(2./3)*csqrt(_t)/_mR*p*(omega(+1,-1)+omega(-1, +1));
        };

        // Typical Regge propagator
        inline complex propagator()
        {
            double alpha = _aP *( _t - _m*_m);

            // Signature calculated from the lowest resonance
            int tau = (_J % 2 == 0) - !(_J % 2 == 0);
            complex signature = (1 + tau*exp(-I*PI*alpha))/2;

            return _aP * signature * cgamma(_J - alpha) * pow(_s/_s0, alpha);
        };

        // Shortcut for product of energy factors that appear from spinors
        // x and y are signs for the target and recoil sqrt(E + \pm m) respetively
        inline complex omega(int x, int y)
        {
            complex ET = (_t + _mT*_mT - _mR*_mR) / csqrt(4*_t);
            complex ER = (_t - _mT*_mT + _mR*_mR) / csqrt(4*_t);
            return sqrt( (ET + x*_mT)*(ER + y*_mR) );
        };
    };
};

#endif