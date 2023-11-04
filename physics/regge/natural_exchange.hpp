// Generic Reggeon exchange amplitude with a natural parity exchange
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#ifndef NATURAL_EXCHANGE_HPP
#define NATURAL_EXCHANGE_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "cgamma.hpp"

namespace jpacPhoto
{
    class natural_exchange : public raw_amplitude
    {
        public: 

        natural_exchange(amplitude_key key, kinematics xkinem, std::string id = "natural_exchange")
        : raw_amplitude(key, xkinem, id), _J(1), _m(M_RHO)
        {
            initialize(6);
        };

        // -----------------------------------------------------------------------
        // REQUIRED VIRTUAL FUNCTIONS

        inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
        {
            store(helicities, s, t);
            return top_coupling() * bottom_coupling() * propagator();
        };

        inline helicity_frame native_helicity_frame(){ return T_CHANNEL; };

        // Specify which final state particles amplitude can acommodate
        inline std::vector<particle> allowed_mesons() { return {AXIALVECTOR}; };
        inline std::vector<particle> allowed_baryons(){ return {THREEPLUS};   };

        inline void allocate_parameters(std::vector<double> pars)
        {
            _gT[0] = pars[0];
            _gT[1] = pars[1];
            _gB[0] = pars[2];
            _gB[1] = pars[3];
            _gB[2] = pars[4];

            _aP = pars[5];
        };

        protected:

        // Coupling constants 
        // gT = gamma - b1    - J 
        // gB = gamma - Delta - J
        std::array<double,2> _gT;
        std::array<double,3> _gB;

        // Exchange trajectory parameters
        double _aP;     // Slope
        double _m;      // Mass of the lowest resonance on trajectory
        int    _J;      // Spin of lowest exchange
        double _s0 = 1; // Fixed scale

        // Photon - B1- J^P coupling
        inline complex top_coupling()
        {
            complex q  = _kinematics->initial_momentum_tframe(_s);
            complex EX = (_t + _mX*_mX) / csqrt(4*_t);

            complex term1 = (abs(_lamX) == 1) * _gT[0] * csqrt(_t);
            complex term2 =                   - _gT[1] * q;

            if ( _lamB == 0 )
            {
                term1 *=         EX / _mX;
                term2 *= - sqrt(_t) / _mX;
            };
            
            return _lamB * (term1 + term2);
        };

        // N - Delta - J^P coupling
        inline complex bottom_coupling()
        {
            // For transparency get the signs of helicities seperately
            int sT = _lamT / abs(_lamT);
            int sR = _lamR / abs(_lamR);

            // Momentum
            complex p  = _kinematics->final_momentum_tframe(_s);
            complex ER = (_t - _mT*_mT + _mR*_mR) / csqrt(4*_t);

            complex term1 = sT * _gB[0] * (omega(+1, +1)+omega(-1, -1))*clebsch(3);
            if (_lamR == _lamT) term1 *= ER / _mR;

            // The following terms only contribute to 1/2 projection of recoil
            complex term2 = 0, term3 = 0;
            if (abs(_lamR) == 1)
            {
                term2 = (_lamR == _lamT) ? _gB[1] *p*p *csqrt(_t)/_mR * sqrt(2/3.) * (sT*omega(+1,+1) + sR*omega(-1,-1)) : 0;
                term3 =               sR * _gB[2] *p   *csqrt(_t)/_mR * sqrt(2)    * (sT*omega(+1,-1) - sR*omega(-1,+1)) * clebsch(1);
            }

            return term1 + term2 + term3;
        };

        // Typical Regge propagator
        inline complex propagator()
        {
            double alpha      = _aP *( _t - _m*_m);

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

        // Clebsch Gordan coefficients for the spin projection spinR
        // clebsch(x) = < 1, m1 ; 1/2, m2 |J, m >
        inline double clebsch(int J)
        {   
            // Helicity projecitons
            int m1 = (_lamR - _lamT) / 2;
            int m2 =  _lamT;
            int M  =  _lamR;

            int phase = 1;
            if (M < 0)
            {
                M *= -1; m1 *= -1; m2 *= -1;
                phase = -(J == 1) + !(J == 1);
            }

            int Jm1m2 = J * 100 + m1 * 10 + m2;
            double result;
            switch (Jm1m2)
            {
                case  311: result =            1; break;
                case -311: result =   sqrt(1/3.); break;
                case  301: result =   sqrt(2/3.); break;
                case -111: result =   sqrt(2/3.); break;
                case  101: result = - sqrt(1/3.); break;
                default : return 0;
            }

            return phase * result;
        }
    };
};

#endif