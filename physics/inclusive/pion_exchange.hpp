// Semi-inclusive production of axial vectors via pion exchange.
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#ifndef PION_EXCHANGE_HPP       
#define PION_EXCHANGE_HPP

#include "constants.hpp"
#include "semi_inclusive.hpp"
#include "inclusive/piN_xsection.hpp"

namespace jpacPhoto
{
    namespace inclusive
    {
        // Actual pion exchange inclusive amplitude
        class pion_exchange : public raw_semi_inclusive
        {
            public: 

            pion_exchange(key k, kinematics kinem, int pm)
            : raw_semi_inclusive(k, kinem, "pion_exchange")
            {
                // For pi+ production we have a pi- exchanged at the bottom vertex
                if (pm == +1) set_option(piN_xsection::kPI_MINUS);
                // likewise pi- final state involves the pi+N xsection
                else          set_option(piN_xsection::kPI_PLUS);
                set_N_pars(1);
            };

            // Minimum mass is the proton 
            inline double minimum_M2(){ return pow(M_PROTON + M_PION, 2); };

            // Only free parameters are the top photocoupling and the form factor cutoff
            inline void allocate_parameters(std::vector<double> pars)
            {
                _g      = pars[0];
            };

            // The invariant cross section used S T and M2 as independent vareiables
            inline double invariant_xsection(double s, double t, double xm2)
            {
                store(s, t, xm2);
                
                double M2, K, P_pi;
                if (_regge)
                {
                    if (is_zero(_x - 1)) return 0;

                    M2    = M2fromTX(s, t, xm2);
                    P_pi  = regge_propagator();
                    K     = (1 - _x);
                }
                else 
                {
                    M2    = xm2;
                    P_pi  = 1 / (M2_PION - t);
                    K     = sqrt(kallen(M2, t, M2_PROTON)/kallen(s, 0., M2_PROTON));

                };

                if (are_equal(M2, minimum_M2())) return 0.;

                // Total cross-section always gets the physical M2 
                double  sigmatot  = _sigma(M2, t) * 1E6; // in nb
                
                return K/(16*PI*PI*PI) * pow(coupling()*P_pi, 2) * sigmatot;
            };

            static const int kNotReggeized = 0;
            static const int kReggeized    = 1;
            inline void set_option (int opt)
            { 
                if (opt == kReggeized || opt == kNotReggeized) _regge = opt;
                else _sigma.set_option(opt); 
            };

            protected:

            inline double coupling()
            {
                // Exponential form factor
                double beta_pi =  exp((_t - TMINfromM2(_s, M2_PROTON))/_lamPi/_lamPi);

                // Scalar coupling
                double T_pi    = (_g/sqrt(_mX2)) * (_mX2 - _t)/2;

                return beta_pi * T_pi ;
            };

            inline double regge_propagator()
            {
                // Cutoff 
                if (std::abs(_t) > 8.3) return 0;

                // Trajectory
                double alpha = _alpha0 + _alphap*_t;

                // Half angle factor
                complex xi = (1. + exp(-I*PI*alpha))/2.;

                complex result = _alphap * xi * cgamma(-alpha) * pow(1 - _x, -alpha);
                return std::abs(result);
            };

            // If we are reggeized we use the high-energy approximation and use t & x
            inline bool use_TX(){ return _regge; };

            private:
            
            piN_xsection _sigma;
            int    _pm     = +1;    // Charge of the produced meson
            double _g      = 0;     // Top coupling
            double _lamPi  = 0.9;   // Exponential cut-off

            // Pion regge trajectory parameters
            double _alphap = 0.7;
            double _alpha0 = -M2_PION*0.7;
        };
    };
};

#endif