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
#include "inclusive_function.hpp"
#include "semi_inclusive/JPAC_piN.hpp"
#include "semi_inclusive/PDG_piN.hpp"
#include "cgamma.hpp"

namespace jpacPhoto
{
    namespace inclusive
    {

        class pion_exchange : public raw_semi_inclusive
        {
            public: 

            pion_exchange(key k, kinematics kinem, int pm, std::string id = "")
            : raw_semi_inclusive(k, kinem, id), _pm(pm),
                _sigma(new_inclusive_function<JPAC_piN>(-pm))
            {
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
                    K     = sqrt(Kallen(M2, t, M2_PROTON)/Kallen(s, 0., M2_PROTON));

                };

                if (are_equal(M2, minimum_M2())) return 0.;

                // Total cross-section always gets the physical M2 
                double  sigmatot  = _sigma->evaluate(M2, t) * 1E6; // in nb
                
                return K/(16*PI*PI*PI) * pow(coupling()*P_pi, 2) * sigmatot;
            };

            // Options are the parameterization of the sigma_tot
            static const int kJPAC  = 0;
            static const int kPDG   = 1;
            static const int kPwave = 2;
            inline void set_option (int opt)
            {
                switch (opt)
                {
                    case kJPAC:  { _sigma = new_inclusive_function<JPAC_piN>(-_pm); 
                                    _option = opt; break; };
                    case kPwave: { _sigma = new_inclusive_function<JPAC_piN>(-_pm, 1); 
                                    _option = opt; break; };
                    case kPDG:   { PDG_total_xsections pdgpm = (_pm == -1) ? pipp : pimp; 
                                    _sigma = new_PDG_sigmatot(pdgpm); 
                                    _option = opt; break; };
                    default: return;
                };
            }

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
            
            inclusive_function _sigma;
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