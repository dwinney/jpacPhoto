// Semi-inclusive production of axial vectors via meson ex using proton structure
// functions.
// This form combines the rho and omega exchanges to incorporate the non-diagonal contirbution 
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
//               Universitat Bonn, HISKP
// Email:        daniel.winney@iu.alumni.edu
//               winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------
// REFERENCES:
//
// [1] - https://arxiv.org/abs/2404.05326
// ------------------------------------------------------------------------------

#ifndef INCLUSIVE_VECTOR_EXCHANGE_HPP       
#define INCLUSIVE_VECTOR_EXCHANGE_HPP

#include "constants.hpp"
#include "semi_inclusive.hpp"
#include "inclusive/structure_functions.hpp"

namespace jpacPhoto
{
    namespace inclusive
    {
        class vector_exchange : public raw_semi_inclusive
        {
            public: 

            vector_exchange(key k, kinematics kin)
            : raw_semi_inclusive(k, kin, "vector_exchange")
            { set_N_pars(3); };

            // Minimum mass is the proton 
            inline double minimum_M2(){ return pow(M_PROTON + M_PION, 2); };

            // Only free parameters are the top photocoupling and the form factor cutoff
            inline void allocate_parameters(std::vector<double> pars)
            {
                _g[0] = pars[0]; _g[1] = pars[1]; _g[2] = pars[2];
            };

            // The invariant cross section used S T and M2 as independent vareiables
            inline double invariant_xsection(double s, double t, double M2)
            {
                store( s, t, M2);  // Sync kinematics
                update(s, t, M2);  // Recalculate form factors

                // Form factors for rho and omega
                double tprime  = t - TMINfromM2(s, M2_PROTON);
                std::array<double,3> beta = {  std::norm(_mX2 / (_mX2 - t)),
                                               exp(tprime/_lam2[1])/pow(1-tprime/0.71,-2), 
                                               exp(tprime/_lam2[2])/pow(1-tprime/0.71,-2)} ;

                // Flux factor
                double flux = 1/(2*sqrt(s)*qGamma(s));

                double TdotW = 0;
                for (int i = 0; i < _gammas.size(); i++) TdotW += _gammas[i] * pow(2*_pdotk/M2, 2 - i);
                
                if (!_regge)
                {
                    double propagators = 0;
                    for (int i = 0; i < 3; i++) propagators += beta[i]*_g[i]*_eta[i]/(_mEx2[i] - t);
                    
                    // in nanobarn!!!!!
                    return flux * E*E/4 * std::norm(propagators) * TdotW / (8*PI*PI) * HBARC; 
                }
                else
                {
                    double photon = std::norm(_g[0]*beta[0]/-t);

                    // Reggeized version of the TdotW tensor contraction
                    double reggeon = 0, TdotW_R = 0, alpha = _alpha0 + _alphaP * t;
                    if (-t < _cutoff)
                    {
                        for (int i = 0; i < _gammas.size(); i++)  TdotW_R += _gammas[i] * pow(s/M2, 2*alpha - i);
                        reggeon  = std::norm(beta[1]*_g[1]*_eta[1]/2. + beta[2]*_g[2]*_eta[2]/2.);
                        reggeon *= std::norm(_alphaP * cgamma(1. - alpha) * (1.-exp(-I*PI*alpha))/2.);
                    };

                    // in nanobarn!!!!!
                    return flux * E*E * (photon*TdotW + reggeon*TdotW_R) / (8*PI*PI) * HBARC; 
                };
            };

            // Options select proton or neutron target
            static const int kNotReggeized = 0, kReggeized = 1;
            inline void set_option (int opt)
            {
                switch (opt)
                {
                    case kNotReggeized:
                    {
                        _regge = false; struct_funcs.set_option(structure_functions::kCB);
                        break;
                    };
                    case kReggeized:
                    {
                        _regge = true; struct_funcs.set_option(structure_functions::kDL);
                        break;
                    };
                    default: struct_funcs.set_option(opt);
                };
            };

            protected:

            // Calculate dot products and form factors
            inline void update(double s, double t, double M2)
            {
                // Dot products of all the relevant momenta
                _pdotk = (s  - M2_PROTON)     / 2;
                _pdotq = (M2 - M2_PROTON - t) / 2;
                _kdotq = (t - _mX2)           / 2;

                // The form factors of the top vertex
                double beta_Qgg = pow(1. - t/_mX2, -2.);
                _prefactors = beta_Qgg/2*t*t/_mX2/_mX2/_mX2;

                _T1 = _prefactors * _kdotq*_kdotq;
                _T2 = _prefactors * _kdotq*(_mX2 - 2*_kdotq);

                // Hadronic structure function
                _F1 = struct_funcs.F1(M2, t);
                _F2 = struct_funcs.F2(M2, t);

                // Expansion in cross section in terms of flip amplitudes
                _gammas[0] =   M2*M2/(4*_pdotq*_kdotq)* _T2*_F2;
                _gammas[1] = - M2/t * _T2*_F2;
                _gammas[2] =   3*_F1*_T1 + (_kdotq/t)*_F1*_T2 + (_pdotq/t - M2_PROTON/_pdotq)*_F2*_T1 + _kdotq*_pdotq/t/t*_F2*_T2;
            };  

            inline bool use_TX(){ return false; };

            private:

            // Free parameters
            std::array<double,3> _g;
            std::array<double,3> _mEx2 = {0., M_RHO*M_RHO, M_OMEGA*M_OMEGA};
            std::array<double,3> _eta  = {1, 16.37,   56.34};
            std::array<double,3> _lam2 = {0., 1.4*1.4, 1.2*1.2  };

            // Internal variables
            double _pdotk, _pdotq, _kdotq;
            double _prefactors, _T1, _T2;
            double _F1, _F2;
            std::array<double,3> _gammas;

            // Proton form factors
            structure_functions struct_funcs;

            // Regge trajectory parameters
            double _alpha0 = 0.5, _alphaP = 0.9;
            double _cutoff = exp(1./_lam2[1]/_alphaP) / _alphaP;
        };
    };
};

#endif