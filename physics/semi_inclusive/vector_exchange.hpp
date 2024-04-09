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
#include "cgamma.hpp"
#include "semi_inclusive/CB_F.hpp"
#include "semi_inclusive/DL_F.hpp"

namespace jpacPhoto
{
    namespace inclusive
    {
        class vector_exchange : public raw_semi_inclusive
        {
            public: 

            vector_exchange(key k, kinematics kin, std::string id = "")
            : raw_semi_inclusive(k, kin, id), 
              F1(new_inclusive_function<CB_F>(1, kProton)),
              F2(new_inclusive_function<CB_F>(2, kProton))
            { set_N_pars(2); };

            // Minimum mass is the proton 
            inline double minimum_M2(){ return pow(M_PROTON + M_PION, 2); };

            // Only free parameters are the top photocoupling and the form factor cutoff
            inline void allocate_parameters(std::vector<double> pars)
            {
                _g[0]     = pars[0]; _g[1]     = pars[1];
            };

            // The invariant cross section used S T and M2 as independent vareiables
            inline double invariant_xsection(double s, double t, double M2)
            {
                if (_regge && -t > _cutoff) return 0;

                store( s, t, M2);  // Sync kinematics
                update(s, t, M2);  // Recalculate form factors

                double TdotW = 0;

                // Form factors for rho and omega
                double tprime  = t - TMINfromM2(s, M2_PROTON);
                std::array<double,2> beta = {  exp(tprime/_lam2[0])/pow(1-tprime/0.71,-2), 
                                               exp(tprime/_lam2[1])/pow(1-tprime/0.71,-2)} ;

                double exch_piece, spiece, alpha;
                if (!_regge)
                {
                    alpha  = 1;
                    exch_piece = std::norm(beta[0]*_g[0]*_eta[0]/(_mEx2[0] - t) + beta[1]*_g[1]*_eta[1]/(_mEx2[1] - t));
                    spiece = 2*_pdotk / M2;
                }
                else
                {
                    alpha = _alpha0 + _alphaP * t;
                    exch_piece = std::norm(beta[0]*_g[0]*_eta[0] + beta[1]*_g[1]*_eta[1]) * std::norm(_alphaP * cgamma(1. - alpha) * (1.-exp(-I*PI*alpha))/2);
                    spiece = s / M2;
                };

                for (int i = 0; i < _gammas.size(); i++)  TdotW += _gammas[i] * pow(spiece, 2*alpha - i);

                // Flux factor
                double flux = 1/(2*sqrt(s)*qGamma(s));

                // in nanobarn!!!!!
                return flux * E*E/4 * exch_piece * TdotW / (8*PI*PI) * HBARC; 
            };

            // Options select proton or neutron target
            static const int kProton = 0, kNeutron = 1;
            inline void set_option (int opt)
            {
                if (_regge)
                {
                    F1 = new_inclusive_function<DL_F>(1); F2 = new_inclusive_function<DL_F>(2);
                    return;
                };

                F1 = new_inclusive_function<CB_F>(1, opt); F2 = new_inclusive_function<CB_F>(2, opt);
                return;
            };

            virtual inline void reggeized(bool x) { _regge = x; set_option(_option); };

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
                _F1 = F1->evaluate(M2, t);
                _F2 = F2->evaluate(M2, t);

                // Expansion in cross section in terms of flip amplitudes
                _gammas[0] =   M2*M2/(4*_pdotq*_kdotq)* _T2*_F2;
                _gammas[1] = - M2/t * _T2*_F2;
                _gammas[2] =   3*_F1*_T1 + (_kdotq/t)*_F1*_T2 + (_pdotq/t - M2_PROTON/_pdotq)*_F2*_T1 + _kdotq*_pdotq/t/t*_F2*_T2;
            };  

            // If we are reggeized we use the high-energy approximation and use t & x
            inline bool use_TX(){ return false; };

            private:

            // Free parameters
            std::array<double,2> _g;
            std::array<double,2> _mEx2 = {M_RHO*M_RHO, M_OMEGA*M_OMEGA};
            std::array<double,2> _eta  = {16.37,   56.34};
            std::array<double,2> _lam2 = {1.4*1.4, 1.2*1.2  };

            // Internal variables
            double _pdotk, _pdotq, _kdotq;
            double _prefactors, _T1, _T2;
            double _F1, _F2;
            std::array<double,3> _gammas;

            // Proton form factors
            inclusive_function F1, F2;

            // Regge trajectory parameters
            double _alpha0 = 0.5, _alphaP = 0.9;
            double _cutoff = exp(1./_lam2[0]/_alphaP) / _alphaP;
        };
    };
};

#endif