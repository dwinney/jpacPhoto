// Semi-inclusive production of axial vectors vector exchange using proton structure
// functions.
// Use this form if we are ignoring cross terms between individual exchanges
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Joint Physics Analysis Center (JPAC)
//               Universitat Bonn, HISKP
// Email:        daniel.winney@iu.alumni.edu
//               winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------
// REFERENCES:
//
// [1] - https://arxiv.org/abs/2404.05326
// ------------------------------------------------------------------------------

#ifndef INCLUSIVE_PHOTON_EXCHANGE_HPP       
#define INCLUSIVE_PHOTON_EXCHANGE_HPP

#include "constants.hpp"
#include "semi_inclusive.hpp"
#include "cgamma.hpp"
#include "semi_inclusive/CB_F.hpp"
#include "semi_inclusive/DL_F.hpp"

namespace jpacPhoto
{
    namespace inclusive
    {
        class photon_exchange : public raw_semi_inclusive
        {
            public: 

            photon_exchange(key k, kinematics kin, double mE, std::string id = "")
            : raw_semi_inclusive(k, kin, id), 
              _photon(is_zero(mE)), _mEx2(mE*mE),
              F1(new_inclusive_function<CB_F>(1, kProton)),
              F2(new_inclusive_function<CB_F>(2, kProton))
            { set_N_pars(3); };

            // Minimum mass is the proton 
            inline double minimum_M2(){ return pow(M_PROTON + M_PION, 2); };

            // Only free parameters are the top photocoupling and the form factor cutoff
            inline void allocate_parameters(std::vector<double> pars)
            {
                _g     = pars[0];
                _eta   = pars[1];
                _lam2  = pars[2]*pars[2];

                // Cutoff for regge amplitude
                if (!is_zero(_lam2)) _cutoff = exp(1./_lam2/_alphaP) / _alphaP;
            };

            // The invariant cross section used S T and M2 as independent vareiables
            inline double invariant_xsection(double s, double t, double M2)
            {
                if (!_photon && _regge && -t > _cutoff) return 0;

                store( s, t, M2);  // Sync kinematics
                update(s, t, M2);  // Recalculate form factors

                double PTdotW = 0;
                for (int i = 0; i < _gammas.size(); i++)
                {
                    double tpiece, spiece, alpha;
                    if (_photon || !_regge)
                    {
                        alpha  = 1;
                        tpiece = pow(_mEx2 - t, -2);
                        spiece = 2*_pdotk / M2;
                    }
                    else
                    {
                        alpha = _alpha0 + _alphaP * t;
                        tpiece = std::norm(_alphaP * cgamma(1. - alpha) * (1.-exp(-I*PI*alpha))/2);
                        spiece = s / M2;
                    };
                    
                    PTdotW += _gammas[i] * tpiece * pow(spiece, 2*alpha - i);
                };

                // Flux factor
                double flux = 1/(2*sqrt(s)*qGamma(s));

                // Form factor in the case of massive vectors
                double tprime  = t - TMINfromM2(s, M2_PROTON);
                double beta_ex = (is_zero(_lam2)) ? 1. : exp(tprime/_lam2)/pow(1-tprime/0.71, -2);

                // The form factor of the top vertex
                double mX2 = pow(_kinematics->get_meson_mass(), 2);
                if (_photon) beta_ex = pow(mX2 / (mX2 - t), 2.);

                // in nanobarn!!!!!
                return flux * pow(beta_ex*_g*_eta*E, 2) * PTdotW / (8*PI*PI) / (2.56819E-6); 
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

            // Everything regarding reggeization is handled in the code
            // except when a photon exchange is evaluated at high energies
            // even through photon doesnt reggeize we have to treat this 
            // case special to avoid probing Q2 > 300 GeV^2!
            virtual inline void reggeized(bool x) { _regge = x; set_option(_option); };

            protected:

            // Calculate dot products and form factors
            inline void update(double s, double t, double M2)
            {
                // Dot products of all the relevant momenta
                _pdotk = (s  - M2_PROTON)     / 2;
                _pdotq = (M2 - M2_PROTON - t) / 2;
                _kdotq = (t - _mX2)           / 2;

                // Production form factors
                _prefactors = 1./2*t*t/_mX2/_mX2/_mX2;
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

            // Whether we consider a photon exchange
            bool _photon = true;

            // Free parameters
            double _g    = 0; // Top couplings
            double _eta  = 1; // VMD coupling for exchange
            double _mEx2 = 0; // Mass (squared) of exchange
            double _lam2 = 1; // Form factor cutoff in GeV

            // Internal variables
            double _pdotk, _pdotq, _kdotq;
            double _prefactors, _T1, _T2;
            double _F1, _F2;
            std::array<double,3> _gammas;

            // Proton form factors
            inclusive_function F1, F2;

            // Regge trajectory parameters
            double _alpha0 = 0.5, _alphaP = 0.9;
            double _cutoff;

        };
    };
};

#endif