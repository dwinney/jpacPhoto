// J/psi photoproduction under the assumption of QCD factorization in a GPD approach
// Adapted from the models in [1] and [2].
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Helmholtz-Instituts für Strahlen- und Kernphysik (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// References:
// [1] - https://arxiv.org/abs/2103.11506
// [2] - https://arxiv.org/abs/2305.06992
// ------------------------------------------------------------------------------

#ifndef JPSI_GPD_HPP
#define JPSI_GPD_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

namespace jpacPhoto { namespace analytic {

    class jpsi_GPD : public raw_amplitude
    {
        public: 

        jpsi_GPD(key k, kinematics xkinem)
        : raw_amplitude(k, xkinem, "GPD")
        { initialize(0); };

        inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
        {
            store(helicities, s, t);

            // Helicity conservation at the top vertex (VMD)
            complex top_vertex;
            if (_schc) { top_vertex = (_lamB == _lamX); }
            else
            {
                _covariants->update(helicities, s, t);
                auto eps_gam = _covariants->eps();
                auto eps_psi = _covariants->eps_prime();
                top_vertex   = contract(eps_gam, eps_psi);
            };

            return top_vertex*E*_eQ*(16*PI*_as)/sqrt(3*M_JPSI*M_JPSI*M_JPSI)*_psi0*G();
        }

                // Covariants are always natively s-channel helicity amplitudes
        inline helicity_frame         native_helicity_frame(){ return S_CHANNEL; };
        inline std::vector<quantum_numbers> allowed_mesons() { return { VECTOR }; };
        inline std::vector<quantum_numbers> allowed_baryons(){ return { HALFPLUS    }; };

        static const int k2021    = 0;
        static const int k2023    = 1;
        static const int kSCHC    = 2;
        static const int kNotSCHC = 3;

        inline void set_option(int opt)
        {
            switch (opt)
            {
                case (k2021)    : { n = 2; _a0 = 0.580; _c0 = -1.00; _ma = 1.13; _mc = 0.48; return; };
                case (k2023)    : { n = 3; _a0 = 0.414; _c0 = -1.21; _ma = 2.07; _mc = 0.91; return; };
                case (kSCHC)    : { _schc = true;  return; };
                case (kNotSCHC) : { _schc = false; return; };
                default: option_error();
            };
            return;
        };

        private:

        // Constants 
        double _as = 0.3, _eQ = 1/3.;
        double _psi0 = sqrt(1.0952/4/PI);
        int n = 2;
        bool _schc = true;

        // Form factors
        inline double A(){ return _a0 / pow(1. - _t/_ma/_ma, n); };
        inline double C(){ return _c0 / pow(1. - _t/_mc/_mc, n); };

        inline complex G()
        {
            double xi2  = pow((_mX*_mX - _t)/(_s - _u), 2);
            double tau2 = _t/(4*M2_PROTON); 

            complex helstruct = (_lamT == _lamR) ? _lamT*csqrt((tau2 - 1.)*xi2 - tau2)*(4*xi2)*C()
                                                : (1-xi2)*A() + 4*xi2*C();
            return (4/xi2)*helstruct/sqrt(1-xi2);
        };

        double _a0 =  0.58, _ma  = 1.13;
        double _c0 = -1.00, _mc  = 0.48;
    };
}; };
#endif