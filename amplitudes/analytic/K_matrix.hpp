// Implementation of a PWA in the scattering-length approximation
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef K_MATRIX_HPP
#define K_MATRIX_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "partial_wave.hpp"

namespace jpacPhoto
{
    namespace analytic
    {
        class K_matrix : public raw_partial_wave
        {
            public: 

            K_matrix(amplitude_key key, kinematics xkinem, int J, std::string id = "K_matrix")
            : raw_partial_wave(key, xkinem, J, "K_matrix", id)
            {
                set_N_pars(2);
                check_QNs(xkinem);
            };

            // -----------------------------------------------------------------------
            // Virtuals 

            // These are projections onto the orbital angular momentum and therefore
            // helicity independent
            inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
            {
                // Arbitrarily pick one of the helicities to evaluate
                if (helicities != _kinematics->helicities(0)) return 0;

                // Save inputes
                store(helicities, s, t);
                
                // Normalization here to get rid of helicity dependence in amplitude::probability_distribution
                // First a 2 removes the factor 1/4 when averaging over initial helicities
                // then a 1/sqrt(2) removes the factor of 2 from the parity relation in amplitude::update_cache
                return sqrt(2) * (2*_J+1) * legendre(_J, cos(_theta)) * evaluate();
            };

            // We can have any quantum numbers
            // but for now explicitly put only pseudo-scalar, vector, and axial-vector
            // and either parity spin-1/2
            inline std::vector<std::array<int,2>> allowed_meson_JP() { return { PSUEDOSCALAR, VECTOR, AXIALVECTOR }; };
            inline std::vector<std::array<int,2>> allowed_baryon_JP(){ return { HALFPLUS, HALFMINUS }; };

            // Parameter names are a[J] and b[J] for scattering length and normalization respectively
            inline std::vector<std::string> parameter_labels()
            {
                std::vector<std::string> labels = { J_label("N"), J_label("A") };
                if (_option == EffectiveRange) labels.push_back(  J_label("B") );
                return labels;
            };
            
            // The amplitude_option allows each partial wave to switch
            // between the effective range and scattering range approximations
            inline void set_option(amplitude_option x)
            {
                switch (x)
                {
                    case ScatteringLength: 
                    {
                        set_N_pars(2); 
                        _B = 0;
                        break;
                    };
                    case EffectiveRange:   
                    { 
                        set_N_pars(3); 
                        _B = 0;
                        break;
                    };
                    default: option_error(); return;
                };

                _option = x;
            };

            // PW is the unitarized K-matrix form
            inline complex evaluate()
            {
                complex K = pow(q2(), _J) * (_A + _B*q2());
                complex M = K / (1 + G()*K);
                complex N = pow(pq(), _J) * _N;

                return N * (1 - G()*M);
            };

            protected:

            inline void allocate_parameters(std::vector<double> pars)
            {
                // Two parameters by default
                _N = pars[0];
                _A = pars[1];

                // But set the third effective range if the option is set
                if (_option == EffectiveRange) _B = pars[2];
                return;
            };

            // -----------------------------------------------------------------------
            private:

            // Free parameters N_L, A_L, and B_L
            double _N = 0, _A = 0, _B = 0; 

            // Chew-Mandelstam phase-space
            inline complex G(double m1, double m2)
            {
                complex rho, xi;
                complex result;

                rho    = csqrt(Kallen(_s, m1*m1, m2*m2)) / _s;
                xi     = 1 - (m1+m2)*(m1+m2)/_s;
                result = (rho*log((xi + rho) / (xi - rho)) - xi*(m2-m1)/(m2+m1)*log(m2/m1)) / PI;
                return result;
            };
            inline complex G(){ return G(_mX, _mR); };

            // Product of momenta for the photoproduction process
            inline complex pq()
            {
                return _kinematics->initial_momentum(_s) * _kinematics->final_momentum(_s);
            };
            
            // Product of momenta of the hadronic rescattering process
            inline complex q2()
            {
                return pow(_kinematics->final_momentum(_s), 2);
            };
        };
    };
};

#endif