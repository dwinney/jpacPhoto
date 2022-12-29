// Implementation of a PWA in the scattering-length approximation with two
// coupled channels
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef COUPLED_K_MATRIX_HPP
#define COUPLED_K_MATRIX_HPP

#include "constants.hpp"
#include "kinematics.hpp"

namespace jpacPhoto
{
    namespace analytic
    {
        class coupled_K_matrix : public raw_partial_wave
        {
            public: 

            coupled_K_matrix(amplitude_key key, kinematics xkinem, int J, std::array<double,2> masses, std::string id = "coupled_K_matrix")
            : raw_partial_wave(key, xkinem, J, "coupled_K_matrix", id),
              _m1(masses[0]), _m2(masses[1])
            {
                // 3 K-matrix parameters and 2 normalizations
                set_N_pars(5);
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
                switch (_option)
                {
                    case Default:         return {J_label("a00"), J_label("a01"), J_label("a11"), J_label("b0"), J_label("b1") };
                    case EffectiveRange:  return {J_label("a00"), J_label("a01"), J_label("a11"), J_label("r00"), J_label("r11"), J_label("b0"), J_label("b1") };
                    default: return {{}};
                };
            };

            // The amplitude_option to choose a limiting case of parameters
            // or Default to reset to the most general parameterization
            inline void set_option(amplitude_option x)
            {
                switch (x)
                {
                    case Default:
                    {
                        set_N_pars(5);
                        _r00 = 0; _r11 = 0;
                        break;
                    };
                    case EffectiveRange:   
                    { 
                        set_N_pars(7); 
                        break;
                    };
                    default: option_error(); return;
                };

                _option = x;
            };

            // PW is the unitarized K-matrix form
            inline complex evaluate()
            {
                // Recalculate momenta
                _q[0] = _kinematics->final_momentum(_s);
                _q[1] = csqrt(Kallen(_s, _m1*_m1, _m2*_m2)) / csqrt(4.*_s);

                // Calculate loop functions
                _G[0] = G(_mX, _mR); 
                _G[1] = G(_m1, _m2);

                // Production amplitude pieces
                _B[0] = pow(pq(0), _J) * _b0;
                _B[1] = pow(pq(1), _J) * _b1;

                // K-matrices
                _K00 = pow(q2(0,0), _J) * (_a00 + _r00*q2(0,0)/2);
                _K01 = pow(q2(0,1), _J) * _a01;
                _K11 = pow(q2(1,1), _J) * (_a11 + _r11*q2(1,1)/2);

                // The A-matrices all share the same denominator
                _D    = (1-_G[0]*_K00)*(1-_G[1]*_K11) - _G[0]*_G[1]*_K01*_K01;

                // Determinant of K-matrix
                _delK = _K00*_K11 - _K01*_K01;

                // Then all we need are the numerators
                _A00 = (_K00 - _G[1]*_delK) / _D;
                _A01 = _K01 / _D;

                return _B[0] * (1 + _G[0]*_A00) + _B[1]*_G[1]*_A01;
            };

            protected:

            inline void allocate_parameters(std::vector<double> pars)
            {
                switch (_option)
                {
                    case Default: 
                    {
                        _a00 = pars[0];
                        _a01 = pars[1];
                        _a11 = pars[2];
                        _b0  = pars[3];
                        _b1  = pars[4];
                        break;
                    }
                    case EffectiveRange: 
                    {
                        _a00 = pars[0];
                        _a01 = pars[1];
                        _a11 = pars[2];
                        _r00 = pars[3];
                        _r11 = pars[4];
                        _b0  = pars[5];
                        _b1  = pars[6];
                        break;
                    }
                    default: return;
                }
            };

            // -----------------------------------------------------------------------
            private:

            // Mass of the second open channel
            double _m1 = 0,  _m2 = 0;
            
            // K-matrix parameters
            double _a00 = 0, _a01 = 0, _a11 = 0;
            double _r00 = 0, _r11 = 0;
            
            // Production amplitude parameters
            double _b0 = 0,  _b1 = 0; 

            // Internal variables for the K and A amplitudes
            complex _K00, _K01, _K11, _A00, _A01, _A11, _D, _delK;
            std::array<complex,2> _q, _G, _B;

            // Chew-Mandelstam phase-space
            inline complex G(double m1, double m2)
            {
                complex rho, xi;
                complex result;

                rho    = csqrt(Kallen(_s, m1*m1, m2*m2)) / _s;
                xi     = 1 - (m1+m2)*(m1+m2)/_s;
                result = - (rho*log((xi + rho) / (xi - rho)) - xi*(m2-m1)/(m2+m1)*log(m2/m1)) / PI;
                return result;
            };

            // Product of momenta for the photoproduction process
            inline complex pq(unsigned i){ return _kinematics->initial_momentum(_s) * _q[i]; };
            
            // Product of momenta of the hadronic rescattering process
            inline complex q2(unsigned i, unsigned j){ return _q[i] * _q[j]; };
        };
    };
};

#endif