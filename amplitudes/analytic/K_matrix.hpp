// Implementation of a PWA in the scattering-length approximation with up to three
// coupled channels
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
#include "partial_wave.hpp"

namespace jpacPhoto
{
    namespace analytic
    {
        class K_matrix : public raw_partial_wave
        {
            public: 

            // Single channel K-matrix
            K_matrix(amplitude_key key, kinematics xkinem, int J, std::string id = "K_matrix")
            : raw_partial_wave(key, xkinem, J, "K_matrix", id)
            {
                // 3 K-matrix parameters and 2 normalizations
                set_N_pars(2);
                check_QNs(xkinem);

                _Nth = 1;

                // Populate the thresholds
                _thresholds[0] = {xkinem->get_meson_mass(), xkinem->get_recoil_mass()};
            };

            // Two-channel K-matrix
            K_matrix(amplitude_key key, kinematics xkinem, int J, std::array<double,2> masses, std::string id = "K_matrix")
            : raw_partial_wave(key, xkinem, J, "K_matrix", id)
            {
                // 3 K-matrix parameters and 2 normalizations
                set_N_pars(5);
                check_QNs(xkinem);

                _Nth = 2;

                // Populate the thresholds
                _thresholds[0] = {xkinem->get_meson_mass(), xkinem->get_recoil_mass()};
                _thresholds[1] = masses;
            };

            // Three-channel K-matrix
            K_matrix(amplitude_key key, kinematics xkinem, int J, std::array<std::array<double,2>,2> masses, std::string id = "K_matrix")
            : raw_partial_wave(key, xkinem, J, "K_matrix", id)
            {
                // 3 K-matrix parameters and 2 normalizations
                set_N_pars(9);
                check_QNs(xkinem);

                _Nth = 3;

                // Populate the thresholds
                _thresholds[0] = {xkinem->get_meson_mass(), xkinem->get_recoil_mass()};
                _thresholds[1] = masses[0];
                _thresholds[2] = masses[1];
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
                std::vector<std::string> labels;
                switch (_Nth)
                {
                    case 1: { labels = { J_label("N"), J_label("A") }; break; }
                    case 2: {
                              labels = { J_label("N0"), J_label("N1"), J_label("A00"), J_label("A01"), J_label("A11") }; break;
                    }
                    case 3: {
                              labels = { J_label("N0"), J_label("N1"), J_label("N2"), J_label("A00"), J_label("A01"), J_label("A02"),
                                                                                                      J_label("A11"), J_label("A12"),
                                                                                                                      J_label("A22") }; break;
                    }
                };

                // If EffectiveRange, add extra labels
                switch (_option)
                {
                    case Default: return labels;
                    case EffectiveRange:  
                    {
                        if (_Nth == 1)
                        {
                            labels.push_back( J_label("B"));
                            return labels;
                        };

                        labels.push_back( J_label("B00") );
                        labels.push_back( J_label("B11") );
                        if (_Nth == 2) return labels;

                        labels.push_back( J_label("B22") );
                        return labels;
                    }
                    default:  return {{}};
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
                        set_N_pars(2+(_Nth>1)*(4*_Nth-5));
                        _B00 = 0; _B11 = 0; _B22 = 0;
                        break;
                    };
                    case EffectiveRange:   
                    { 
                        set_N_pars(3+(_Nth>1)*(5*_Nth-6)); 
                        break;
                    };
                    default: option_error(); return;
                };

                _option = x;
            };

            // PW is the unitarized K-matrix form
            inline complex evaluate()
            {
                // Recalculate production related quantities for each threshold
                for (int i = 0; i < _Nth; i++)
                {
                    _q[i] = q(i);  // Break up momenta
                    _G[i] = G(i);  // Phase-space factor

                    // Production amplitude
                    _N[i] = pow(pq(i), _J) * _Nhat[i]; 
                }
                // K-matrices elements
                _K00 = pow(q2(0,0), _J) * (_A00 + _B00*q2(0,0));

                _K01 = pow(q2(0,1), _J) *  _A01;
                _K11 = pow(q2(1,1), _J) * (_A11 + _B11*q2(1,1));

                _K02 = pow(q2(0,2), _J) *  _A02;
                _K12 = pow(q2(1,2), _J) *  _A12;
                _K22 = pow(q2(2,2), _J) * (_A22 + _B22*q2(2,2));

                // Determinant of K-matrix
                _detK = _K00*_K11*_K22 + 2*_K01*_K02*_K12 - _K02*_K02*_K11 - _K12*_K12*_K00 - _K01*_K01*_K22;

                // Determinants of the sub-matrices
                _detK01 = _K00*_K11 - _K01*_K01;
                _detK02 = _K00*_K22 - _K02*_K02;
                _detK12 = _K11*_K22 - _K12*_K12;

                // The M-matrices all share the same denominator    
                _D    = (1+_G[0]*_K00)*(1+_G[1]*_K11)*(1+_G[2]*_K22) - _G[0]*_G[1]*_K01*_K01
                                                                     - _G[0]*_G[2]*_K02*_K02
                                                                     - _G[1]*_G[2]*_K12*_K12 - _G[0]*_G[1]*_G[2]*(_K00*_K11*_K22 - _detK);

                // Then all we need are the numerators
                _M00 = (_K00 + _G[1]*_detK01 + _G[2]*_detK02 + _G[1]*_G[2]*(_detK12*_K00 + _detK02*_K11 + _detK01*_K22 - 2*_K00*_K11*_K22 + 2*_K01*_K02*_K12)) / _D;
                _M01 = (_K01 + _G[2]*(_K01*_K22 - _K02*_K12)) / _D;
                _M02 = (_K02 + _G[1]*(_K02*_K11 - _K01*_K12)) / _D;

                complex T = _N[0] * (1 - _G[0]*_M00) - _N[1]*_G[1]*_M01 - _N[2]*_G[2]*_M02;

                return T;
            };

            protected:

            inline void allocate_parameters(std::vector<double> pars)
            {
                switch (_Nth)
                {
                    default: return;
                    case 1: 
                    {
                        _Nhat[0] = pars[0];
                        _A00     = pars[1];
                        break;
                    };
                    case 2:
                    {
                        _Nhat[0] = pars[0];
                        _Nhat[1] = pars[1];
                        _A00     = pars[2];
                        _A01     = pars[3];
                        _A11     = pars[4];
                        break;
                    };
                    case 3:
                    {
                        _Nhat[0] = pars[0];
                        _Nhat[1] = pars[1];
                        _Nhat[2] = pars[2];
                        _A00     = pars[3];
                        _A01     = pars[4];
                        _A02     = pars[5];
                        _A11     = pars[6];
                        _A12     = pars[7];
                        _A22     = pars[8];
                        break;
                    };
                };

                if ( _option == EffectiveRange )
                {
                    _B00 = pars[2+(_Nth>1)*(4*_Nth-5)];
                    _B11 = (_Nth > 1)*pars[4*_Nth-2];
                    _B22 = (_Nth > 1)*pars[4*_Nth-1];
                }
            };

            // -----------------------------------------------------------------------
            private:

            // If we use three channels
            int _Nth = 0;

            // Mass of intermediate coupled channels
            // we can have up to two additional channels
            std::array<std::array<double,2>,3> _thresholds;
            
            // Production parameters
            std::array<complex,3> _q, _G, _N;

            // K-matrix parameters
            double _A00 = 0, _A01 = 0, _A02 = 0;
            double           _A11 = 0, _A12 = 0;
            double                     _A22 = 0;

            // Second order parameters (effective range)
            double _B00 = 0, _B11 = 0, _B22 = 0;
            
            // Production amplitude parameters
            std::array<double,3> _Nhat;

            // Internal variables for the K matrices
            complex _K00 = 0, _K01 = 0, _K02 = 0;
            complex           _K11 = 0, _K12 = 0;
            complex                     _K22 = 0;
            complex  _D = 0,      _detK = 0;
            complex  _detK01 = 0, _detK02 = 0, _detK12 = 0;

            // As well as the final M matrixes
            complex _M00 = 0, _M01 = 0, _M02 = 0;

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
            inline complex G(unsigned i){ return G(_thresholds[i][0], _thresholds[i][1]); };

            // Outgoing break-up momentum
            inline complex q(double m1, double m2)
            {
                return csqrt(Kallen(_s, m1*m1, m2*m2)) / csqrt(4.*_s);
            };
            inline complex q(unsigned i){ return q(_thresholds[i][0], _thresholds[i][1]); };

            // Product of momenta for the photoproduction process
            inline complex pq(unsigned i){ return _kinematics->initial_momentum(_s) * _q[i]; };
            
            // Product of momenta of the hadronic rescattering process
            inline complex q2(unsigned i, unsigned j){ return _q[i] * _q[j]; };
        };
    };
};

#endif