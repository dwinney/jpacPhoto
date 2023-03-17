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
                initialize(2);

                _Nth = 1;

                // Populate the thresholds
                _thresholds[0] = {xkinem->get_meson_mass(), xkinem->get_recoil_mass()};
            };

            // Two-channel K-matrix
            K_matrix(amplitude_key key, kinematics xkinem, int J, std::array<double,2> masses, std::string id = "K_matrix")
            : raw_partial_wave(key, xkinem, J, "K_matrix", id)
            {
                initialize(5);

                _Nth = 2;

                // Populate the thresholds
                _thresholds[0] = {xkinem->get_meson_mass(), xkinem->get_recoil_mass()};
                _thresholds[1] = masses;
            };

            // Three-channel K-matrix
            K_matrix(amplitude_key key, kinematics xkinem, int J, std::array<std::array<double,2>,2> masses, std::string id = "K_matrix")
            : raw_partial_wave(key, xkinem, J, "K_matrix", id)
            {
                initialize(9);

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
                    case 1: { labels = { J_label("n"), J_label("a") }; break; }
                    case 2: {
                              labels = { J_label("n0"), J_label("n1"), J_label("a00"), J_label("a01"), J_label("a11") }; break;
                    }
                    case 3: {
                              labels = { J_label("n0"), J_label("n1"), J_label("n2"), J_label("a00"), J_label("a01"), J_label("a02"),
                                                                                                      J_label("a11"), J_label("a12"),
                                                                                                                      J_label("a22") }; break;
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
                            labels.push_back( J_label("b"));
                            return labels;
                        };

                        labels.push_back( J_label("b00") );
                        labels.push_back( J_label("b11") );
                        if (_Nth == 2) return labels;

                        labels.push_back( J_label("b22") );
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
                        _b00 = 0; _b11 = 0; _b22 = 0;
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
                }

                // K-matrices elements
                _K00 = pow(q2(0,0), _J) * (_a00 + _b00*q2(0,0));

                _K01 = pow(q2(0,1), _J) *  _a01;
                _K11 = pow(q2(1,1), _J) * (_a11 + _b11*q2(1,1));

                _K02 = pow(q2(0,2), _J) *  _a02;
                _K12 = pow(q2(1,2), _J) *  _a12;
                _K22 = pow(q2(2,2), _J) * (_a22 + _b22*q2(2,2));

                // Determinant of K-matrix
                _detK = _K00*_K11*_K22 + 2*_K01*_K02*_K12 - _K02*_K02*_K11 - _K12*_K12*_K00 - _K01*_K01*_K22;

                // The M-matrices all share the same denominator    
                _D    = (1+_G[0]*_K00)*(1+_G[1]*_K11)*(1+_G[2]*_K22) - _G[0]*_G[1]*_K01*_K01
                                                                     - _G[0]*_G[2]*_K02*_K02
                                                                     - _G[1]*_G[2]*_K12*_K12 - _G[0]*_G[1]*_G[2]*(_K00*_K11*_K22 - _detK);

                // Then all we need are the numerators
                _N[0] =   pow(pq(0), _J) * _n[0] * ((1 + _G[1]*_K11)*(1 + _G[2]*_K22) - _G[1]*_G[2]*_K12*_K12);
                _N[1] = - pow(pq(1), _J) * _n[1] * _G[1]*(_K01 + _G[2]*(_K01*_K22 - _K02*_K12));
                _N[2] = - pow(pq(2), _J) * _n[2] * _G[2]*(_K02 + _G[1]*(_K02*_K11 - _K01*_K12));

                complex T = (_N[0] + _N[1] + _N[2]) / _D;

                return T;
            };

            // Isolate the production amplitude (numerator) 
            inline double production(double Egam)
            {
                double w = W_cm(Egam);
                store(_kinematics->helicities(0), w*w, 0);
                evaluate();
                return std::real(_N[0] + _N[1] + _N[2]);
            };

            // Isolate the production amplitude (numerator) for a specific channel
            inline double production(int i, double Egam)
            {
                double w = W_cm(Egam);
                store(_kinematics->helicities(0), w*w, 0);
                evaluate();
                return std::real(_N[i]);
            };

            // Ratio of elastic to total production
            inline double inelasticity(double Egam)
            {
                double w = W_cm(Egam);
                store(_kinematics->helicities(0), w*w, 0);
                evaluate();

                double el = std::abs(_N[0]);
                double in = std::abs(_N[1] + _N[2]);

                return in / (el + in);
            };

            inline double scattering_length()
            {
                store(_kinematics->helicities(0), _kinematics->sth(), 0);
                evaluate();

                // Determinants of the sub-matrices
                complex detK01, detK02, detK12;
                detK01 = _K00*_K11 - _K01*_K01;
                detK02 = _K00*_K22 - _K02*_K02;
                detK12 = _K11*_K22 - _K12*_K12;

                // Use this to calculate the elastic T-matrix element
                double T00_th = std::real((_K00 + _G[1]*detK01 + _G[2]*detK02 + _G[1]*_G[2]*(detK12*_K00 + detK02*_K11 + detK01*_K22 - 2*_K00*_K11*_K22 + 2*_K01*_K02*_K12)) / _D);

                // Then calculate the scattering length
                double SL = - 2*T00_th / _kinematics->Wth();
                return SL / 5.068; // IN UNITS OF fm!!!!
            };

            inline double VMD_test()
            {
                store(_kinematics->helicities(0), _kinematics->sth(), 0);
                evaluate();
                
                double VMD  = 0.0273;
                double fit = std::abs(2*_n[0] / (_kinematics->Wth()*scattering_length()*5.068));
                return fit / VMD;
            };

            protected:

            inline void allocate_parameters(std::vector<double> pars)
            {
                switch (_Nth)
                {
                    default: return;
                    case 1: 
                    {
                        _n[0] = pars[0];
                        _a00  = pars[1];
                        break;
                    };
                    case 2:
                    {
                        _n[0] = pars[0];
                        _n[1] = pars[1];
                        _a00  = pars[2];
                        _a01  = pars[3];
                        _a11  = pars[4];
                        break;
                    };
                    case 3:
                    {
                        _n[0] = pars[0];
                        _n[1] = pars[1];
                        _n[2] = pars[2];
                        _a00  = pars[3];
                        _a01  = pars[4];
                        _a02  = pars[5];
                        _a11  = pars[6];
                        _a12  = pars[7];
                        _a22  = pars[8];
                        break;
                    };
                };

                if ( _option == EffectiveRange )
                {
                    _b00 = pars[2+(_Nth>1)*(4*_Nth-5)];
                    _b11 = (_Nth > 1)*pars[4*_Nth-2];
                    _b22 = (_Nth > 1)*pars[4*_Nth-1];
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
            double _a00 = 0, _a01 = 0, _a02 = 0;
            double           _a11 = 0, _a12 = 0;
            double                     _a22 = 0;

            // Second order parameters (effective range)
            double _b00 = 0, _b11 = 0, _b22 = 0;
            
            // Production amplitude parameters
            std::array<double,3> _n;

            // Internal variables for the K matrices
            complex _K00 = 0, _K01 = 0, _K02 = 0;
            complex           _K11 = 0, _K12 = 0;
            complex                     _K22 = 0;
            complex  _D = 0,  _detK = 0;

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