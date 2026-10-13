// Implementation of a HPWA with spin-1/2 intermediate state
//
// This amplitude explicitly uses Eigen C++ to manipulate complex matrices,
// make sure an environment variable EIGEN points to the top level directory
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Instituto de Ciencias Nuclears (ICN),
//               Universidad Nacional Autonoma de México (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#ifndef ONE_HALF_HPP
#define ONE_HALF_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "utilities.hpp"
#include "partial_wave.hpp"

#include <Eigen/Dense>

namespace jpacPhoto
{
    namespace kmatrix
    {
        class one_half : public raw_partial_wave
        {
            public: 

            // Single channel K-matrix (no arguments needed at initialization)
            one_half(key k, kinematics xkinem)
            : raw_partial_wave(k, xkinem, 1, "kmatrix::one_half")
            {
                // 2 Production parameters and 3 elastic rescattering
                initialize(5);
            };

            // -----------------------------------------------------------------------
            // Virtuals 

            // Specify vector quantum numbers
            inline std::vector<quantum_numbers> allowed_mesons() { return {VECTOR};   };
            // and proton target
            inline std::vector<quantum_numbers> allowed_baryons(){ return {HALFPLUS}; };
            // We have explicit helicity dependence in s-channel
            inline helicity_frame native_helicity_frame(){ return S_CHANNEL; };

            // Partial wave comes from K-matrix unitarized form
            inline complex partial_wave(std::array<int,4> helicities, double s)
            {
                store(helicities, s, _t);

                bool need_to_recalculate =  !are_equal(sqrt(_s), _cached_sqrts, _cache_tolerance) || _pars_changed;
                if ( need_to_recalculate ) recalculate();

                // Many years ago, I didnt think very hard about the order the helicities should be arranged
                // Because of that we need to remap them to match our matrices....
                int hi = _kinematics->helicity_index(helicities);
                switch (hi)
                {   
                    // {+1, +1, +1, +1}
                    case 0: return _Fp(0) + _Fm(0);
                    // {+1, +1, -1, -1}
                    case 5: return _Fp(0) - _Fm(0);
                    // {+1, +1,  0, -1}
                    case 3: return _Fp(1) + _Fm(1);
                    // {+1, +1,  0, +1}
                    case 2: return _Fp(1) - _Fp(1);
                    default: return 0;
                };
            };

            protected:

            // Save the different parameters
            Eigen::Vector2cd _Q;
            std::array<double,3> _elastic_pars;
            inline void allocate_parameters(std::vector<double> pars)
            {
                for (int i = 0; i < pars.size(); i++)
                {
                    if (i < 2) _Q(i)               = pars[i]; 
                    else       _elastic_pars[i-2]  = pars[i];
                };
                _pars_changed = true;
            };

            // Small imaginary part to help control cut structures
            complex _ieps = I*1E-10;

            bool _pars_changed = false;
            double _cached_sqrts = -1, _cache_tolerance = 1E-4;

            // Partial wave matrices
            Eigen::Vector2cd _Fp, _Fm;

            // Recalculate cached partial-waves if needed
            void recalculate()
            {
                // Work in terms of \pm \sqrt{s}
                _cached_sqrts = sqrt(_s);
                
                Eigen::Matrix2cd One  = Eigen::Matrix2cd::Identity();
                // Unitarized K-matrix form for elastic amplitude
                Eigen::Matrix2cd K_p  = k_matrix(+1);
                Eigen::Matrix2cd T_p  = (One - i_rho()*K_p).inverse() *K_p;
                // Mixing matrix for production amplitude
                Eigen::Matrix2cd M_p  = mixing_matrix(+1);
                // Unitarized with the Q-vector of production parameters
                _Fp = (M_p + i_rho()*T_p*M_p)*_Q;

                // Now we do the exact same but with \eta = -1
                Eigen::Matrix2cd K_m  = k_matrix(-1);
                Eigen::Matrix2cd T_m  = (One - i_rho()*K_m).inverse() *K_m;
                Eigen::Matrix2cd M_m  = mixing_matrix(+1);
                _Fm = (M_m + i_rho()*T_m*M_m)*_Q;
            };  

            // Initial fermion factor as function of sqrts
            // NOTE: MASS OF PROTON IS HARD-CODED HERE
            inline complex omega_i(int pm, double sqrt_s)
            { 
                double Ei = (sqrt_s*sqrt_s + M_PROTON*M_PROTON)/(2*sqrt_s);
                return csqrt(Ei+pm*M_PROTON+_ieps);
            };
            // Recoil fermion factor as function of sqrts
            // NOTE: MASS OF PROTON IS HARD-CODED HERE
            inline complex omega_f(int pm, double sqrt_s)
            { 
                double Ef = (sqrt_s*sqrt_s + M_PROTON*M_PROTON - _mX*_mX)/(2*sqrt_s);
                return csqrt(Ef+pm*M_PROTON+_ieps);
            };

            // Calculate 2x2 K-matrices
            Eigen::Matrix2cd k_matrix(uint n, int eta)
            {              
                double sqrt_s = (eta == +1) ? + _cached_sqrts : - _cached_sqrts;
                Eigen::Matrix2cd K = Eigen::Matrix2cd::Zero();

                // Factors of (2 p_s \sqrt{s}) 
                double lamf = sqrt(kallen(sqrt_s*sqrt_s,_mX*_mX, _mR*_mR));
               
                switch (n)
                {
                    case 0: 
                    {
                        K(1,1) = 1;
                        K *= lamf *lamf;
                        K *= omega_f(+1, sqrt_s) * omega_f(+1, sqrt_s);
                        return K;
                    };
                    case 1:
                    {
                        K(0,1) = 1; K(1,0) = 1;
                        K *= lamf;
                        K *= omega_f(+1, sqrt_s) * omega_f(-1, sqrt_s);
                        return K;
                    };
                    case 2:
                    {
                        double ep = M_PROTON + sqrt_s;
                        K(0,0) = 2*_mX*_mX;       K(0,1) = sqrt(2)*_mX*ep;
                        K(1,0) = sqrt(2)*_mX*ep;  K(1,1) = ep*ep;   
                        K *= omega_f(-1, sqrt_s) * omega_f(-1, sqrt_s);
                        return K;
                    };
                    default: return K;
                    // No eta dependence because macdowell phase is trivial
                };
            };
            Eigen::Matrix2cd k_matrix(int eta)
            {
                Eigen::Matrix2cd K = Eigen::Matrix2cd::Zero();
                for (int n = 0; n < 3; n++) K += k_matrix(n, eta) * _elastic_pars[n];
                return K;
            };
            
            // Mixing matrix between helicities and the 2 possible helicity structures
            Eigen::Matrix2cd mixing_matrix(int eta)
            {
                double sqrt_s = (eta == +1) ? + _cached_sqrts : - _cached_sqrts;

                Eigen::Matrix2cd M = Eigen::Matrix2cd::Zero();
                double ep = M_PROTON + sqrt_s, em = M_PROTON - sqrt_s;

                M(0,0) = 0;                M(0,1) = 2*sqrt(2)*_mX;
                M(1,0) = ep*ep - _mX*_mX;  M(1,1) = 2*ep;
                M *= omega_i(+1, sqrt_s) * omega_f(-1, sqrt_s) * em;

                // No eta dependence because macdowell phase is trivial
                return M;
            };
        };
    };
};

#endif