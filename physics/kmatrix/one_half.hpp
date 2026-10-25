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
        using matrix = Eigen::MatrixXcd;

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

                auto hels = helicities;
                bool barred = (2*helicities[2] - helicities[3] < 0);
                if (barred) { hels[2] *= -1; hels[3] *= -1; };

                // Many years ago, I didnt think very hard about the order the helicities should be arranged
                // Because of that we need to remap them to match our matrices....
                // is there a more elegant way to do this? probably but idc
                int h = _kinematics->helicity_index(hels), i;
                switch (h)
                {   
                    // {+1, +1, +1, +1}
                    case 0: i = 0; break;
                    // {+1, +1,  0, -1}
                    case 3: i = 1; break;
                    default: return 0;
                };

                return (!barred) ? (_Fp(i)+_Fm(i))/2 : (_Fp(i)-_Fm(i))/2; 
            };

            protected:

            // Save the different parameters
            matrix _Q = Eigen::Vector2cd::Zero();
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
            inline std::vector<std::string> parameter_labels()
            {
                return {"f_1/2[0]", "f_1/2[1]", 
                        "a_1/2[0]", "a_1/2[1]", "a_1/2[2]"};
            };

            // Small imaginary part to help control cut structures
            complex _ieps = I*1E-10;

            bool _pars_changed = false;
            double _cached_sqrts = -1, _cache_tolerance = 1E-4;

            // Partial wave matrices
            matrix _Fp = Eigen::Vector2cd::Zero(), _Fm = Eigen::Vector2cd::Zero();

            // Recalculate cached partial-waves if needed
            void recalculate()
            {
                // Work in terms of \pm \sqrt{s}
                _cached_sqrts = sqrt(_s);

                matrix One  = Eigen::Matrix2cd::Identity();
                // Unitarized K-matrix form for elastic amplitude
                matrix K_p  = K_matrix(+1);
                matrix T_p  = (One - i_rho()*K_p).inverse() *K_p;
                // Mixing matrix for production amplitude
                matrix M_p  = M_matrix(+1);
                // Unitarized with the Q-vector of production parameters
                _Fp = (M_p + i_rho()*T_p*M_p)*_Q;

                // Now we do the exact same but with \eta = -1
                matrix K_m  = K_matrix(-1);
                matrix T_m  = (One - i_rho()*K_m).inverse() *K_m;
                matrix M_m  = M_matrix(-1);
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

            // Calculate 2x2 N-matrices
            inline matrix N_matrix(uint n, int eta)
            {              
                double sqrt_s = (eta == +1) ? + _cached_sqrts : - _cached_sqrts;
                matrix N = Eigen::Matrix2cd::Zero();

                // Factors of (2 p_s \sqrt{s}) 
                double ps = sqrt(kallen(sqrt_s*sqrt_s,_mX*_mX, _mR*_mR))/2/sqrt_s;

                switch (n)
                {
                    case 0: 
                    {
                        N(1,1) = 1;
                        N *= pow(ps*sqrt_s, 2);
                        N *= omega_f(+1, sqrt_s) * omega_f(+1, sqrt_s);
                        return N;
                    };
                    case 1:
                    {
                        N(0,1) = 1; N(1,0) = 1;
                        N *= pow(ps*sqrt_s, 2);
                        N *= omega_f(+1, sqrt_s) * omega_f(-1, sqrt_s);
                        return N;
                    };
                    case 2:
                    {
                        double ep = M_PROTON + sqrt_s;
                        N(0,0) = 2*_mX*_mX;       N(0,1) = sqrt(2)*_mX*ep;
                        N(1,0) = sqrt(2)*_mX*ep;  N(1,1) = ep*ep;   
                        N *= omega_f(-1, sqrt_s) * omega_f(-1, sqrt_s);
                        return N;
                    };
                    default: return N;
                    // No eta dependence because macdowell phase is trivial
                };
            };

            // 2x2 K-matrix is sum of N's multiplied by the scattering lengths
            inline matrix K_matrix(int eta)
            {
                matrix K = Eigen::Matrix2cd::Zero();
                for (int n = 0; n < 3; n++) K += N_matrix(n, eta) * _elastic_pars[n];
                return K;
            };
            
            // Mixing matrix between helicities and the 2 possible helicity structures
            inline matrix M_matrix(int eta)
            {
                double sqrt_s = (eta == +1) ? + _cached_sqrts : - _cached_sqrts;

                matrix M = Eigen::Matrix2cd::Zero();
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