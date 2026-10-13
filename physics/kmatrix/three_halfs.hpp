// Implementation of a HPWA with spin-3/2 intermediate state
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

#ifndef THREE_HALFS_HPP
#define THREE_HALFS_HPP

#include "constants.hpp"
#include "kinematics.hpp"
#include "utilities.hpp"
#include "partial_wave.hpp"

#include <Eigen/Dense>

namespace jpacPhoto
{
    namespace kmatrix
    {
        class three_halfs : public raw_partial_wave
        {
            public: 

            // Single channel K-matrix (no arguments needed at initialization)
            three_halfs(key k, kinematics xkinem)
            : raw_partial_wave(k, xkinem, 3, "kmatrix::three_halfs")
            {
                // 6 Production parameters and 6 elastic rescattering
                initialize(12);
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
                bool barred = (2*hels[2] - hels[3] < 0);
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
                    // {+1, +1, +1, -1}
                    case 1: i = 2; break;
                    // {+1, -1, +1, +1}
                    case 6: i = 3; break;
                    // {+1, -1,  0, -1}
                    case 9: i = 4; break;
                    // {+1, -1, +1, -1}
                    case 7: i = 5; break;
                    default: return 0;
                };

                return (!barred) ? (_Fp(i)+_Fm(i))/2 : (_Fp(i)-_Fm(i))/2; 
            };

            protected:

            // Save the different parameters
            Eigen::VectorXcd _Q = Eigen::VectorXcd::Zero(6);
            std::array<double,6> _elastic_pars;
            inline void allocate_parameters(std::vector<double> pars)
            {
                for (int i = 0; i < pars.size(); i++)
                {
                    if (i < 6) _Q(i)               = pars[i]; 
                    else       _elastic_pars[i-6]  = pars[i];
                };
                _pars_changed = true;
            };

            // Small imaginary part to help control cut structures
            complex _ieps = I*1E-10;

            bool _pars_changed = false;
            double _cached_sqrts = -1, _cache_tolerance = 1E-4;

            // Partial wave matrices
            Eigen::VectorXcd _Fp = Eigen::VectorXcd::Zero(6), _Fm = Eigen::VectorXcd::Zero(6);

            // Recalculate cached partial-waves if needed
            void recalculate()
            {
                // Work in terms of \pm \sqrt{s}
                _cached_sqrts = sqrt(_s);
                
                Eigen::MatrixXcd One  = Eigen::MatrixXcd::Identity(6,6);
                // Unitarized K-matrix form for elastic amplitude
                Eigen::MatrixXcd K_p  = K_matrix(+1);
                Eigen::MatrixXcd T_p  = (One - i_rho()*K_p).inverse() *K_p;
                // Mixing matrix for production amplitude
                Eigen::MatrixXcd M_p  = M_matrix(+1);
                // Unitarized with the Q-vector of production parameters
                _Fp = (M_p + i_rho()*T_p*M_p)*_Q;

                // Now we do the exact same but with \eta = -1
                Eigen::MatrixXcd K_m  = K_matrix(-1);
                Eigen::MatrixXcd T_m  = (One - i_rho()*K_m).inverse() *K_m;
                Eigen::MatrixXcd M_m  = M_matrix(+1);
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

            // Calculate 3x3 (little) n-matrices
            Eigen::Matrix3cd n_matrix(uint i, int eta)
            {              
                double sqrt_s = (eta == +1) ? + _cached_sqrts : - _cached_sqrts;
                double mV = _mX, mN = M_PROTON;
                double ps = sqrt(kallen(sqrt_s*sqrt_s, mV*mV, mN*mN))/2/sqrt_s;
                
                Eigen::Matrix3cd n = Eigen::Matrix3cd::Zero();
                switch (i)
                {
                    case 0: 
                    {
                        n(1,1) = 1;
                        n *= pow(ps, 4) * sqrt_s*sqrt_s;
                        n *= omega_f(+1, sqrt_s) * omega_f(+1, sqrt_s);
                        break;
                    };
                    case 1:
                    {
                        n(0,1) = 1; n(1,0) = 1;
                        n *= pow(ps, 3) * sqrt_s;
                        n *= omega_f(+1, sqrt_s) * omega_f(-1, sqrt_s);
                        break;
                    };
                    case 2:
                    {
                        n(0,1) = 1;       n(1,0) = n(0,1);
                        n(1,2) = sqrt(3); n(2,1) = n(1,2);
                        n *= pow(ps, 2) * sqrt_s; 
                        n *= omega_f(+1, sqrt_s) * omega_f(+1, sqrt_s);
                        break;
                    };
                    case 3:
                    {
                        n(0,0) = 2*mV*mV; 
                        n(1,1) = pow(mN+sqrt_s, 2); 
                        n(0,1) = sqrt(2)*mV*(mN+sqrt_s);  n(1,0) = n(0,1);  
                        n *= pow(ps, 2); 
                        n *= omega_f(-1, sqrt_s) * omega_f(-1, sqrt_s);
                        break;
                    };
                    case 4:
                    {
                        double EV = (sqrt_s*sqrt_s+mV*mV-mN*mN)/2/sqrt_s;
                        n(0,1) = 2*EV-mN-sqrt_s;       n(1,0) = n(0,1);
                        n(0,2) = sqrt(6)*mV;           n(2,0) = n(0,2);
                        n(1,2) = -sqrt(3)*(mN+sqrt_s); n(2,1) = n(1,2);
                        n *= ps; 
                        n *= omega_f(+1, sqrt_s) * omega_f(-1, sqrt_s);
                        break;
                    };
                    case 5:
                    {
                        double EV = (sqrt_s*sqrt_s+mV*mV-mN*mN)/2/sqrt_s;
                        n(0,0) = mV*mV;
                        n(1,1) = 2*EV*EV;
                        n(2,2) = 3*mV*mV;

                        n(0,1) = sqrt(2)*EV*mV;   n(1,0) = n(0,1);
                        n(0,2) = sqrt(3)*mV*mV;   n(2,0) = n(0,2);
                        n(1,2) = sqrt(6)*EV*mV;   n(2,1) = n(1,2);
                        n *= omega_f(+1, sqrt_s) * omega_f(+1, sqrt_s);
                        break;
                    };
                    default: return n;
                };
                
                // If unnatural, multiply by (-1)^{\lam-lamp}
                if (eta == -1)
                {
                    n(0,2) *= -1; n(2,0) *= -1;
                    n(1,2) *= -1; n(2,1) *= -1;
                };

                return n;
            };

            // Assemble the 3x3 (little) k-matrix into a 6x6 block diagonal (big) K-matrix
            Eigen::MatrixXcd K_matrix(int eta)
            {
                Eigen::Matrix3cd k = Eigen::Matrix3cd::Zero();
                for (int n = 0; n < 6; n++) k += n_matrix(n, eta) * _elastic_pars[n];

                Eigen::MatrixXcd K = Eigen::MatrixXcd::Zero(6,6);
                K.block<3,3>(0,0) = k; K.block<3,3>(3,3) = k;
                return K;
            };
            
            // Mixing matrix between helicities and the 2 possible helicity structures
            Eigen::MatrixXcd M_matrix(int eta)
            {
                double sqrt_s = (eta == +1) ? + _cached_sqrts : - _cached_sqrts;
                double mV = _mX, mN = M_PROTON;
                double qs = sqrt(kallen(sqrt_s*sqrt_s, 0.,    mN*mN))/2/sqrt_s;
                double ps = sqrt(kallen(sqrt_s*sqrt_s, mV*mV, mN*mN))/2/sqrt_s;
                double EV = (sqrt_s*sqrt_s+mV*mV-mN*mN)/2/sqrt_s;

                Eigen::MatrixXcd A = Eigen::MatrixXcd::Zero(6,6);
                A(0,0) = ps*ps*qs*pow(sqrt_s,3);
                A(1,1) = qs*sqrt_s*(pow(mN-sqrt_s, 2) - mV*mV);
                A(2,2) = qs*sqrt_s*sqrt_s;
                A(3,3) = ps*ps*sqrt_s*sqrt_s;
                A(4,4) = (pow(mN-sqrt_s, 2) - mV*mV);
                A(5,5) = sqrt_s;

                Eigen::MatrixXcd C = Eigen::MatrixXcd::Zero(6,6);
                C(0,0) = mV;
                C(1,1) = 1;
                C(2,2) = sqrt(3)*mV;
                C(3,3) = sqrt(3)*mV*sqrt_s;
                C(4,4) = sqrt(3)*sqrt_s;
                C(5,5) = 3*mV*sqrt_s;

                Eigen::MatrixXcd B = Eigen::MatrixXcd::Zero(6,6);
                B(0,1) = 1;
                B(1,0) = sqrt(2); B(1,1) = (mN+sqrt_s);
                B(2,0) = 1; B(2,1) = sqrt(2)*EV; B(2,2) = 1;
                B(3,1) = mN; B(3,4) = 1;
                B(4,0) = sqrt(2)*mN; B(4,1) = mN*(mN+sqrt_s); B(4,3) = sqrt(2); B(4,4) = (mN+sqrt_s);
                B(5,0) = mN; B(5,1) = sqrt(2)*EV*mN; B(5,2) = mN; B(5,3) = 1; B(5,4) = sqrt(2)*EV; B(5,5) = 1; 

                Eigen::MatrixXcd M = (A*B*C).adjoint();
                M *= omega_i(+1, sqrt_s) * omega_f(+1, sqrt_s) * (mN - sqrt_s)/sqrt_s/sqrt_s;

                if (eta == -1) 
                {
                    // MacDowell phase
                    for (int i = 0; i < 6; i++){ M(2,i) *= -1; M(3,i) *= -1; M(4,i) *= -1; };
                };
                
                return M;
            };
        };
    };
};

#endif