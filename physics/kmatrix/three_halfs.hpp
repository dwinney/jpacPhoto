// Implementation of a HPWA with spin-3/2 intermediate state. This tries to be
// the most human readable by putting everything in 6x6 matrices but as a result
// its way too slow for fitting....
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

                return (!barred) ? (_FN(i)+_FU(i))/2 : (_FN(i)-_FU(i))/2; 
            };

            protected:

            // Save the different parameters
            // In the J=3/2 case we decompose the 6x6 problem into two 3x3's
            Eigen::Vector3cd _fp = Eigen::Vector3cd::Zero(), _fm = Eigen::Vector3cd::Zero();
            std::array<double,6> _a;
            inline void allocate_parameters(std::vector<double> pars)
            {

                for (int i = 0; i < 3; i++) _fp(i) = pars[i];
                for (int i = 0; i < 3; i++) _fm(i) = pars[3+i];
                for (int i = 0; i < 6; i++) _a[i]  = pars[6+i];
                _pars_changed = true;
            };
            inline std::vector<std::string> parameter_labels()
            {
                std::vector<std::string> labels;
                for (int i = 0; i < 6; i++) labels.push_back("f_3/2[" + to_string(i) +"]");
                for (int i = 0; i < 6; i++) labels.push_back("a_3/2[" + to_string(i) +"]");
                return labels;
            };

            // Small imaginary part to help control cut structures
            complex _ieps = I*1E-10;

            bool _pars_changed = false;
            double _cached_sqrts = -1, _cache_tolerance = 1E-4;

            // Partial wave vectors
            Eigen::VectorXcd _FN = Eigen::VectorXcd::Zero(6); // Natural 
            Eigen::VectorXcd _FU = Eigen::VectorXcd::Zero(6); // Unnatural

            // Recalculate cached partial-waves if needed
            void recalculate()
            {
                // Work in terms of \pm \sqrt{s}
                _cached_sqrts = sqrt(_s);
                
                Eigen::Matrix3cd One = Eigen::Matrix3cd::Identity(), one_minus_irhot;
                std::array<Eigen::Matrix3cd,4> m;
                
                // First we calculate the Natural case (\eta_s = +1 or \sqrt{s}>1)
                one_minus_irhot = (One+i_rho()*t_matrix(+1)); // Unitarized elastic amplitudes
                m = m_matrices(+1);                           // Production kinematic matrices
                _FN.head<3>() = one_minus_irhot*(m[0]*_fp+m[1]*_fm);
                _FN.tail<3>() = one_minus_irhot*(/* m[2]*_fp+ */m[3]*_fm);

                // Now we do the exact same but with \eta = -1
                one_minus_irhot = (One+i_rho()*t_matrix(-1));                   
                m = m_matrices(-1); 
                _FU.head<3>() = one_minus_irhot*(m[0]*_fp+m[1]*_fm);
                _FU.tail<3>() = one_minus_irhot*(/* m[2]*_fp+ */m[3]*_fm);

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
                complex omegafp = omega_f(+1, sqrt_s), omegafm = omega_f(-1, sqrt_s);

                Eigen::Matrix3cd n = Eigen::Matrix3cd::Zero();
                switch (i)
                {
                    case 0: 
                    {
                        n(1,1) = 1;
                        n *= pow(ps, 4) * sqrt_s*sqrt_s;
                        n *= omegafp * omegafp;
                        break;
                    };
                    case 1:
                    {
                        n(0,1) = 1; n(1,0) = 1;
                        n *= pow(ps, 3) * sqrt_s;
                        n *= omegafp * omegafm;
                        break;
                    };
                    case 2:
                    {
                        n(0,1) = 1;       n(1,0) = n(0,1);
                        n(1,2) = sqrt(3); n(2,1) = n(1,2);
                        n *= pow(ps, 2) * sqrt_s; 
                        n *= omegafp * omegafp;
                        break;
                    };
                    case 3:
                    {
                        n(0,0) = 2*mV*mV; 
                        n(1,1) = pow(mN+sqrt_s, 2); 
                        n(0,1) = sqrt(2)*mV*(mN+sqrt_s);  n(1,0) = n(0,1);  
                        n *= pow(ps, 2); 
                        n *= omegafm * omegafm;
                        break;
                    };
                    case 4:
                    {
                        double EV = (sqrt_s*sqrt_s+mV*mV-mN*mN)/2/sqrt_s;
                        n(0,1) = 2*EV-mN-sqrt_s;       n(1,0) = n(0,1);
                        n(0,2) = sqrt(6)*mV;           n(2,0) = n(0,2);
                        n(1,2) = -sqrt(3)*(mN+sqrt_s); n(2,1) = n(1,2);
                        n *= ps; 
                        n *= omegafp * omegafm;
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
                        n *= omegafp * omegafp;
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

            // Assemble the 3x3 (little) k-matrix
            Eigen::Matrix3cd k_matrix(int eta)
            {
                Eigen::Matrix3cd k = Eigen::Matrix3cd::Zero();
                for (int n = 0; n < 6; n++) k += n_matrix(n, eta) * _a[n];
                return k;
            };

            // Assemble the 3x3 t-matrix
            // i.e. the unitarized elastic rescattering amplitudes
            Eigen::Matrix3cd t_matrix(int eta)
            {
                Eigen::Matrix3cd One = Eigen::Matrix3cd::Identity();
                Eigen::Matrix3cd k   = k_matrix(eta);

                return (One-i_rho()*k).inverse()*k;
            };
            
            // Diagonal 3x3 blocks of the 6x6 (big) M-matrix
            std::array<Eigen::Matrix3cd,4> m_matrices(int eta)
            {
                double sqrt_s = (eta == +1) ? + _cached_sqrts : - _cached_sqrts;
                double mV = _mX, mN = M_PROTON;
                double qs = sqrt(kallen(sqrt_s*sqrt_s, 0.,    mN*mN))/2/sqrt_s;
                double ps = sqrt(kallen(sqrt_s*sqrt_s, mV*mV, mN*mN))/2/sqrt_s;
                double EV = (sqrt_s*sqrt_s+mV*mV-mN*mN)/2/sqrt_s;
                double x  = pow(mN-sqrt_s, 2) - mV*mV;
                double y  = (mN+sqrt_s);

                // m+ (i.e. the upper 3x3 diagonal block)
                Eigen::Matrix3cd mp = Eigen::Matrix3cd::Zero();
                mp(0,1) = sqrt(2)*mV*qs*x/sqrt_s;
                mp(0,2) = mV*qs;
                mp(1,0) = ps*ps*qs*sqrt_s;
                mp(1,1) = qs*x*y/sqrt_s;
                mp(1,2) = sqrt(2)*EV*qs;
                mp(2,2) = sqrt(3)*mV*qs;

                // m- (i.e. the lower 3x3 diagonal block)
                Eigen::Matrix3cd mm = Eigen::Matrix3cd::Zero();
                mm(0,0) = 0;
                mm(0,1) = sqrt(6)*mV*x/sqrt_s;
                mm(0,2) = sqrt(3)*mV;
                mm(1,0) = sqrt(3)*ps*ps*sqrt_s;
                mm(1,1) = sqrt(3)*x*y/sqrt_s;
                mm(1,2) = sqrt(6)*EV;
                mm(2,2) = 3*mV;

                // m'+ (i.e. the upper 3x3 off-diagonal block)
                Eigen::Matrix3cd mprimep = Eigen::Matrix3cd::Zero();
                mprimep(0,1) = sqrt(2)*mN*mV*x/sqrt_s/sqrt_s;
                mprimep(0,2) = mN*mV/sqrt_s;
                mprimep(1,0) = mN*ps*ps;
                mprimep(1,1) = mN*x*y/sqrt_s/sqrt_s;
                mprimep(1,2) = sqrt(2)*EV*mN/sqrt_s;
                mprimep(2,2) = sqrt(3)*mN*mV/sqrt_s;

                // m'- (i.e. the lower 3x3 off-diagonal block)
                Eigen::Matrix3cd mprimem = Eigen::Matrix3cd::Zero();

                complex prefactor = omega_i(+1, sqrt_s)*omega_f(+1, sqrt_s)*(mN-sqrt_s);
                mp *= prefactor; mprimep *= prefactor;
                mm *= prefactor; mprimem *= prefactor;
                
                if (eta == -1) 
                {
                    // MacDowell phase
                    for (int i = 0; i < 3; i++)
                    { 
                        // Multiply rows 2-4 of the 6x6 by -1
                        mp(2,i) *= -1;  mprimep(2,i) *= -1; 
                        mm(0,i) *= -1;  mprimem(0,i) *= -1;
                        mm(1,i) *= -1;  mprimem(1,i) *= -1; 
                    };
                };
                
                return {mp, mprimep, mprimem, mm};
            };
        };
    };
};

#endif