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

                if ( !are_equal(s, _cached_s, _cache_tolerance) ) recalculate();

                double rs = sqrt(s);
                print("sqrt_s =", rs);
                print("omega+(-sqrt_s) =", omega_f(+1, -rs));
                print("i*omega-(+sqrt_s) =", I*omega_f(-1, +rs));
                return 1;
            };

            protected:

            // Small imaginary part to help control cut structures
            complex _ieps = I*1E-10;

            // Since we calculate the full matrix every time, we cache and recalculate only if 
            // s changed.
            double _cached_s, _cache_tolerance = 1E-4;
            Eigen::Matrix2cd _Fp, _Fm;

            // Save the different parameters
            std::vector<double> _production_pars, _elastic_pars;
            inline void allocate_parameters(std::vector<double> pars)
            {
                _production_pars.clear(); _elastic_pars.clear();
                for (int i = 0; i < pars.size(); i++)
                {
                    if (i < 2) _production_pars.push_back(pars[i]); 
                    else       _elastic_pars.push_back(   pars[i]);
                };
            };
            
            // Initial fermion factor as function of sqrts
            inline complex omega_i(int pm, double sqrt_s)
            { 
                double Ei = (sqrt_s*sqrt_s + _mT*_mT - _mB*_mB)/(2*sqrt_s);
                return csqrt(Ei+pm*_mT+_ieps);
            };
            // Recoil fermion factor as function of sqrts
            inline complex omega_f(int pm, double sqrt_s)
            { 
                double Ef = (sqrt_s*sqrt_s + _mR*_mR - _mX*_mX)/(2*sqrt_s);
                return csqrt(Ef+pm*_mR+_ieps);
            };

            void recalculate()
            {

            };
        };
    };
};

#endif