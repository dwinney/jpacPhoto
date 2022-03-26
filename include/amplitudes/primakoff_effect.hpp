// Axial-vector meson photoproduction proceeding through the primakoff effect
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _PRIMAKOFF_
#define _PRIMAKOFF_

#include "amplitude.hpp"

namespace jpacPhoto
{
    class primakoff_effect : public amplitude
    {
        public:
        // Constructor 
        primakoff_effect(reaction_kinematics * xkinem, std::string amp_id = "primakoff_effect")
        : amplitude(xkinem, "primakoff_effect", amp_id)
        {
            set_nParams(4);
            check_JP(xkinem->_jp);
        };

        void set_params(std::vector<double> params)
        {
            check_nParams(params); 
            _atomicZ        = params[0];
            _atomicRadius   = params[1];
            _skinThickness  = params[2];
            _photonCoupling = params[3];

            calculate_norm();
        };

        inline void set_LT(int LT)
        {
            if (LT > 1 || LT < 0)
            {
                std::cout << "error! invalid parameter in set_LT(). \n";
                std::cout << "LT = 0 for longitudinal and 1 for transverse photon.\n";
            };

            _helProj = LT;
        };

        // individual helicity amplitudes not supported but need to provide definition for virtual class.
        inline std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t)
        {
            std::cout << "Warning! Individual helicity amplitudes not supported by primakoff_effect!\n";
            return 0.;
        }

        // instead we override the definition of differential_xsection in amplitude.hpp
        double differential_xsection(double s, double t);
        double integrated_xsection(double s);

        // only axial-vector kinematics allowed
        inline std::vector<std::array<int,2>> allowedJP()
        {
            return {{1, 1}};
        };

        private:

        // Parameters
        int    _helProj              = 0 ;  // longitudinal (0) or transverse (1) photon
        int    _atomicZ              = 0 ;  // atomic number
        double _atomicRadius         = 0.;  // radius parameter
        double _skinThickness        = 0.;  // skin thickness parameter
        double _photonCoupling       = 0.;  // X -> gamma gamma* coupling

        // Fermi model nuclear charge distribution
        inline double charge_distribution(double r)
        {
            return 1. / ( 1. + exp((r - _atomicRadius) / _skinThickness) );
        };

        // Normalized fourier transform of the above charge_distributions 
        double form_factor(double x);
        double _formFactor; // Form factor at energy t
        
        void calculate_norm();     
        double _rho0 = 0.;  // normalizaton
        inline double W_00()
        {
            return 64. * _atomicZ*_atomicZ * _mA2 * _mA2 * _mA2 * _formFactor * _formFactor / ((_t - 4.*_mA2) * (_t - 4.*_mA2));
        };

        // Kinematic quantities   
        long double _mX2 =  _kinematics->_mX2;
        long double _mA2 =  _kinematics->_mT2;
        long double _mQ2  = -_kinematics->_mB2;

        long double _cosX, _sinX2;
        long double _pGam, _pX;
        long double _nu, _enX;
        inline void update_kinematics()
        {
            // lab frame momentum transfer
            _nu = (_s - _mA2 + _mQ2) / (2. * sqrt(_mA2));

            // Momentum of photon
            _pGam = sqrt(_nu*_nu + _mQ2);

            // // Momentum of the X
            _pX  = sqrt(_t*_t + 4.*sqrt(_mA2)*_t*_nu + 4.*_mA2*(_nu*_nu - _mX2));
            _pX /= 2. * sqrt(_mA2);

            // Energy of the X
            _enX = sqrt(_pX*_pX + _mX2);

            // Cosine of scattering angle of the X in the lab frame
            _cosX  = _t + _mQ2 - _mX2 + 2.*_nu*_enX;
            _cosX /= 2. * _pX * _pGam;

            // // Sine of the above 
            _sinX2 = 1. - _cosX * _cosX;
        };

        // Spin summed amplitude squared
        long double amplitude_squared();
    };
};

#endif