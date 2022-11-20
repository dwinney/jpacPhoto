// Amplitude assuming a simple effective-range expansion near-threshold
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef SCATTERING_LENGTH
#define SCATTERING_LENGTH

#include "amplitude.hpp"

#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace jpacPhoto
{
    class scattering_length : public amplitude
    {
        // ---------------------------------------------------------------------------
        
        public:

        scattering_length(reaction_kinematics * xkinem, int lmax, std::string id = "scattering_length")
        : amplitude(xkinem, "scattering_length", id), _lmax(lmax)
        {
            set_nParams(_sc_pars);
            check_JP(xkinem);
        };

        // Set the free parameters
        inline void set_params(std::vector<double> params)
        {
            check_nParams(params);
    
            pars.clear();
            for (int n = 0; n < get_nParams(); n++)
            {
                pars.push_back(params[n]);
            }
        };

        inline void add_threshold(double m1, double m2)
        {
            if (_thresholds.size() == 2)
            {
                std::cout << "Only two additional thresholds can be added! Returning without change..." << std::endl;
                return;
            }

            // Save the masses of the new threshold
            _thresholds.push_back({m1, m2});

            // and add two to the number of free parameters
            // a_{0i} and r_i
            // set_nParams( get_nParams() + 2 );
        };  

        // Assemble the helicity amplitude by summing with appropriate angular funxtions
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double xs, double xt);

        // Expansion is directly in the S channel so helicities should always be defined in this frame
        helicity_channel helicity_CM_frame(){ return helicity_channel::S; };

        // only axial-vector, vector, and pseudo-scalar available
        inline std::vector<std::array<int,2>> allowed_meson_JP(){  return { {1, -1} }; };
        inline std::vector<std::array<int,2>> allowed_baryon_JP(){ return { {1, +1} }; };

        // ---------------------------------------------------------------------------

        protected:

        // Number of parameters in the single-channel amplitude
        static const int _sc_pars = 9;

        // Number of partial-wave terms to consider 
        int _lmax = 1;

        // Free parameters
        std::vector<double> pars;

        // A vector containing all open thresholds
        std::vector<std::array<double,2>> _thresholds; 

        // Legendre functions
        double P_l(int l, double z);
        inline std::complex<double> q2() { return Kallen(_s * XR, _mX*_mX * XR, _mR*_mR * XR) / (4.*_s); };
      
        // Redefine momenta here instead of using the reaction_kinematics versions
        // so that they can be appropriately analytically continued below threshold
        std::complex<double> barrier_factor(int l, double m1, double m2);
        inline std::complex<double> barrier_factor(int l){ return barrier_factor(l, _mX, _mR); };

        // Phase-space factor, this is defined with respect to a threshold to allow multiple thresholds to be considered
        std::complex<double> rho(double m1, double m2, double s);
        inline std::complex<double> rho(){ return rho(_mX, _mR, _s); };
        
        // Once-subtracted dispersion relation
        std::complex<double> rhoCM(double m1, double m2, double s);
        inline std::complex<double> rhoCM(std::array<double,2> x){ return rhoCM(x[0], x[1], _s); };
        inline std::complex<double> rhoCM(){ return rhoCM(_mX, _mR, _s); };
        
        // Partial wave amplitude
        inline std::complex<double> f_l(int l)
        {
            if (l == 0)
            {
                _rho0 = rhoCM();
                _rho1 = rhoCM(_thresholds[0]);
                _rho2 = rhoCM(_thresholds[1]);

                _a00 = pars[0];
                _a01 = pars[1];
                _a02 = pars[2];
                _a11 = pars[3];
                _a12 = pars[4];
                _a22 = pars[5];
                _b0  = pars[6];
                _c   = pars[7];
                _s0  = pars[8];

                // _denominator =  _a02*_a02*_rho1 + _a01*_a01*_rho2 + _a00*_rho1*_rho2 - _rho0*_rho1*_rho2;
                _denominator = 2.*_a01*_a02*_a12 - _a00*_a12*_a12 - _a01*_a01*_a22 + _a12*_a12*_rho0 - _a11*_a22*_rho0 + _a22*_rho0*_rho1 
                                + _a02*_a02*(_rho1-_a11) + _a00*(_a11-_rho1)*(_a22-_rho2) + _a01*_a01*_rho2 + _a11*_rho0*_rho2 - _rho0*_rho1*_rho2;

                return (_b0 + _rho0*T00() + _rho1*T10() + _rho2*T20()) ; 
            }

            return barrier_factor(l) / ( _c*pow(_s0, l) - barrier_factor(l)*rhoCM() );
        };  
    

        // Scattering length parameters
        double _a00 = 0., _a01 = 0., _a02 = 0., _a12 = 0., _a11 = 0., _a22 = 0., _b0 = 0., _b1 = 0., _c = 0., _s0 = 0.;
        double _r1, _r2;

        // Phase-space factors
        std::complex<double> _rho0, _rho1, _rho2;
        
        // Denominator of the T-matrixes
        std::complex<double> _denominator;

        // I'm going to hardcode the expressions for the three-channel T-matrix 
        // is there a better way to do this? probably but I'll refactor later
        inline std::complex<double> T00(){ return (_rho1 * _rho2 - _a11*_rho2 - _a22*_rho1 - _a12*_a12 + _a11*_a22) / _denominator; };
        inline std::complex<double> T10(){ return  (_a02*_a12 - _a01*_a22 + _a01 * _rho2) / _denominator; };
        inline std::complex<double> T20(){ return  (_a01*_a12 - _a02*_a11 + _a02 * _rho1) / _denominator; };
        
    };
};

#endif