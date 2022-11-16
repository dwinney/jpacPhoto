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
            set_nParams(NumPars);
            check_JP(xkinem);
        };

        // Set the free parameters
        inline void set_params(std::vector<double> params)
        {
            check_nParams(params);
    
            pars.clear();
            for (int n = 0; n < get_nParams() - _extra_thresholds.size(); n++)
            {
                pars.push_back(params[n]);
            }

            _extra_couplings.clear();
            for (int n_th = 0; n_th < _extra_thresholds.size(); n_th++)
            {
                _extra_couplings.push_back(params[NumPars+n_th]);
            };
        };

        inline void add_threshold(double m1, double m2)
        {
            // Save the masses of the new threshold
            _extra_thresholds.push_back({m1, m2});

            // and add one to the number of free parameters
            set_nParams( get_nParams() + 1 );
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

        // Number of partial-wave terms to consider 
        int _lmax = 1;

        // Free parameters
        std::vector<double> pars;

        // A vector containing all open thresholds
        std::vector<std::array<double,2>> _extra_thresholds; 
        std::vector<double>               _extra_couplings;

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
        inline std::complex<double> rhoCM(){ return rhoCM(_mX, _mR, _s); };

        std::complex<double> inelastic();
        
        static const int NumPars = 3;
        // Partial wave amplitude
        inline std::complex<double> f_l(int l)
        {
            std::complex<double> P, Kinv;
            std::complex<double> BF = barrier_factor(l);

            double a, b, s0;
            a  = pars[0];
            b  = pars[1];
            s0 = pars[2];

            switch (l)
            {
                case 0:
                {
                    P =  b;
                    Kinv = -1./a + inelastic();
                    break;
                }
                default:
                {
                    P = b;
                    Kinv = -pow(s0, l) / a;
                    break;
                };
            }

            return BF * P / ( Kinv - BF * rhoCM() );
        };
    };
};

#endif