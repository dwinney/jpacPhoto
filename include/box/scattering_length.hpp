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

        scattering_length(reaction_kinematics * xkinem, std::string id = "scattering_length")
        : amplitude(xkinem, "scattering_length", id)
        {
            set_Lmax(0);
            check_JP(xkinem);
        };

        scattering_length(reaction_kinematics * xkinem, int lmax, std::string id = "scattering_length")
        : amplitude(xkinem, "scattering_length", id)
        {
            set_Lmax(lmax);
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

        inline void set_Lmax(int l)
        {
            _lmax = l;
            set_nParams(2+2*_lmax);
        }

        // Assemble the helicity amplitude by summing with appropriate angular funxtions
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double xs, double xt);

        // Expansion is directly in the S channel so helicities should always be defined in this frame
        helicity_channel helicity_CM_frame(){ return helicity_channel::S; };

        // only vector and proton available
        inline std::vector<std::array<int,2>> allowed_meson_JP(){  return { {1, -1} }; };
        inline std::vector<std::array<int,2>> allowed_baryon_JP(){ return { {1, +1} }; };

        // ---------------------------------------------------------------------------

        protected:

        // Number of partial-wave terms to consider 
        int _lmax = 1;

        // Free parameters
        std::vector<double> pars;

        inline double q2(){ return Kallen(_s, _mX*_mX, _mR*_mR ) / (4.*_s); };
      
        // Redefine momenta here instead of using the reaction_kinematics versions
        // so that they can be appropriately analytically continued below threshold
        std::complex<double> barrier_factor(int l, double m1, double m2);
        inline std::complex<double> barrier_factor(int l){ return barrier_factor(l, _mX, _mR); };
        
        // Phase-space factor, this is defined with respect to a threshold to allow multiple thresholds to be considered
        std::complex<double> rho(double m1, double m2, double s);
        inline std::complex<double> rho(){ return rho(_mX, _mR, _s); };
        
        // Once-subtracted dispersion relation (Chew-Mandelstam)
        std::complex<double> rhoCM(double m1, double m2, double s);
        inline std::complex<double> rhoCM(std::array<double,2> x){ return rhoCM(x[0], x[1], _s); };
        inline std::complex<double> rhoCM(){ return rhoCM(_mX, _mR, _s); };

        inline std::complex<double> A_L(int l)
        {
            double a, b;

            if ((2*l) > pars.size()) return 0.;
            a = pars[2*l]; b = pars[2*l+1];

            std::complex<double> K = pow(q2(), l) * a;
            std::complex<double> T = K / (1. -  rhoCM() * K);

            return barrier_factor(l) * b * (1. - rhoCM()*T);
        };
    };
};

#endif        