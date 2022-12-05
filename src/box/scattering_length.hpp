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
            set_nParams(2*(1+_lmax) + 3);
        };

        inline void add_threshold(double m1, double m2)
        {
            // only allow 2 thresholds
            if (_thresholds.size() > 2) return;

            _thresholds.push_back({m1,m2});
            set_nParams(2*(1+_lmax) + 3);
        };

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

        std::vector<std::array<double,2>> _thresholds;

        inline double q2(){ return Kallen(_s, _mX*_mX, _mR*_mR ) / (4.*_s); };
      
        // Redefine momenta here instead of using the reaction_kinematics versions
        // so that they can be appropriately analytically continued below threshold
        std::complex<double> barrier_factor(int l, double m1, double m2);
        inline std::complex<double> barrier_factor(int l){ return barrier_factor(l, _mX, _mR); };
        
        // Phase-space factor, this is defined with respect to a threshold to allow multiple thresholds to be considered
        std::complex<double> rho(double m1, double m2, double s);
        inline std::complex<double> rho(){ return rho(_mX, _mR, _s); };
        
        // Once-subtracted dispersion relation (Chew-Mandelstam) subtracted at (m1 + m2)^2
        std::complex<double> rhoCM(double m1, double m2, double s);

        // with no arguments this is the primary two-body threshold saves in _kinematics
        inline std::complex<double> rhoCM(){ return rhoCM(_mX, _mR, _s); };

        // For arbitrary extra open threshold we subtract at the main (lowest) threshold
        inline std::complex<double> rhoCM(std::array<double,2> x){ return rhoCM(x[0], x[1], _s) - rhoCM(x[0], x[1], _kinematics->sth()); };

        inline std::complex<double> A_L(int l)
        {
            double a, b;

            if (l > 0 || _thresholds.size() == 0)
            {
                if ((2*l) > pars.size()) return 0.;
                a = pars[2*l]; b = pars[2*l+1];

                std::complex<double> K = pow(q2(), l) * a;
                std::complex<double> T = K / (1. -  rhoCM() * K);

                return barrier_factor(l) * b * (1. - rhoCM()*T);
            }

            return coupled_swave();
        };

        inline std::complex<double> coupled_swave()
        {
            // three different phase-space factors
            std::complex<double> rho0, rho1, rho2;  
            rho0 = rhoCM();
            rho1 = (_thresholds.size() > 0) ? rhoCM(_thresholds[0]) : 1.;
            rho2 = (_thresholds.size() > 1) ? rhoCM(_thresholds[1]) : 1.;

            // Couplings from the 3x3 coupled-channel K-matrix
            double a00 = 0., a01 = 0., a02 = 0.;
            double a11 = 0., a12 = 0., a22 = 0.;

            // Normalization / scattering length
            double b = 0.;
            double r = 0.;

            // a00 comes from the from of the parameter vector 
            a00 = pars[0]; b = pars[1];

            // The rest from the back
            // At first we only assume only direct 01 and 02 parameters are important
            a01 = pars[(2*_lmax+1) + 1];
            a02 = pars[(2*_lmax+1) + 2];
            r   = pars[(2*_lmax+1) + 3];
            // a11 = pars[(2*_lmax+1) + 3];

            // Denominator common to all amplitudes
            std::complex<double> denominator = 2.*a01*a02*a12 - a00*a12*a12 - a01*a01*a22 + a12*a12*rho0 - a11*a22*rho0 + a22*rho0*rho1 
                                                        + a02*a02*(rho1-a11) + a00*(a11-rho1)*(a22-rho2) + a01*a01*rho2 + a11*rho0*rho2
                                                                                                                       - rho0*rho1*rho2;
            
            std::complex<double> T00 = 0., T01 = 0., T02 = 0.;
            T00 = (rho1*rho2 - a11*rho2 - a22*rho1 - a12*a12 + a11*a22) / denominator;
            T01 = (a02*a12 - a01*a22 + a01*rho2) / denominator;
            T02 = (a01*a12 - a02*a11 + a02*rho1) / denominator;

            return b * (1. - rho0*T00) - r*(rho1*T01 - rho2*T02);
        };
    };
};

#endif        