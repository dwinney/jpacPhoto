// Phenomenological expressions for the total cross-sections using SAID PW's at low energies
// and Regge poles at high energies
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef SIGMA_TOT_JPAC
#define SIGMA_TOT_JPAC

#include "total_xsection.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Total cross-sections based off of the 2015 JPAC publication
    // V. Mathieu et al. [Phys. Rev. D 92, 074004]
    class JPAC_parameterization : public total_xsection
    {
        public:

        // Currently only avialable for pion-nucleon 
        JPAC_parameterization(int iso, bool resonances = true)
        : total_xsection(M_PION, M_PROTON, 6.5), _iso(iso),
          _interp(0, ROOT::Math::Interpolation::kCSPLINE)
        {
            if (abs(iso) != 1)
            {
                std::cout << "Error! JPAC_parameterization argument must be +1 or -1! Results may vary..." << std::endl;
            };

            if (resonances)
            {
                import_data();
            }
            else
            {
                _cutoff = 0.;
            }
        };

        protected:
        // ----------------------------------------------------------------------
        // Low-energy pieces

        // Below the cutoff we do a simple interpolation of the SAID amplitude
        inline double resonances(double s, double q2)
        {
            return _interp.Eval( pLab(s) );
        };
        
        class SAID_PW
        {
            SAID_PW(double L, double I, double J)
            {
                import_data();
            }
            
            // Read in the data file and interpolate the imaginary part of the amplitude
            void import_data();
            ROOT::Math::Interpolator _interp;
            std::vector<double> _plab, _imagPart;

            double _L, _I, _J; // 2 x quantum number
        };


        // If theres data available in the low enegergy, 
        // these methods read data files and store them
        ROOT::Math::Interpolator _interp;
        std::vector<double> _plab, _sigma;
        void import_data(std::string datfile = "SAID.dat");

        // ----------------------------------------------------------------------
        // High-energy pieces

        // Above the cutoff, we use a Regge pole parameterization
        inline double regge(double s)
        {
            return 0.389352 * std::imag(isoscalar(s, 0.) + double(_iso) * isovector(s, 0.)) / pLab(s);
        };

        // Parameter selecting pi + or pi - scattering
        int _iso;  // plus or minus 1

        // Internal class object for a Regge pole
        class regge_pole
        {
            public: 

            regge_pole(std::array<double, 2> pars, regge_trajectory * traj)
            : _c(pars[0]), _b(pars[1]), _alpha(traj)
            {};

            inline std::complex<double> operator()(double s, double t)
            {
                double nu = crossing_nu(s, t);
                double alpha = std::real(_alpha->eval(t));
                std::complex<double> signature_factor = 0.5 * (double(_alpha->_signature) - exp(- XI * M_PI * alpha));

                return _c * exp(_b * t) * signature_factor * pow(nu, alpha) * cgamma(1. - double(_alpha->_minJ) - alpha);
            };

            private:

            // Crossing variable
            inline double crossing_nu(double s, double t)
            {
                double u = 2.* (M_PROTON + M_PION) - s - t;
                return (s - u) / (4.* M_PROTON);
            };

            double _b; // t-slope
            double _c; // Overall coupling

            // Trajectory object
            regge_trajectory * _alpha;
        };

        // Trajectories for dominant Regge poles: Pomeron, F2, and Rho
        linear_trajectory alphaPom = linear_trajectory(-1, 0, 1.075, 0.434, "Pomeron");
        linear_trajectory alphaF   = linear_trajectory(-1, 0, 0.490, 0.943, "f2");
        linear_trajectory alphaRho = linear_trajectory(+1, 1, 0.490, 0.943, "rho");

        // Instances of the poles themselves
        regge_pole _f    = regge_pole( {71.35, 3.18}, &alphaF);
        regge_pole _pom  = regge_pole( {23.89, 2.21}, &alphaPom);
        regge_pole _rho1 = regge_pole( { 5.01*1.573, 10.1}, &alphaRho);
        regge_pole _rho2 = regge_pole( {-5.01*0.573,  0.0}, &alphaRho);

        // t-channel amplitude contributions are just sums of the above Regge poles
        inline std::complex<double> isovector(double s, double t)
        {
            return _rho1(s,t) + _rho2(s,t);
        };
        inline std::complex<double> isoscalar(double s, double t)
        {
            return _pom(s,t) + _f(s,t);
        };
    };
};

#endif

    