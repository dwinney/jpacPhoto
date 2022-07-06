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
        : total_xsection(M_PION, M_PROTON), _iso(iso), _resonances(resonances),
          _interp(0, ROOT::Math::Interpolation::kCSPLINE)
        {
            if (abs(iso) != 1)
            {
                std::cout << "Error! JPAC_parameterization argument must be +1 or -1! Results may vary..." << std::endl;
            };

            if (resonances) initialize_PWAs();
            this->set_debug(3);
        };

        ~JPAC_parameterization()
        {
            for (int i = 0; i <= _Lmax; i++)
            {
                delete _pw_1p[i]; delete _pw_1m[i];
                delete _pw_3p[i]; delete _pw_3m[i];
            };
        };

        double eval(double s, double q2)
        {
            double _cutoff = 4.;
            double x1 = _cutoff;
            double x2 = _cutoff + 0.12;

            if ( sqrt(s) < sqrt(_sth) + 1000.*EPS ) return 0.;
            else if ( _resonances && s <= x1 ) 
            {
                return sigma_resonance(s, q2);
            }
            else if ( _resonances && (s > x1) && (s < x2) )
            {
                double fx1 = sigma_resonance(x1, M2_PION);
                double fx2 = sigma_regge(x2);
                double y  = (s - x1) / (x2 - x1);
                return fx1 * (1. - y) + fx2 * y;
            }
            else 
            {
                return sigma_regge(s);
            }

            return 0.;
        };

        protected:
        // Parameter selecting pi + or pi - scattering
        int  _iso = 1;  // plus or minus 1
        bool _resonances = false;

        // ----------------------------------------------------------------------
        // Low-energy pieces

        // Calculate the total cross-section be evaluating the forward imaginary part from individual partial waves
        class SAID_PWA
        {
            public:

            SAID_PWA(double L, double I, double J)
            : _L(L), _I(I), _J(J)
            {
                if (J >= 0) import_data(); 
            }
            
            // Lab beam momentum
            inline double pLab(double s) const
            { 
                double Elab = (s - M2_PION - M2_PROTON) / (2.*M_PROTON); 
                return sqrt(Elab*Elab - M2_PION);
            };
            
            double imaginary_part(double s)
            {
                if ( _J < 0) return 0.;
                return _interpImag.Eval( pLab(s) );
            };

            double real_part(double s)
            {
                if ( _J < 0) return 0.;
                return _interpReal.Eval( pLab(s) );
            };

            std::complex<double> eval(double s)
            {
                if ( _J < 0) return 0.;
                return XR * real_part(s) + XI * imaginary_part(s);
            };

            int orbital_spin(){ return _L; };
            
            private:

            // Read in the data file and interpolate the imaginary part of the amplitude
            void import_data();
            ROOT::Math::Interpolator _interpReal, _interpImag;
            std::vector<double> _plab, _realPart, _imagPart;

            int _L, _I, _J; 
            std::string get_filename();
        };

        // If theres data available in the low enegergy, 
        // these methods read data files and store them
        ROOT::Math::Interpolator _interp;
        std::vector<double> _plab, _sigma;
        void import_data(std::string datfile = "SAID.dat");

        // For given isospin set up the individual partial waves and store them in a vector
        int _Lmax = 7.;
        void initialize_PWAs();
        std::vector<SAID_PWA*> _pw_1p, _pw_1m, _pw_3p, _pw_3m;

        // First derivative of Legendre function evaluated at z = 1
        double coeffP(int l){ return double(l*(l+1))/2.; };

        void update_amplitudes(double s);
        std::vector<double> _CpL, _CmL;
        double sigma_resonance(double s, double q2)
        {
            double result = 0.;
            double bfpi  = Kallen(s, q2, M2_PROTON) / Kallen(s, M2_PION, M2_PROTON);

            update_amplitudes(s);
            for (int L = 0; L <= _Lmax; L++)
            {
                int LL;
                (L < _debug) ? (LL = L) : (LL = _debug);
                
                if (L >= 3 && s < 1.2) continue;

                double sigmaL = 0.389352 * ( _CpL[L] - double(_iso) * _CmL[L] ) / pLab(s);
                if (sigmaL < 0) continue;
                result += pow( bfpi, double(LL) ) * sigmaL ;
            };  

            return result;
        };

        // ----------------------------------------------------------------------
        // High-energy pieces

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
        regge_pole _rho1 = regge_pole( {-5.01*1.573, 10.1}, &alphaRho);
        regge_pole _rho2 = regge_pole( { 5.01*0.573,  0.0}, &alphaRho);

        // t-channel amplitude contributions are just sums of the above Regge poles
        inline double isovector(double s, double t)
        {
            return std::imag(_rho1(s,t) + _rho2(s,t));
        };
        inline double isoscalar(double s, double t)
        {
            return std::imag(_pom(s,t) + _f(s,t));
        };
        inline double sigma_regge(double s) 
        {
            double Ap = isoscalar(s, 0.); double Am = isovector(s, 0.);
            return 0.389352 * ( Ap - double(_iso) * Am ) / pLab(s);
        };
    };
};

#endif

    