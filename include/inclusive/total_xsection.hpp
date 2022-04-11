// Phenomenological expressions for the total cross-sections.
// We use a generic class callable by double total_xsection(double) to select different
// parameterizations or reactions
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef SIGMA_TOT
#define SIGMA_TOT

#include "constants.hpp"
#include "regge_trajectory.hpp"
#include "misc_math.hpp"

#include <vector>
#include <array>
#include <fstream>
#include <sstream>

#include <Math/Interpolator.h>

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Generic template class

    class total_xsection
    {
        public:

        // Default constructor 
        total_xsection(double mb, double mt, double cutoff = 0.)
        : _mBeam(mb), _mTarget(mt), 
          _sth((mb + mt)*(mb + mt)),
          _cutoff(cutoff)
        {};
        
        // Only thing that is needed is a way to evaluate the cross-section
        double eval(double s);

        protected:
        
        // Lab beam momentum
        inline double pLab(double s) const
        { 
            double Elab = (s - _mBeam*_mBeam - _mTarget*_mTarget) / (2.*_mTarget); 
            return sqrt(Elab*Elab - _mBeam*_mBeam);
        };

        double _mBeam, _mTarget;
        double _sth;
        double _cutoff;

        // To evaluate the cross-section we use two regions
        // some way to handle the resonances at low s < _cutoff
        // and the regge dominated behavior at high s > _cutoff
        virtual double resonances(double s) = 0;
        virtual double regge(double s) = 0;
    };

    // ---------------------------------------------------------------------------
    // Total cross-sections based off of the 2016 PDG parameterizations

    class PDG_parameterization : public total_xsection
    {
        public:

        // Constructor with masses and a filename to look for data
        PDG_parameterization(double m1, double m2, std::array<double, 5> pdgparams, std::string datfile = "")
        : total_xsection(m1, m2, 0.), _interp(0, ROOT::Math::Interpolation::kLINEAR)
        {
            _iso   = pdgparams[0];  // isospin sign for pi+ or pi- dependence
            _delta = pdgparams[1];  // VMD/Quark counting paramter for photon 
            _R1    = pdgparams[2];  // Regge terms
            _R2    = pdgparams[3];
            _P     = pdgparams[4];  // Constant pomeron term
            
            if (datfile != "")
            {
                _cutoff = 6.;
                import_data(datfile);
            };
        };

        protected:

        // Above the cutoff, we use the PDG parameterization of Regge behavior
        inline double regge(double s)
        {
            // These params are universal for all collisions
            double M = 2.1206, H = 0.2720, eta1 = 0.4473, eta2 = 0.5486;

            double sab;
            if (_mBeam < 1.E-4) sab = pow(M_RHO + _mTarget + M, 2);
            else sab = pow(_mBeam + _mTarget + M, 2);

            return _delta*(H*pow(log(s/sab), 2.) + _P) + _R1*pow(s/sab, -eta1) - _iso*_R2*pow(s/sab, -eta2);
        };

        // If theres data available in the low enegergy, 
        // these methods read data files and store them
        ROOT::Math::Interpolator _interp;
        std::vector<double> _plab, _sigma;
        void import_data(std::string datfile);

        // Below we do a simple linear interpolation of data in resonance region
        inline double resonances(double s)
        {
            return _interp.Eval( pLab(s) );
        };

        int _iso = 0;
        double _delta = 1.;
        double _R1, _R2, _P;
    };

    
    // ---------------------------------------------------------------------------
    // Total cross-sections based off of the 2015 JPAC publication
    // V. Mathieu et al. [Phys. Rev. D 92, 074004]
    class JPAC_parameterization : public total_xsection
    {
        public:

        // Currently only avialable for pion-nucleon 
        JPAC_parameterization(int iso)
        : total_xsection(M_PION, M_PROTON, 6.5), _iso(iso),
          _interp(0, ROOT::Math::Interpolation::kCSPLINE)
        {
            if (abs(iso) != 1)
            {
                std::cout << "Error! JPAC_parameterization argument must be +1 or -1! Results may vary..." << std::endl;
            };

            import_data();
        };

        protected:
        // ----------------------------------------------------------------------
        // Low-energy pieces

        // Below the cutoff we do a simple interpolation of the SAID amplitude
        inline double resonances(double s)
        {
            return _interp.Eval( pLab(s) );
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