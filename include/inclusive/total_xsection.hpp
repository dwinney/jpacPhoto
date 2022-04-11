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
    // Coded up PDG parameterization for the asymptotic cross-sections

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

        double eval(double s);

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

    // // Pi+ Proton 
    // total_xsection_PDG total_xsection_pipp(M_PION, M_PROTON, {+1., 1., 9.56, 1.767, 18.75}, "rpp2020-pipp_total.dat");

    // // Pi- Proton
    // total_xsection_PDG total_xsection_pimp(M_PION, M_PROTON, {-1., 1., 9.56, 1.767, 18.75}, "rpp2020-pimp_total.dat");

    // // Gamma Proton
    // total_xsection_PDG total_xsection_gamp(0.,     M_PROTON, { 0., 3.065E-3, 0.0139, 0., 34.41}, "rpp2020-gammap_total.dat");
};

#endif