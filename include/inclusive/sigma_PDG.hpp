// Phenomenological expressions for the total cross-sections from the PDG
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef SIGMA_TOT_PDG
#define SIGMA_TOT_PDG

#include "total_xsection.hpp"

namespace jpacPhoto
{
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

        // If theres data available in the low energy, 
        // these methods read data files and store them
        ROOT::Math::Interpolator _interp;
        std::vector<double> _plab, _sigma;
        void import_data(std::string datfile);

        // Below we do a simple linear interpolation of data in resonance region
        inline double resonances(double s, double q2)
        {
            return _interp.Eval( pLab(s) );
        };

        int _iso = 0;
        double _delta = 1.;
        double _R1, _R2, _P;
    };
};

#endif

    