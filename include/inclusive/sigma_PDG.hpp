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
        PDG_parameterization(double m1, double m2, std::array<double, 5> pdgparams)
        : total_xsection(m1, m2)
        {
            _iso   = pdgparams[0];  // isospin sign for pi+ or pi- dependence
            _delta = pdgparams[1];  // VMD/Quark counting paramter for photon 
            _R1    = pdgparams[2];  // Regge terms
            _R2    = pdgparams[3];
            _P     = pdgparams[4];  // Constant pomeron term
        };

        // Above the cutoff, we use the PDG parameterization of Regge behavior
        inline double eval(double s, double q2)
        {
            // These params are universal for all collisions
            double M = 2.1206, H = 0.2720, eta1 = 0.4473, eta2 = 0.5486;

            double sab;
            if (_mBeam < 1.E-4) sab = pow(M_RHO + _mTarget + M, 2);
            else sab = pow(_mBeam + _mTarget + M, 2);

            return _delta*(H*pow(log(s/sab), 2.) + _P) + _R1*pow(s/sab, -eta1) - _iso*_R2*pow(s/sab, -eta2);
        };

        protected:

        int _iso = 0;
        double _delta = 1.;
        double _R1, _R2, _P;
    };
};

#endif

    