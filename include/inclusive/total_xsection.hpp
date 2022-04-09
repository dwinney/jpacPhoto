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
        total_xsection(){};
        
        // Only thing that is needed is a way to evaluate the cross-section
        virtual double eval(double s) = 0;
    };

    // ---------------------------------------------------------------------------
    // Coded up PDG parameterization for the asymptotic cross-sections

    class total_xsection_PDG : public total_xsection
    {
        public:

        // Constructor with masses and a filename to look for data
        total_xsection_PDG(double m1, double m2, std::array<double, 5> pdgparams, std::string datfile = "")
        : _mBeam(m1), _mTarget(m2), _threshold((m1+m2)*(m1+m2)), interp(0, ROOT::Math::Interpolation::kLINEAR)
        {
            _iso   = pdgparams[0];
            _delta = pdgparams[1];
            _R1    = pdgparams[2];
            _R2    = pdgparams[3];
            _P     = pdgparams[4];  
            
            if (datfile != "")
            {
                import_data(datfile);
            };
        };

        double eval(double s);

        private:

        double _mBeam, _mTarget;
        double _threshold;
        double _cutoff = 10.;

        // Lab beam momentum
        inline double pLab(double s) const
        { 
            double Elab = (s - _mBeam*_mBeam - _mTarget*_mTarget) / (2.*_mTarget); 
            return sqrt(Elab*Elab - _mBeam*_mBeam);
        };

        // If theres data available interpolate between this first
        ROOT::Math::Interpolator interp;
        std::vector<double> _plab, _sigma;
        void import_data(std::string datfile);

        // Else use the fits from the PDG 
        inline double PDG_parameterization(double s)
        {
            // These params are universal for all collisions
            double M = 2.1206, H = 0.2720, eta1 = 0.4473, eta2 = 0.5486;

            double sab;
            if (_mBeam < 1.E-4) sab = pow(M_RHO + _mTarget + M, 2);
            else sab = pow(_mBeam + _mTarget + M, 2);

            return _delta*(H*pow(log(s/sab), 2.) + _P) + _R1*pow(s/sab, -eta1) - _iso*_R2*pow(s/sab, -eta2);
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