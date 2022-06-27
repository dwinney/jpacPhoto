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
        // A second argument can be passed in case the method has a 
        // way to incorporate pion virtuality effects
        double eval(double s, double q2 = M2_PION)
        {
            double result = 0.;

            if (s < _sth + 10.*EPS)
            {
                return 0.;
            }
            else if (s <= _cutoff)
            {
                return resonances(s, q2);
            }
            else
            {
                return regge(s);  
            };
        };

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

        // Pion virtualilty is assumed to not matter in the Regge region
        virtual double resonances(double s, double q2) = 0;
        virtual double regge(double s) = 0;
    };

    // ---------------------------------------------------------------------------
    // "Zero" instance of the total_xsection to use as a default object
    class zero_xsection : public total_xsection
    {
        public: 
        zero_xsection()
        : total_xsection(0.,0.,0.)
        {};

        inline double resonances(double s, double q2)   { return 0.; };
        inline double regge(double s)                   { return 0.; };
    };

};

#endif