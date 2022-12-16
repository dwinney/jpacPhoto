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
#include "cgamma.hpp"

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
        total_xsection(double mb, double mt)
        : _mBeam(mb), _mTarget(mt), 
          _sth((mb + mt)*(mb + mt))
        {};
        
        virtual ~total_xsection() = default;

        // Only thing that is needed is a way to evaluate the cross-section
        // A second argument can be passed in case the method has a 
        // way to incorporate pion virtuality effects
        virtual double eval(double s, double q2) = 0;

        void set_debug(int i){ _debug = i; };

        protected:
        
        // Lab beam momentum
        inline double pLab(double s) const
        { 
            double Elab = (s - _mBeam*_mBeam - _mTarget*_mTarget) / (2.*_mTarget); 
            return sqrt(Elab*Elab - _mBeam*_mBeam);
        };

        int _debug = 0;
        double _mBeam, _mTarget;
        double _sth;
    };

    // ---------------------------------------------------------------------------
    // "Zero" instance of the total_xsection to use as a default object
    class zero_xsection : public total_xsection
    {
        public: 
        zero_xsection()
        : total_xsection(0.,0.)
        {};

        inline double eval(double s, double q2)   { return 0.; };
    };

};

#endif