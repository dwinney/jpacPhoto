// Simplified version of the 2->2 box diagram with contact S-wave vertices
// i.e. S-wave angular behavior and energy depedence given by scalar bubble integral

// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef BUBBLE
#define BUBBLE

#include "amplitude.hpp"

#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace jpacPhoto
{
    class bubble_amplitude : public amplitude
    {   
        public:

        // Constructor should specify the masses of the intermediate hadrons in the bubble
        bubble_amplitude(reaction_kinematics * xkinem, std::array<double,2> masses, std::string id = "bubble_amplitude")
        : amplitude(xkinem, "bubble_amplitude", id), _m1(masses[0]), _m2(masses[1])
        {
            set_nParams(2);
            check_JP(xkinem);
        };

        // Constructor should specify the masses of the intermediate hadrons in the bubble
        bubble_amplitude(reaction_kinematics * xkinem, std::array< std::array<double,2>, 2> masses, std::string id = "bubble_amplitude")
        : amplitude(xkinem, "bubble_amplitude", id), _m1(masses[0][0]), _m2(masses[0][1]), _m3(masses[1][0]), _m4(masses[1][1]),
          _twochannel(true)
        {
            set_nParams(4);
            check_JP(xkinem);
        };


        // Setting utility for free parameters
        void set_params(std::vector<double> params)
        {
            check_nParams(params);
            _norm   = params[0];
            _eps    = params[1];
            _n      = params[2];

            if (_twochannel) _r = params[3];
        };

        // Scalar bubble integral function
        std::complex<double> G0(double s, double m1, double m2);
        std::complex<double>  G(double s, double m1, double m2);
      
        // Assemble the helicity amplitude by contracting the lorentz indices for contact interaction
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t);

        // only vector and proton quantum numbers available currently
        inline std::vector<std::array<int,2>> allowed_meson_JP()
        {
            return { {1, -1} };
        };
        inline std::vector<std::array<int,2>> allowed_baryon_JP()
        {
            return { {1,  1} };
        };
        
        // Box always works with s-channel helicity projections
        inline helicity_channel helicity_CM_frame(){ return S; };

        private:

        // 
        bool _twochannel = false;

        // Short-hand function for (ubar . gamma_mu . u) 
        std::complex<double> baryon_current(int mu);

        // Masses of the internal propagating hadrons
        double _m1, _m2;
        // if we have a second open channel 
        double _m3, _m4;

        // cutoff and subtraction (two free parameters)
        double _norm;
        // If considering two channels simultaneously we have an additional parameter
        double _r, _eps, _n, _s0 = 16.;
    };
};

#endif 