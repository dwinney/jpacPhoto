// Discontinuity across unitarity cut for vector production via a one-loop box diagram
// Used as a container class in integration processes
//
// box_discontinuity is abstract to allow different parameterizations
//
// brute_force_discontinuity will try to calculate the intermediate phase-space integration numerically
// We consider the OVERALL process: B T -> X R (Beam + Target -> Meson + Recoil)
// through an intermediate state:   Xp Rp 
// 
// The "left"  amplitude corresponds to B T -> Xp Rp 
// the "right" amplitude corresponds to X R -> Xp Rp 
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef BOX_DISC
#define BOX_DISC

#include "constants.hpp"
#include "amplitude.hpp"
#include "reaction_kinematics.hpp"

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "Math/IntegratorMultiDim.h"

namespace jpacPhoto
{
    // Abstract class for a generic discontinuity
    class box_discontinuity
    {
        public: 

        box_discontinuity(reaction_kinematics * xkinem)
        : _kinematics(xkinem)
        {};

        virtual ~box_discontinuity() = default;

        // Alternatively evaluate the full dispersion relation
        virtual std::complex<double> helicity_amplitude(std::array<int,4> helicities, double s, double t) = 0;

        // Number of parameters needed
        // By default we require no additional parameters
        virtual int  get_nParams() = 0;
        virtual void set_params(std::vector<double> params) { _params.clear(); _params = params; };
        
        inline void set_intermediate_threshold(double wth){ _intermediateThreshold = wth*wth; };

        protected:

        // All kinematic quantities of the overall scattering process
        reaction_kinematics * _kinematics;
        double _intermediateThreshold;
        std::vector<double> _params;
    };

    // -------------------------------------------------------------------
    // Attempt to do calculate the phase-space integral analytically

    class brute_force_discontinuity : public box_discontinuity
    {
        public:

        // Constructor requires a kinematics object for the OVERALL process
        // Then two amplitudes for the left and right subprocesses
        brute_force_discontinuity(reaction_kinematics * xkinem, amplitude * left, amplitude * right)
        : box_discontinuity(xkinem), _left_amp(left), _right_amp(right)
        {
            compatibility_check(xkinem, left, right);
        };

        // Only need one parameter 
        int get_nParams(){ return 1; };

        // Try doing the double integral over the intermediate phasepace
        // and sum over intermediate helicities numerically
        double eval(double s);

        // Evaluate the dispersion relation and spit out the helicity amplitude
        std::complex<double> helicity_amplitude(std::array<int,4> helicities, double s, double t);

        // -------------------------------------------------------------------
        protected:

        // External quantities saved for easy access
        double             _external_theta;
        std::array<int, 4> _external_helicities;

        // The sub-process amplitudes
        amplitude * _left_amp;
        amplitude * _right_amp;

        // Save an array of the helicities for easier indexing sums later
        std::vector<std::array<int,4>> _intermediate_helicities;

        // Two-body phase-space
        inline double phase_space(double s)
        {
            if (s < _left_amp->_kinematics->sth()) return 0.;
            return _left_amp->_kinematics->final_momentum(s) / (8. * PI * sqrt(s));
        };

        // At initialization, check that the left and right amplitudes are kinematically compatible
        void compatibility_check(reaction_kinematics * kinem, amplitude * left, amplitude * right);
        void qn_error(amplitude * left, amplitude * right);
        void mass_error(reaction_kinematics * kinem, amplitude * left, amplitude * right);
        bool _match_error = false;
    };

    // ---------------------------------------------------------------------------
    // Testing discontinuity which return a helicity-conserving delta function
    
    class flat_discontinuity : public box_discontinuity
    {
        public:
        // Do not require any amplitudes
        flat_discontinuity(reaction_kinematics * xkinem)
        : box_discontinuity(xkinem)
        {};

        inline std::complex<double> helicity_amplitude(std::array<int,4> helicities, double s, double t)
        {
            // if below threshold return 0
            if (s < _kinematics->sth()) return 0.;
            
            // make sure the external helicites get passed correctly
            int lam_gam = helicities[0];
            int lam_tar = helicities[1];
            int lam_vec = helicities[2];
            int lam_rec = helicities[3];

            if (lam_gam != lam_vec || lam_rec != lam_tar) return 0.;
            return 1.;
        };
    };
};

#endif