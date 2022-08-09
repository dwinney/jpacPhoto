// Charged axial-vector meson photoproduction proceeding through a pseudoscalar (pion) exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// References:
// [1] arXiv:1503.02125 [hep-ph]
// ---------------------------------------------------------------------------

#ifndef _PSCALAR_
#define _PSCALAR_

#include "amplitude.hpp"
#include "regge_trajectory.hpp"

// ---------------------------------------------------------------------------
// pseudoscalar_exchange class describes the amplitude for a spin-0 exchange
// in the t-channel. Derived in terms of simple feynman rules at tree level.
//
// Initialization required a reaction_kinematics object.
// Then either: the mass (in GeV) of the exchange (for fixed-spin exchange),
//          or: a pointer to a linear_trajectory object (for Reggeize exchange).
// and an optional string to identify the amplitude with.
//
// Evaluation requires two couplings:
// photon coupling, gGamma, and nucleon coupling, gNN respectively.
//
// Set couplings with amp.set_params({gGamma, gNN});
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
    class pseudoscalar_exchange : public amplitude
    {
        public:
        // constructor for fixed meson exchange
        pseudoscalar_exchange(reaction_kinematics * xkinem, double mass, std::string name = "pseudoscalar_exchange")
        : amplitude(xkinem, "pseudoscalar_exchange", name),
          _mEx2(mass*mass)
        {
            set_nParams(2);
            check_JP(xkinem);
        };

        // constructors for regge exchange
        pseudoscalar_exchange(reaction_kinematics * xkinem, linear_trajectory * traj, std::string name = "pseudoscalar_exchange")
        : amplitude(xkinem, "pseudoscalar_exchange", name), 
          _alpha(traj)
        {
            set_nParams(2);
            check_JP(xkinem);
            _reggeized = true;
        };

        // Setting utility
        void set_params(std::vector<double> params)
        {
            check_nParams(params);
            _gT = params[0];
            _gB = params[1];
        };

        // Whether or not to include an exponential form factor (default false)
        void set_formfactor(int FF, double bb = 0.)
        {
            _useFormFactor = FF;
            _cutoff = bb;
        }

        // Assemble the helicity amplitude by contracting the spinor indices
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double xs, double xt);

        // Covariant quantities define lambda in S-channel
        // Analytic ones in the T-channel
        helicity_channel helicity_CM_frame()
        {
            if (_useCovariant)
            {
                return helicity_channel::S;
            }
            else
            {
                return helicity_channel::T;
            }
        };

        // only axial-vector, vector, and pseudo-scalar available
        inline std::vector<std::array<int,2>> allowed_meson_JP()
        {
            if (!_reggeized)
            {
                return { {1, +1}, {1, -1}, {0, -1} };
            }
            else
            {
                return { {1, +1}, {1, -1} };
            }
        };
        inline std::vector<std::array<int,2>> allowed_baryon_JP()
        {
            if (!_reggeized) return { {1, +1}, {3, +1} };
            else return { {1, +1} };
        };


        // Accessor functions for private memebers 
        inline bool if_reggeized(){ return _reggeized; };
        inline double get_mEx2(){ return _mEx2; };
        inline linear_trajectory * get_trajectory(){ return _alpha; };
        inline double get_cutoff(){ return _cutoff; };
        inline double get_coupling(){ return _gT; };

        // return the coupling function for the top vertex
        // For use with inclusive
        inline double top_coupling(double t)
        { 
            update({0,0,0,0}, 0, t);
            return std::real(top_residue()); 
        };

        // --------------------------------------------------------------------

        private:

        // Whether to use fixed-spin propagator (false) or regge (true)
        std::complex<double> _qt, _pt; // Momentum in t-channel

        // Mass of the exchanged pseudo-scalar (if REGGE = false)
        // ignored otherwise
        double _mEx2;

        // Regge trajectory for the pion (if REGGE = true)
        // ignored otherwise
        linear_trajectory * _alpha = NULL;

        // Coupling constants
        double _gT = 0.; // Gamma - Axial - Pseudoscalar coupling 
        double _gB = 0.;    // Pseudoscalar - Nucleon coupling

        // Exchange form-factor
        int _useFormFactor = 0;   // Whether to include the exponential form factor
        double _cutoff = 0.;      // "t-slope" parameter in the FF
        double form_factor();

        // Analytic residues
        std::complex<double> top_residue();
        std::complex<double> bottom_residue();

        // Covariant vertices
        std::complex<double> top_vertex();
        std::complex<double> axialvector_coupling();
        std::complex<double> vector_coupling();
        std::complex<double> pseudoscalar_coupling();

        std::complex<double> bottom_vertex();
        std::complex<double> halfplus_coupling();
        std::complex<double> threehalvesplus_coupling();

        // Simple pole propagator
        std::complex<double> scalar_propagator();
    };
};

#endif
