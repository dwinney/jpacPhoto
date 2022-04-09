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
        : amplitude(xkinem, "pseudoscalar_exchange", name), _mEx2(mass*mass), _reggeized(false)
        {
            set_nParams(2);
            check_JP(xkinem->_jp);
        };

        // constructors for regge exchange
        pseudoscalar_exchange(reaction_kinematics * xkinem, linear_trajectory * traj, std::string name = "pseudoscalar_exchange")
        : amplitude(xkinem, "pseudoscalar_exchange", name), _alpha(traj), _reggeized(true)
        {
            set_nParams(2);
            check_JP(xkinem->_jp, true);
        };

        // Setting utility
        void set_params(std::vector<double> params)
        {
            check_nParams(params);
            _gGamma = params[0];
            _gNN = params[1];
        };

        // Whether or not to include an exponential form factor (default false)
        void set_formfactor(int FF, double bb = 0.)
        {
            _useFormFactor = FF;
            _cutoff = bb;
        }

        // Assemble the helicity amplitude by contracting the spinor indices
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double xs, double xt);

        inline int parity_phase(std::array<int, 4> helicities)
        {
            if (_useCovariant || _debug >= 1)
            {
                return _kinematics->parity_phase(helicities, HELICITY_CHANNEL::S);
            }
            else
            {
                return _kinematics->parity_phase(helicities, HELICITY_CHANNEL::T);
            }

            return 0.;
        };

        // only axial-vector, vector, and pseudo-scalar available
        inline std::vector<std::array<int,2>> allowedJP()
        {
            return { AXIAL_VECTOR, VECTOR, PSEUDO_SCALAR };
        };

        // only axial-vector, vector, and pseudo-scalar available
        inline std::vector<std::array<int,2>> allowedJP_Regge()
        {
            return { AXIAL_VECTOR, VECTOR };
        };

        // Accessor functions for private memebers 
        inline bool if_reggeized(){ return _reggeized; };
        inline double get_mEx2(){ return _mEx2; };
        inline linear_trajectory * get_trajectory(){ return _alpha; };
        inline double get_cutoff(){ return _cutoff; };
        inline double get_coupling(){ return _gGamma; };

        // return the coupling function for the top vertex
        // For use with inclusive
        inline double top_coupling(double t)
        { 
            _t = t;
            return std::real(top_residue(0,0)); 
        };

        private:

        // Whether to use fixed-spin propagator (false) or regge (true)
        bool _reggeized = false;

        // Mass of the exchanged pseudo-scalar (if REGGE = false)
        // ignored otherwise
        double _mEx2;

        // Regge trajectory for the pion (if REGGE = true)
        // ignored otherwise
        linear_trajectory * _alpha = NULL;

        // Coupling constants
        double _gGamma = 0.; // Gamma - Axial - Pseudoscalar coupling 
        double _gNN = 0.;    // Pseudoscalar - Nucleon coupling

        int _useFormFactor = 0;   // Whether to include the exponential form factor
        double _cutoff = 0.;      // "t-slope" parameter in the FF
        double form_factor();

        // Whether to switch to using the feynman rules
        bool _useCovariant = false; 

        // Photon - pseudoscalar - Axial vertex
        std::complex<double> top_residue(int lam_gam, int lam_vec);
        std::complex<double> top_vertex( int lam_gam, int lam_vec);

        // Pseudoscalar - Nucleon - Nucleon vertex
        std::complex<double> bottom_residue(int lam_targ, int lam_rec);
        std::complex<double> bottom_vertex( int lam_targ, int lam_rec);

        // Simple pole propagator
        std::complex<double> scalar_propagator();
    };
};

#endif
