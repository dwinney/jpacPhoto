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
        : amplitude(xkinem, name), _mEx2(mass*mass), _reggeized(false)
        {
            set_nParams(2);
            check_JP(xkinem->_jp);
        };

        // constructors for regge exchange
        pseudoscalar_exchange(reaction_kinematics * xkinem, linear_trajectory * traj, std::string name = "pseudoscalar_exchange")
        : amplitude(xkinem, name), _alpha(traj), _reggeized(true)
        {
            set_nParams(2);
            check_JP(xkinem->_jp);

            if ((xkinem->_jp != AXIAL_VECTOR) && (xkinem->_jp != VECTOR))
            {
                std::cout << "Error! Only reggeized axial-vector and vector production implemented so far!" << std::endl;
                exit(0);
            }
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

        // only axial-vector, vector, and pseudo-scalar available
        inline std::vector<std::array<int,2>> allowedJP()
        {
            return { AXIAL_VECTOR, VECTOR, PSEUDO_SCALAR };
        };

        private:

        // Whether to use fixed-spin propagator (false) or regge (true)
        bool _reggeized = false;


        // Mass of the exchanged pseudo-scalar (if REGGE = false)
        // ignored otherwise
        double _mEx2;

        // Regge trajectory for the pion (if REGGE = true)
        // ignored otherwise
        linear_trajectory * _alpha;

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
