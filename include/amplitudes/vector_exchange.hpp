// Axial-vector meson photoproduction proceeding through a vector meson exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _AXIAL_
#define _AXIAL_

#include "amplitude.hpp"
#include "regge_trajectory.hpp"

// ---------------------------------------------------------------------------
// vector_exchange class describes the amplitude for a fixed-spin-1 exchange
// in the t-channel. Derived in terms of simple feynman rules at tree level
//
// Initialization required a reaction_kinematics object, the mass of the exchange,
// and an optional string to identify the amplitude with.
//
//  Evaluation requires three couplings photon coupling, gGamma, and vector/tensor
// nucleon couplings, gV and gT respectively.
//
// Set couplings with amp.set_params({gGamma, gV, gT});
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
    class vector_exchange : public amplitude
    {
        public:

        // Constructor for fixed spin
        vector_exchange(reaction_kinematics * xkinem, double mass, std::string id = "vector_exchange")
        : amplitude(xkinem, id), _mEx2(mass*mass), _ifReggeized(false)
        {
            set_nParams(3);
            check_JP(xkinem->_jp);

            // Analytical residues only available for axial-vector production
            if (xkinem->_jp == SCALAR) _useCovariant = true;
        };

        // Constructor for the reggized)
        vector_exchange(reaction_kinematics * xkinem, linear_trajectory * traj, std::string id = "vector_exchange")
        : amplitude(xkinem, id), _alpha(traj), _ifReggeized(true)
        {
            set_nParams(3);
            check_JP(xkinem->_jp);

            if (xkinem->_jp != AXIAL_VECTOR)
            {
                std::cout << "Error! Reggeized vector_exchange only implemented for axial vector production so far.\n";
                exit(0);
            }
        };

        // Setting utility
        inline void set_params(std::vector<double> params)
        {
            check_nParams(params); // make sure the right amout of params passed
            _gGam = params[0];
            _gV = params[1];
            _gT = params[2];
        };

        // Whether or not to include an exponential form factor (default false)
        inline void set_formfactor(int FF, double bb = 0.)
        {
            _useFormFactor = FF;
            _cutoff = bb;
        }

        // Assemble the helicity amplitude by contracting the lorentz indices
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t);

        // axial vector and scalar kinematics allowed
        inline std::vector<std::array<int,2>> allowedJP()
        {
            return {AXIAL_VECTOR, VECTOR, SCALAR, PSEUDO_SCALAR};
        };

        private:

        // if using reggeized propagator
        bool _ifReggeized;
        // or the regge trajectory of the exchange
        linear_trajectory * _alpha;
        double _zt;

        // Whether using analytic or covariant expression
        bool _useCovariant = false;

        // Form factor parameters
        int _useFormFactor = 0;
        double _cutoff = 0.;
        double form_factor();

        // Couplings to the axial-vector/photon and vector/tensor couplings to nucleon
        double _gGam = 0., _gpGam = 0., _gV = 0., _gT = 0.;

        // ---------------------------------------------------------------------------
        // Covariant evaluation

        // Mass of the exchange
        double _mEx2 = 0.;

        // Full covariant amplitude
        std::complex<double> covariant_amplitude(std::array<int, 4> helicities);

        // Photon - Axial Vector - Vector vertex
        std::complex<double> top_vertex(int mu, int lam_gam, int lam_vec);

        // Nucleon - Nucleon - Vector vertex
        std::complex<double> bottom_vertex(int nu, int lam_targ, int lam_rec);

        // Vector propogator
        std::complex<double> vector_propagator(int mu, int nu);

        // ---------------------------------------------------------------------------
        // Analytic evaluation

        // Photon - Axial - Vector
        std::complex<double> top_residue(int lam_gam, int lam_vec);

        // Nucleon - Nucleon - Vector
        std::complex<double> bottom_residue(int lam_targ, int lam_rec);

        // Reggeon propagator
        std::complex<double> regge_propagator(int j, int lam, int lamp);

        // Half angle factors
        std::complex<double> half_angle_factor(int lam, int lamp);

        // Angular momentum barrier factor
        std::complex<double> barrier_factor(int j, int M);
    };
};

#endif
