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
        : amplitude(xkinem, "vector_exchange", id),
          _mEx2(mass*mass)
        {
            set_nParams(3);
            check_JP(xkinem);

            // Analytical residues not available for scalar quantum numbers yet
            std::array<int,2> jp = xkinem->get_meson_JP();
            if ( jp[0] == 0 && jp[0] == +1 ) _useCovariant = true;

            // If massive beam use covariant to avoid 1/t instabilities
            if (!xkinem->is_photon())        _useCovariant = true;
        };

        // Constructor for the reggized)
        vector_exchange(reaction_kinematics * xkinem, linear_trajectory * traj, std::string id = "vector_exchange")
        : amplitude(xkinem, "vector_exchange", id),
          _alpha(traj)
        {
            set_nParams(3);
            check_JP(xkinem);
            this->_reggeized = true;
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

        // axial vector and scalar kinematics allowed
        inline std::vector<std::array<int,2>> allowed_meson_JP()
        {
            if (!_reggeized)
            {
                return { {1, +1}, {1, -1}, {0, +1}, {0, -1} };
            }
            else
            {
                return { {1, +1} };
            }
        };
        inline std::vector<std::array<int,2>> allowed_baryon_JP()
        {
            return { {1, +1} };
        }; 

        private:

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
        std::complex<double> covariant_amplitude();

        // Beam -- Vector -- Produced meson vertex
        std::complex<double> top_vertex(int mu);

        // Depends on quantum numbers and carries one lorentz index
        std::complex<double> axialvector_coupling(int mu);
        std::complex<double> vector_coupling(int mu);
        std::complex<double> pseudoscalar_coupling(int mu);
        std::complex<double> scalar_coupling(int mu);

        // Target -- Vector -- Produced baryon vertex
        std::complex<double> bottom_vertex(int nu);

        // Vector propogator, with two Lorentz indices
        std::complex<double> vector_propagator(int mu, int nu);

        // ---------------------------------------------------------------------------
        // Analytic evaluation
        
        // if using reggeized propagator
        linear_trajectory * _alpha;  // regge trajectory of the exchange

        // Net helicities
        int _lam, _lamp, _M;

        double _zt; // Scattering angle in the t-channel
        std::complex<double> _qt; // momentum in t-channel
        
        // Full analytic amplitude
        std::complex<double> analytic_amplitude();

        // Photon - Axial - Vector
        std::complex<double> top_residue();
        std::complex<double> axialvector_residue();
        std::complex<double> vector_residue();
        std::complex<double> pseudoscalar_residue();

        // Nucleon - Nucleon - Vector
        std::complex<double> bottom_residue();

        // Reggeon propagator
        std::complex<double> regge_propagator();

        // Half angle factors
        std::complex<double> half_angle_factor();

        // Angular momentum barrier factor
        std::complex<double> barrier_factor();
    };
};

#endif
