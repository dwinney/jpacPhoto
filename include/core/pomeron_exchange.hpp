// Vector meson photoproduction dynamics proceeding through a pomeron exchange
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
// REFERENCES:
// [1] 1907.09393
// [2] 1606.08912
// [3] 1904.11706
// ---------------------------------------------------------------------------

#ifndef _POMERON_
#define _POMERON_

#include "amplitude.hpp"
#include "regge_trajectory.hpp"

// ---------------------------------------------------------------------------
// The pomeron_exchange class describes the amplitude correspinding to
// vector meson photoproduction via a vector pomeron coupling.
//
// As a reggeon exchange model it is parameterized in terms three functions
// 1. top_vertex() coupling the vector meson to the incoming photon
// 2. bottom_vertex() coupling the two proton dirac spinors
// 3. regge_factor() the function describing the energy dependence of the amplitude
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
    class pomeron_exchange : public amplitude
    {
        public:

        // Constructor
        // need a pointer to kinematic object, pointer to trajectory.
        pomeron_exchange(reaction_kinematics * xkinem, regge_trajectory * alpha, int model = 0, std::string name = "pomeron_exchange")
        : amplitude(xkinem, "pomeron_exchange", name), _traj(alpha), _model(model)
        {
            set_nParams(2);
            check_JP(xkinem);
        };

        // Setting utility
        void set_params(std::vector<double> params)
        {
            check_nParams(params);
            _norm = params[0];
            _b0 = params[1];
        };

        // Assemble the helicity amplitude by contracting the lorentz indices
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t);

        inline int parity_phase(std::array<int, 4> helicities)
        {
            return _kinematics->parity_phase(helicities, HELICITY_CHANNEL::S);
        };

        // only vector kinematics allowed
        inline std::vector<std::array<int,2>> allowed_meson_JP()
        {
            return {{1, -1}};
        };
        inline std::vector<std::array<int,2>> allowed_baryon_JP()
        {
            return {{1,  1}};
        };

        private:
        
        // Which model to use. 
        int _model = 0; 
        // 0 - model in [1] Lesniak-Szcepaniak
        // 1 - model in [2] Helicity conserving
        // 2 - model in [3] Wang et al.

        double _norm = 0., _b0 = 0.; // Regge factor parameters: normalization and t-slope
        regge_trajectory * _traj;

        // Photon - Vector - Pomeron vertex
        std::complex<double> top_vertex(int mu);

        // Nucleon - Nucleon - Pomeron vertex
        std::complex<double> bottom_vertex(int mu);

        // Energy dependence from Pomeron propogator
        std::complex<double> regge_factor();
    };
};

#endif
