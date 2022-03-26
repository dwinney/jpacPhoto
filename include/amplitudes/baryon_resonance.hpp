// Parameterization of a resonant amplitude in the s-channel
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _RESONANCE_
#define _RESONANCE_

#include "amplitude.hpp"

// ---------------------------------------------------------------------------
// baryon_resonance class describes the amplitude corresponding to a narrow
// (Breit-Wigner) resonance in the channel. It is parameterized in terms of
// 3 functions:
//
// 1. Breit-Wigner pole for given mass and width
// 2. Hadronic decay coupling to J/psi p final state
// 3. Photo-excitation coupling to the gamma p initial state
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
    class baryon_resonance : public amplitude
    {
        public:
        // Constructor
        baryon_resonance(reaction_kinematics * xkinem, int j, int p, double mass, double width, std::string name = "baryon_resonance")
        : amplitude(xkinem, "baryon_resonance", name),
          _mRes(mass), _gamRes(width), 
          _resJ(j), _resP(p), _naturality(p * pow(-1, (j-1)/2))
        {
            set_nParams(2);
            check_JP(xkinem->_jp);

            // save momentum and other J^P dependent quantities
            _pibar = real(xkinem->_initial_state->momentum(mass * mass));
            _pfbar = real(xkinem->_final_state->momentum(mass * mass));

            if (abs(p) != 1)
            {
                std::cout << "Invalid parity " << p << " passed to " << name << ". Quitting...\n";
                exit(0);
            };
            switch (p * j)
            {
                case  1: {_lmin = 0; _pt = 2./3.; break;}
                case -1: {_lmin = 1; _pt = 3./5.; break;}
                case  3: {_lmin = 1; _pt = 3./5.; break;}
                case -3: {_lmin = 0; _pt = 2./3.; break;}
                case  5: {_lmin = 1; _pt = 3./5.; break;}
                case -5: {_lmin = 2; _pt = 1./3.; break;}
            
                default:
                {
                std::cout << "\nbaryon_resonance: spin-parity combination for J = " << j << "/2 and p = " << p << " not available. ";
                std::cout << "Quiting... \n";
                exit(0);
                }
            };
        };

        // Setting utility
        void set_params(std::vector<double> params)
        {
            check_nParams(params);
            _xBR = params[0];
            _photoR = params[1];
        };

        // Combined total amplitude including Breit Wigner pole
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t);

        // only vector kinematics allowed
        inline std::vector<std::array<int,2>> allowedJP()
        {
            return {{1, -1}};
        };
        
        inline int parity_phase(std::array<int, 4> helicities)
        {
            return _kinematics->parity_phase(helicities, HELICITY_CHANNEL::S);
        };

        private:

        // Photoexcitation helicity amplitude for the process gamma p -> R
        std::complex<double> photo_coupling(int lam_i);

        // Hadronic decay helicity amplitude for the R -> J/psi p process
        std::complex<double> hadronic_coupling(int lam_f);

        // Ad-hoc threshold factor to kill the resonance at threshold
        double threshold_factor(double beta);

        int _resJ, _resP, _naturality; // (2xSpin) and parity of the resonance
        double _mRes, _gamRes; // Resonant mass and width

        int _lmin; // lowest allowed relative angular momentum
        double _pt; // Combinatorial factor due to only transverse polarized J/psi contribute

        // Couplings
        double _xBR; // Hadronic banching fraction to j/psi p
        double _photoR; // Photocoupling ratio

        // Initial and final CoM momenta evaluated at resonance energy.
        double _pibar, _pfbar;

        // saved energies and angle
        double _s, _t, _theta;
    };
};
#endif
