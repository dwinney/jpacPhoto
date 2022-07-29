// Spin-1/2 exchange ampltiude from perturbation theory
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _PROTON_
#define _PROTON_

#include <string>
#include <vector>
#include <iostream>

#include <iomanip>

#include "amplitude.hpp"

namespace jpacPhoto
{
    class dirac_exchange : public amplitude
    {
        public:
        
        // constructor
        dirac_exchange(reaction_kinematics * xkinem, double mass, std::string name = "dirac_exchange")
        : amplitude(xkinem, "dirac_exchange", name),
            _mEx(mass), _mEx2(mass*mass)
        {
            set_nParams(2);
            check_JP(xkinem);
        };

        // Setting utility
        void set_params(std::vector<double> params)
        {
            check_nParams(params);
            _gG = params[0];
            _gN = params[1];
        };

        // Whether or not to include an form factor (default false)
        // FF = 0 (none), 1 (exponential), 2 (monopole)
        inline void set_formfactor(int FF, double bb = 0.)
        {
            _useFF = FF;
            _cutoff = bb;
        }

        // Assemble the helicity amplitude by contracting the spinor indices
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t);

        // Helicities are always assumed to be in the s-channel cm frame
        inline helicity_channel helicity_CM_frame(){ return S; };

        // only Vector and psuedo-scalar kinematics
        inline std::vector<std::array<int,2>> allowed_meson_JP()
        {
            return { {1, -1}, {0, -1} };
        };
        inline std::vector<std::array<int,2>> allowed_baryon_JP()
        {
            return { {1, +1}, {1, -1} };
        };

        protected:
    
        // Exchange nucleon mass
        double _u;
        double _mEx, _mEx2;

        // Form factor parameters
        int _useFF = 0;
        double _cutoff = 0.;
        double form_factor();

        // couplings
        double _gG, _gN;

        // We only have the covariant amplitude evaulation here
        std::complex<double> covariant_amplitude();

        // Beam -- Exchange -- Produced Baryon
        std::complex<double> top_vertex(int i);
        std::complex<double> halfplus_coupling(int i);
        std::complex<double> halfminus_coupling(int i);

        // Targer -- Exchange -- Produced Meson
        std::complex<double> bottom_vertex(int j);
        std::complex<double> vector_coupling(int j);
        std::complex<double> pseudoscalar_coupling(int j);

        // Spin-1/2 propagator
        std::complex<double> dirac_propagator(int i, int j);
    };
};
#endif
