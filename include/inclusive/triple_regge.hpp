// Form of the triple Regge interaction using JPAC's parameterization
// i.e. using the t dependence from properly normalized Regge propagators
// and M2 dependence from the total hadronic cross-section of the bottom vertex.
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef TRIPLE_JPAC
#define TRIPLE_JPAC

#include <complex>
#include <tuple>
#include <functional>

#include "misc_math.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "inclusive/inclusive_production.hpp"
#include "inclusive/inclusive_kinematics.hpp"
#include "inclusive/total_xsection.hpp"

namespace jpacPhoto
{
    class triple_regge : public inclusive_production
    {
        public:

        triple_regge(pseudoscalar_exchange * exclusive)
        : inclusive_production(exclusive->_kinematics->_mX, exclusive->_identifier),
         _exclusive(exclusive),
         _useRegge(exclusive->if_reggeized()),
         _exchange_mass2(exclusive->get_mEx2()), _trajectory(exclusive->get_trajectory()),
         _b(exclusive->get_cutoff()), _g(exclusive->get_coupling())
        {
            initialize(exclusive->amplitude_name());
        };

        ~triple_regge()
        {
            delete _sigma_tot;
        };

        // Set flag to use the simplified high-energy approximation kinematics
        inline void set_high_energy_approximation(bool ifuse)
        {
            // Why have two flags if we only use (t, x) in the high energy approxiation?
            // I envision setting up the full (t, x) kinematics one day... will this happen? who knows. 
            _useTX = ifuse;
        };

        // Evaluate the invariant amplitude
        double d3sigma_d3p(double s, double t, double mm);

        // Change the stored sigma_total to another one
        void set_sigma_total(sigma_option opt);

        // ---------------------------------------------------------------------------
        // Individual pieces that go into the amplitude

        protected:

        // Pointer to an already set-up exclusive amplitude object
        amplitude * _exclusive;

        // Initialize the amplitude appropriate to the given exclusive amplitude
        // This sets the coupling function for the top vertex
        // and the default choice of sigma_tot for the bottom vertex
        void initialize(std::string exclusive_amp_name);

        // Coupling function associated with the t dependence of the top vertex
        std::function<double(double)> _coupling;

        // Whether or not to use regge or fixed-spin propagators        
        bool _useRegge = true;

        // Total cross-section associated with the bottom vertex
        total_xsection * _sigma_tot = NULL;

        // Regge trajectory for the exchange if Reggeized
        regge_trajectory * _trajectory = NULL;

        // Else we have a fixed pole-mass
        double _exchange_mass2 = 0.;

        double _s0 = 1.; // scale factor
        double _b  = 0.; // t-slope parameter in form-factor
        double _g  = 0.; // gamma coupling constant
    };
};

#endif