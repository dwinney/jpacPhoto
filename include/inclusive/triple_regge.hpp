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
#include "inclusive/sigma_tot.hpp"

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
         _b(exclusive->get_cutoff())
        {
            _coupling = [&](double t)
            {
                return  (0.24 / M_B1) * (t - M2_B1);
            };

            _sigma_tot = new sigma_tot_PDG(M_PION, M_PROTON, {-1., 1., 9.56, 1.767, 18.75}, "rpp2020-pimp_total.dat");
        };

        ~triple_regge()
        {
            delete _sigma_tot;
        };

        // Fully parameterized constructor for the usual triple regge form
        // // Requires:
        // // 1. a kinematics object
        // // 2. a regge trajectory (assumed to be the same for the top vertices)
        // // Optional string identifier to associate with the object
        // triple_regge(double mass, regge_trajectory * alpha, std::string id = "")
        // : inclusive_production(mass, id), _trajectory(alpha)
        // {
        //     _useRegge = true;
        // };

        // triple_regge(double mass, double exmass, std::string id = "")
        // : inclusive_production(mass, id), _exchange_mass2(exmass*exmass)
        // {
        //     _useRegge = false;
        // };

        // // Set the top vertex coupling
        // inline void set_coupling(const std::function<double(double)> coupling)
        // {
        //     _coupling = coupling;
        //     _couplingSet = true;
        // };

        // // Set the total cross-section associated with the bottom vertex
        // inline void set_sigma_tot(sigma_tot * sigma)
        // {
        //     _sigma_tot = sigma;
        // };
        
        // // Set the t-slope associated with the exponential form-factor
        // inline void set_form_factor(double b)
        // {
        //     _b = b;
        // };

        // Set flag to use the simplified high-energy approximation kinematics
        inline void set_high_energy_approximation(bool ifuse)
        {
            // Why have two flags if we only use (t, x) in the high energy approxiation?
            // I envision setting up the full (t, x) kinematics one day... will this happen? who knows. 
            _useTX = ifuse;
        };

        // Evaluate the invariant amplitude
        double d3sigma_d3p(double s, double t, double mm);

        // ---------------------------------------------------------------------------
        // Individual pieces that go into the amplitude

        protected:

        // Pointer to an already set-up exclusive amplitude object
        amplitude * _exclusive;

        // Coupling function associated with the t dependence of the top vertex
        std::function<double(double)> _coupling;

        // Whether or not to use regge or fixed-spin propagators        
        bool _useRegge = true;

        // Total cross-section associated with the bottom vertex
        sigma_tot * _sigma_tot = NULL;

        // Regge trajectory for the exchange if Reggeized
        regge_trajectory * _trajectory = NULL;
        // Else we have a fixed pole-mass
        double _exchange_mass2 = 0.;

        double _s0 = 1.; // scale factor
        double _b  = 0.; // t-slope parameter in form-factor
    };
};

#endif