// This is a template file to illustrate basic requirements to define a new amplitude 
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef TEMPLATE_HPP
#define TEMPLATE_HPP

// Basic includes
#include "constants.hpp"
#include "kinematics.hpp"
#include "amplitude.hpp"

// Amplitudes should be added to the jpacPhoto namespace
namespace jpacPhoto
{
    // Additionally a second namespace may be used to differentiate different
    // implementations of the same amplitude
    // e.g. covariant::my_amplitude and analytic::my_amplitude

    // All amplitude implementations must derive from raw_amplitude
    class my_amplitude : public raw_amplitude
    {
        public: 

        // Basic constructor
        my_amplitude(amplitude_key key, kinematics xkinem, std::string id = "my_amplitude's id")
        : raw_amplitude(key, xkinem, "my_amplitude's name", id)
        {
            // Constructor should initialize with number of parameters N_par
            initialize(N_par);
        };

        // Additional constructors can be used with up to three additional parameters 
        // used as local data of arbitrary type
        my_amplitude(amplitude_key key, kinematics xkinem, arbitrarytype local_data, std::string id = "my_amplitude's id")
        : raw_amplitude(key, xkinem, "my_amplitude's name", id)
        {
            initialize(N_par);

            // Do something with local_data
        };

        // -----------------------------------------------------------------------
        // REQUIRED VIRTUAL FUNCTIONS

        // Provide a way to compute the amplitude at fixed helicities, s, and t
        inline complex helicity_amplitude(std::array<int,4> helicities, double s, double t)
        {
            // jpacPhoto::raw_amplitude has protected member data
            // _helicities, _s, and _t
            // which can be accessed from other functions without needed to be passed as argument
            // the store function saves them;
            store(helicities, s, t);

            return /* compute amplitude */;
        };

        // Must specifiy which center-of-mass scattering channel the helicities are defined
        // This is required to correctly give parity-relations and polarization observables
        inline helicity_channel native_helicity_frame(){ return helicity_channel::S_CHANNEL; };

        // Specify which final state particles amplitude can acommodate
        inline std::vector<std::array<int,2>> allowed_meson_JP() { return { VECTOR }; };
        inline std::vector<std::array<int,2>> allowed_baryon_JP(){ return { HALFPLUS }; };

        // Given a std::vector<double> of size N_pars
        // save parameters to local variables with which to calculate amplitude
        // The parameters are set by user script with raw_amplitude::set_parameters() which 
        // will do a check the input is of appropriate size and then call 
        // raw_amplitude::allocate_parameters internally
        inline void allocate_parameters(std::vector<double> pars)
        {
            /* save pars somewhere */
        };

        // -----------------------------------------------------------------------
        // OPTIONAL VIRTUAL FUNCTIONS

        // Assign each parameter a name, useful for fitting utlities
        inline std::vector<std::string> parameter_labels()
        {
            return {"p1", "p2", etc... };
        };
        
        // Use amplitude_option's to turn on and off features of your amplitude
        inline void set_option(amplitude_option x)
        {
            if (x == ExponentialFF) /* do something*/;
            else if (x == MonopoleFF) /* do something else */;
        };
    };
};

#endif