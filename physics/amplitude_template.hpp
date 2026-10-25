// This is a template file to illustrate basic requirements to define a new amplitude 
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef AMPLITUDE_TEMPLATE_HPP
#define AMPLITUDE_TEMPLATE_HPP

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
    // Every amplitude is carries a string id in the constructor to identify it from other amplitudes when plotting/fitting/error handling
    class my_amplitude : public raw_amplitude
    {
        public: 

        // Basic constructor
        my_amplitude(key k, kinematics xkinem)
        : raw_amplitude(k, xkinem, "default_id")
        {
            // Constructor should initialize with number of parameters
            // this specifies 2 free parameters
            initialize(2);

            // Do anything else you want at construction
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

            complex result;
            /* compute amplitude */

            return result;
        };

        // Must specifiy which center-of-mass scattering channel the helicities are defined
        // This is required to correctly give parity-relations and polarization observables
        inline helicity_channel native_helicity_frame(){ return helicity_channel::S_CHANNEL; };

        // Specify which final state particles amplitude can acommodate
        inline std::vector<particle> allowed_mesons() { return { VECTOR, AXIALVECTOR }; };
        inline std::vector<particle> allowed_baryons(){ return { HALFPLUS, HALFMINUS }; };

        // Given a std::vector<double> of size N_pars
        // save parameters to local variables with which to calculate amplitude
        // The parameters are set by user script with raw_amplitude::set_parameters() which 
        // will do a check the input is of appropriate size and then call 
        // raw_amplitude::allocate_parameters internally
        inline void allocate_parameters(std::vector<double> pars)
        {
            /* save pars somewhere */
            _p1 = pars[0];
            _p2 = pars[1];
        };

        // -----------------------------------------------------------------------
        // OPTIONAL VIRTUAL FUNCTIONS

        // Assign each parameter a name, useful for fitting utlities
        inline std::vector<std::string> parameter_labels()
        {
            return {"p1", "p2" /*, etc...*/ };
        };
        
        // The set_option method allows you to use flags to change different things of 
        // the amplitude. For human-readability its useful to define static consts for
        // each relevant option
        static const int kOption1 = 0;
        static const int kOption2 = 1;
        inline void set_option(int x)
        {
            switch (x)
            {
                case kOption1: { /* do something */      break; };
                case kOption2: { /* do something else */ break; };
                default: return; // do nothing
            };
        };

        // -----------------------------------------------------------------------
        // Other than the required virtual functions structure can be whatever you want

        inline void new_function()
        {
            /* do something not incorporated in the parent raw_amplitude class */
        };

        // Note that amplitudes can only be created through the make_amplitude<T> method
        // This means that public non-virtual functions arent accessable from an amplitude object
        // To access class-specific member data you need to downcast, e.g.:

        // amplitude amp = make_amplitude<my_amplitude>(kinem);
        // auto downcasted_amp = std::dynamic_pointer_cast<my_amplitude>(amp);
        // downcasted_amp->new_function();

        protected:

        // Local data
        double _p1 = 0, _p2 = 0;
    };
};

#endif