// Wrapper for a jpacPhoto amplitude to be used in AmpTools
// Because all member data needs to be provided at construction, we sneak in all 
// the customization of the jpacPhoto component through a template which provides
// static functions!
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#ifndef PRODUCTION_AMPLITUDE_HPP
#define PRODUCTION_AMPLITUDE_HPP

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/AmpParameter.h"

#include "constants.hpp"
#include "amplitude.hpp"
#include <string>

namespace jpacPhoto
{
    namespace one_meson
    {
        // Here the template A must be a struct with the following STATIC
        // functions:
        //  
        // jpacPhoto::amplitude initialize_amplitude()
        //   - Return an amplitude pointer of choice already set up with kinematics, sums, etc
        //
        // unsigned int N_variables();
        //   - Return number of kinematic variables needed to calculate intensity
        //
        // void calculate_variables(GDouble** pKin, GDouble* userVars)
        //   - Take in the 4-vector and calculate the N_variables quantities into userVars
        //
        //  double intensity( GDouble* userVars, jpacPhoto::amplitude)
        //    - from the precalculated variables userVars and the initialized amplitude
        //      calculate the intensity function

        std::vector<double> convertParameters( std::vector<AmpParameter> input)
        {
            std::vector<double> output;
            for (auto par : input)
            {
                output.push_back( par );
            }
            return output;
        };

        template<class A>
        class jpacAmplitude : public UserAmplitude< jpacAmplitude<A> >
        {
            public:

            // -----------------------------------------------------------------
            // Required AmpTools Methods

            jpacAmplitude()
            : UserAmplitude< jpacAmplitude<A> >()
            {};

            jpacAmplitude( const vector< std::string>& args ) 
            : UserAmplitude< jpacAmplitude<A> >(args)
            {
                _amplitude = A::initialize_amplitude();

                // All others are assmed to be amplitude parameters
                for (int i = 0; i < args.size(); i++)
                {
                    _pars.push_back( AmpParameter(args[i]) );
                    this->registerParameter( _pars[i] );
                };

                // Pass parameters to amplitude
                _amplitude->set_parameters( convertParameters(_pars) );
            };

            ~jpacAmplitude(){};

            // For simplicity the name always correspond to jpacAmplitude
            // Details of sub amplitudes can be gleamed from jpacPhoto::amplitude::get_id()
            inline std::string name() const { return "jpacAmplitude"; };

            // We explicitly want to only calculate amplitudes using the userVars 
            inline bool needsUserVarsOnly() const { return true; }

            // A single jpacAmplitude can already have a sum included so we only need to consider
            // so we only need to consider a single AmpTools ampltiude. This means the intensity
            // can be calculated on the jpacPhoto side and square rooted 
            inline std::complex<GDouble> calcAmplitude( GDouble** pKin, GDouble* userVars ) const 
            { 
                return sqrt( A::intensity(userVars, _amplitude) );               
            };

            inline unsigned int numUserVars()                                     const { return A::N_variables(); };
            inline void         calcUserVars( GDouble** pKin, GDouble* userVars ) const { A::calculate_variables(pKin, userVars); };
            inline double       Intensity(    GDouble** pKin, GDouble* userVars ) const { return A::intensity( userVars, _amplitude); };

            protected:

            // This function must be overwritten to provide exact details 
            // of the dynamical model provided by jpacPhoto
            jpacPhoto::one_meson::amplitude  _amplitude;

            // Parameters which are fed to the amplitude object
            std::vector<AmpParameter> _pars;
        };
    };
};

#endif