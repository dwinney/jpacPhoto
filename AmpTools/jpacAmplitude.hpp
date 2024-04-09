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

#ifndef JPAC_AMPLITUDE_HPP
#define JPAC_AMPLITUDE_HPP

#include "IUAmpTools/Amplitude.h"
#include "IUAmpTools/UserAmplitude.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/AmpParameter.h"

#include "constants.hpp"
#include "amplitude.hpp"
#include <string>

namespace jpacPhoto
{
    std::vector<double> convertParameters( std::vector<AmpParameter> input)
    {
        std::vector<double> output;
        for (auto par : input) output.push_back( double(par) );
        return output;
    };

    // Here the template A must be a struct with the following STATIC
    // functions:
    //  
    // jpacPhoto::amplitude initialize_amplitude()
    //   - Return an amplitude pointer of choice already set up with kinematics, sums, etc

    template<class A>
    class jpacAmplitude : public UserAmplitude< jpacAmplitude<A> >
    {
        public:

        // -----------------------------------------------------------------
        // Required AmpTools Methods

        jpacAmplitude()
        : UserAmplitude< jpacAmplitude<A> >()
        {};

        jpacAmplitude( const std::vector<std::string>& args ) 
        : UserAmplitude< jpacAmplitude<A> >(args)
        {
            _amplitude = A::initAmplitude();

            for (int i = 0; i < args.size(); i++) _pars.push_back( AmpParameter(args[i]) );
            for (int i = 0; i < args.size(); i++) this->registerParameter( _pars[i] );
        };

        ~jpacAmplitude(){};

        // For simplicity the name always correspond to jpacAmplitude
        // Details of sub amplitudes can be gleamed from jpacPhoto::amplitude::get_id()
        inline std::string name() const { return A::ampName(); };

        // If we explicitly want to only calculate amplitudes using the userVars 
        inline bool   needsUserVarsOnly() const { return A::needsUserVarsOnly(); }
        inline unsigned int numUserVars() const { return A::numUserVars(); };

        inline void calcUserVars( GDouble** pKin, GDouble* userVars ) const { A::calcUserVars(pKin, userVars); };

        // A single jpacAmplitude can already have a sum included so we only need to consider
        // so we only need to consider a single AmpTools ampltiude. This means the intensity
        // can be calculated on the jpacPhoto side and square rooted 
        inline std::complex<GDouble> calcAmplitude( GDouble** pKin, GDouble* userVars ) const 
        { 
            // Make sure the amplitude knows about the most up-to-date parameters in AmpTools
            _amplitude->set_parameters( convertParameters(_pars) );
            return A::calcAmplitude(pKin, userVars, _amplitude);
        };

        protected:

        // This function must be overwritten to provide exact details 
        // of the dynamical model provided by jpacPhoto
        amplitude  _amplitude;

        // Parameters which are fed to the amplitude object
        std::vector<AmpParameter> _pars;
    };
};

#endif