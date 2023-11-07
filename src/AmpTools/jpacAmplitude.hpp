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
#include "IUAmpTools/Kinematics.h"
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
        //  double intensity( GDouble* userVars, jpacPhoto::amplitude)
        //    - from the precalculated variables userVars and the initialized amplitude
        //      calculate the intensity function

       std::vector<double> convertParameters( std::vector<AmpParameter> input)
        {
            std::vector<double> output;
            for (auto par : input)
            {
                output.push_back( double(par) );
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

            jpacAmplitude( const std::vector<std::string>& args ) 
            : UserAmplitude< jpacAmplitude<A> >(args)
            {
                _amplitude = A::initialize_amplitude();

                for (int i = 0; i < args.size(); i++)
                {
                    _pars.push_back( AmpParameter(args[i]) );
                };

                for (int i = 0; i < args.size(); i++)
                {
                    this->registerParameter( _pars[i] );
                }
            };

            ~jpacAmplitude(){};

            // For simplicity the name always correspond to jpacAmplitude
            // Details of sub amplitudes can be gleamed from jpacPhoto::amplitude::get_id()
            inline std::string name() const { return "jpacAmplitude"; };

            // We explicitly want to only calculate amplitudes using the userVars 
            inline bool   needsUserVarsOnly() const { return true; }
            inline unsigned int numUserVars() const { return 2; };
            enum UsersVars { kVar_s = 0, kVar_t = 1 };

            inline void calcUserVars( GDouble** pKin, GDouble* userVars ) const
            {
                TLorentzVector meson( pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0]);
                TLorentzVector baryon(pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0]);

                double s = (meson + baryon).M2();
                userVars[kVar_s] = s;

                double Egam = jpacPhoto::E_beam(sqrt(s));
                TLorentzVector beam( 0, 0, Egam, Egam);
                double t = (beam - meson).M2();
                userVars[kVar_t] = t;
            };

            // A single jpacAmplitude can already have a sum included so we only need to consider
            // so we only need to consider a single AmpTools ampltiude. This means the intensity
            // can be calculated on the jpacPhoto side and square rooted 
            inline std::complex<GDouble> calcAmplitude( GDouble** pKin, GDouble* userVars ) const 
            { 
                // Make sure the amplitude knows about the most up-to-date parameters in AmpTools
                _amplitude->set_parameters( convertParameters(_pars) );
                return sqrt( A::intensity(userVars[kVar_s], userVars[kVar_t], _amplitude) );               
            };
            inline double Intensity(GDouble** pKin, GDouble* userVars ) const { return A::intensity(userVars[kVar_s], userVars[kVar_t], _amplitude); };

            protected:

            // This function must be overwritten to provide exact details 
            // of the dynamical model provided by jpacPhoto
            jpacPhoto::one_meson::amplitude  _amplitude;

            // Parameters which are fed to the amplitude object
            std::vector<AmpParameter> _pars;
        };
    };

    namespace two_meson
    {
        // Here the template A must be a struct with the following STATIC
        // functions:
        //  
        // std::string name()
        // - Return a name identifier which will be used throughout AmpTools
        // - equivalent of UserAmplitude<A>::name()
        //
        // jpacPhoto::amplitude initialize_amplitude()
        //   - Return an amplitude pointer of choice already set up with kinematics, sums, etc
        //
        //  double intensity( GDouble* userVars, jpacPhoto::amplitude)
        //    - from the precalculated variables userVars and the initialized amplitude
        //      calculate the intensity function

       std::vector<double> convertParameters( std::vector<AmpParameter> input)
        {
            std::vector<double> output;
            for (auto par : input)
            {
                output.push_back( double(par) );
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

            jpacAmplitude( const std::vector<std::string>& args ) 
            : UserAmplitude< jpacAmplitude<A> >(args)
            {
                _amplitude = A::initialize_amplitude();

                for (int i = 0; i < args.size(); i++)
                {
                    _pars.push_back( AmpParameter(args[i]) );
                };

                for (int i = 0; i < args.size(); i++)
                {
                    this->registerParameter( _pars[i] );
                }
            };

            ~jpacAmplitude(){};

            // For simplicity the name always correspond to jpacAmplitude
            // Details of sub amplitudes can be gleamed from jpacPhoto::amplitude::get_id()
            inline std::string name() const { return A::name(); };

            // We explicitly want to only calculate amplitudes using the userVars 
            inline bool   needsUserVarsOnly() const { return true; }
            inline unsigned int numUserVars() const { return 5; };
            enum UsersVars { kVar_s = 0, kVar_t = 1, kVar_s12 = 2, kVar_thetaGJ = 3, kVar_phiGJ = 4 };
         
            inline void calcUserVars( GDouble** pKin, GDouble* userVars ) const
            {
                // Inital state reconstructed from the internally calcualted _Ebeam
                TLorentzVector m1(    pKin[0][1], pKin[0][2], pKin[0][3], pKin[0][0]);
                TLorentzVector m2(    pKin[1][1], pKin[1][2], pKin[1][3], pKin[1][0]);
                TLorentzVector recoil(pKin[2][1], pKin[2][2], pKin[2][3], pKin[2][0]);
                userVars[kVar_s]   = (m1 + m2 + recoil).M2();

                double Ebeam = E_beam( sqrt(userVars[kVar_s]) );
                TLorentzVector beam(  0., 0., Ebeam, Ebeam);
                TLorentzVector target(0., 0., 0., M_PROTON);

                userVars[kVar_t]   = (recoil - target).M2();
                userVars[kVar_s12] = (m1 + m2).M2();

                // Boost only the initalstate and particle 1 to GJ frame to calculate angles
                beam.Boost(   -(m1 + m2).BoostVector());
                recoil.Boost( -(m1 + m2).BoostVector());
                m1.Boost(     -(m1 + m2).BoostVector());

                TVector3 z = beam.Vect().Unit();
                TVector3 y = (beam.Vect().Cross(recoil.Vect())).Unit();
                TVector3 x = y.Cross(z);

                TVector3 GJ( m1.Vect().Dot(x), m1.Vect().Dot(y), m1.Vect().Dot(z));
                userVars[kVar_thetaGJ] = GJ.Theta(); userVars[kVar_phiGJ] = GJ.Phi(); 
            };

            // A single jpacAmplitude can already have a sum included so we only need to consider
            // so we only need to consider a single AmpTools ampltiude. This means the intensity
            // can be calculated on the jpacPhoto side and square rooted 
            inline std::complex<GDouble> calcAmplitude( GDouble** pKin, GDouble* userVars ) const 
            { 
                // Make sure the amplitude knows about the most up-to-date parameters in AmpTools
                _amplitude->set_parameters( convertParameters(_pars) );
                return sqrt( A::intensity(userVars[kVar_s], userVars[kVar_t], userVars[kVar_s12], userVars[kVar_thetaGJ], userVars[kVar_phiGJ], _amplitude) );               
            };
            inline double Intensity(GDouble** pKin, GDouble* userVars ) const 
            {
                return A::intensity(userVars[kVar_s], userVars[kVar_t], userVars[kVar_s12], userVars[kVar_thetaGJ], userVars[kVar_phiGJ], _amplitude); 
            };

            protected:

            // This function must be overwritten to provide exact details 
            // of the dynamical model provided by jpacPhoto
            jpacPhoto::two_meson::amplitude  _amplitude;

            // Parameters which are fed to the amplitude object
            std::vector<AmpParameter> _pars;
        };
    };
};

#endif