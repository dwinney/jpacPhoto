#include "EventGenerator.hpp"
#include "Experiment.hpp"
#include "DataWriter.hpp"
#include "jpacAmplitude.hpp"
#include "analytic/pomeron_exchange.hpp"
#include <string>
#include <cstdlib>

namespace jpacPhoto
{
    struct decayless_jpsip
    {
        inline static std::string ampName(){ return "jpacAmplitude<pomeron_exchange>"; };

        inline static bool needsUserVarsOnly()  { return true; };
        inline static unsigned int numUserVars(){ return 2;    };

        enum UsersVars { kVar_s = 0, kVar_t = 1 };

        inline static amplitude initAmplitude()
        {
            kinematics kpsi = new_kinematics( M_JPSI, M_PROTON );
            kpsi->set_meson_JP( 1, -1);
            kpsi->set_baryon_JP(1, +1);

            amplitude pomeron  = new_amplitude<analytic::pomeron_exchange>(kpsi, "Pomeron exchange");
            return pomeron;
        };
        
        // For 2->2 we only have s and t independent variables
        inline static void calcUserVars( GDouble** pKin, GDouble* userVars)
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

        // We have a single amplitude, so just use the square-root of spin-summed amplitude
        inline static std::complex<GDouble> calcAmplitude( GDouble ** pKin, GDouble* userVars, amplitude amp)
        {
            return sqrt(amp->probability_distribution( userVars[kVar_s], userVars[kVar_t]));
        };
    };
};

void generator()
{
    using namespace jpacPhoto;

    EventGenerator<TwoBodyWriter,2> generator( {M_JPSI, M_PROTON}, {"jpsi", "p"} );

    Monoenergetic egam(10);
    generator.setExperiment(&egam);
    
    generator.generatePhaseSpace(1E5, "phasespace.root");
    generator.generatePhysics< jpacAmplitude<decayless_jpsip> >(1E5, "pomeron.cfg", "physics.root");
};