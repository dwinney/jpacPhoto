#include "EventGenerator.hpp"
#include "jpacAmplitude.hpp"
#include "ThreeBodyWriter.hpp"
#include "double_regge/model_30.hpp"

#include <string>
#include <cstdlib>

namespace jpacPhoto
{
    struct double_regge
    {
        inline static std::string name(){ return "jpacAmplitude<Model_3.0>"; };

        inline static two_meson::amplitude initialize_amplitude()
        {
            using namespace jpacPhoto::two_meson;

            kinematics kEtaPi = new_kinematics( M_ETA, M_PION, M_PROTON );
            amplitude amp     = new_amplitude<model_30>(kEtaPi);
            amp->set_option(model_30::kFast1);
            return amp;
        };

        // Use the differential cross-section as unpolarized intensity
        inline static double intensity(double s, double t, double s12, double thetaGJ, double phiGJ, two_meson::amplitude amp) 
        { 
            return amp->differential_xsection(s, t, s12, thetaGJ, phiGJ);               
        };
    };
};

void generator()
{
    using namespace jpacPhoto;
    using namespace jpacPhoto::two_meson;

    EventGenerator<ThreeBodyWriter,3> generator( {M_ETA, M_PION, M_PROTON}, {"eta", "pi", "p"} );
    generator.setBeamEnergy(9.);

    // Pure phase space for gam p -> eta pi p
    // generator.generatePhaseSpace(1E5, "unweighted.root"); 
    // generator.generateWeightedPhaseSpace(1E5, "weighted.root"); 

    // Hit or miss MC using model
    // generator.generatePhysics<jpacAmplitude<double_regge>>(1E5, "etapi.cfg", "unweighted.root");      
    // Phase space but with saved weights calculated from amplitude
    generator.generateWeightedPhysics<jpacAmplitude<double_regge>>(1E6, "etapi.cfg", "weighted.root"); 
};