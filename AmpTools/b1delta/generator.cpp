#include "EventGenerator.hpp"
#include "jpacAmplitude.hpp"
#include "TwoBodyWriter.hpp"
#include "regge/unnatural_exchange.hpp"
#include "regge/natural_exchange.hpp"
#include <string>
#include <cstdlib>

namespace jpacPhoto
{
    struct stable_b1delta
    {
        inline static one_meson::amplitude initialize_amplitude()
        {
            using namespace jpacPhoto::one_meson;

            kinematics kB1Delta = new_kinematics( M_B1, M_DELTA );
            kB1Delta->set_meson_JP( 1, +1);
            kB1Delta->set_baryon_JP(3, +1);

            amplitude pi  = new_amplitude<regge::unnatural_exchange>(kB1Delta, "pi  exchange");
            amplitude rho = new_amplitude<regge::natural_exchange>(  kB1Delta, "rho exchange");
            return pi + rho;
        };

        // Use the differential cross-section as unpolarized intensity
        inline static double intensity( double s, double t, jpacPhoto::one_meson::amplitude amp) 
        { 
            return amp->differential_xsection(s, t);               
        };
    };
};

void generator()
{
    using namespace jpacPhoto;

    std::string cfgname = "b1delta.cfg";

    double Egam = 9.;

    EventGenerator<TwoBodyWriter,2> generator( {M_B1, M_DELTA}, {"b1", "delta"} );
    generator.setBeamEnergy(Egam);

    generator.generatePhaseSpace(1E5, "phase_space.root");
    generator.generatePhysics< one_meson::jpacAmplitude<stable_b1delta> >(1E5, cfgname, "physics.root");
};