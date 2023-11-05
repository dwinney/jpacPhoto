#include "EventGenerator.hpp"
#include "jpacAmplitude.hpp"
#include "TwoBodyWriter.hpp"
#include <string>
#include <cstdlib>

#include "b1delta/stable_b1delta.hpp"

void generator()
{
    using namespace jpacPhoto;

    std::string cfgname = "b1delta.cfg";

    double Egam = 9.;

    EventGenerator<TwoBodyWriter,2> generator( {M_B1, M_DELTA}, {"b1", "delta"} );
    generator.setBeamEnergy(Egam);

    generator.generatePhaseSpace(1E5, "phase_space.root");
    generator.generatePhysics< one_meson::jpacAmplitude<stable_b1delta> >(1E7, cfgname, "physics.root");
};