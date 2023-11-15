#include "EventGenerator.hpp"
#include "Experiment.hpp"
#include "GlueX.hpp"
#include "jpacAmplitude.hpp"
#include "ThreeBodyWriter.hpp"
#include "DataReader.hpp"
#include "plotter.hpp"
#include "histogram.hpp"
#include "cgamma.hpp"

#include "double_regge/model_30.hpp"

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
    };
};

void generator()
{
    using namespace jpacPhoto;
    using R = ThreeBodyReader;
    
    std::string outfile = "results";

    // ------------------------------------------------------------------------------
    // Generate phasespace

    GlueX gluex(8.0, 9.0);
    Monoenergetic mono(8.5);
    
    EventGenerator<ThreeBodyWriter,3> generator( {M_ETA, M_PION, M_PROTON}, {"eta", "pi", "p"} );
    generator.setExperiment(&gluex);

    // Weighted MC
    // generator.generatePhaseSpace(1E7, outfile + ".root");
    generator.generateWeightedPhysics<two_meson::jpacAmplitude<double_regge>>(1E5, "etapi.cfg", outfile + ".root"); 

    // Hit-or-miss MC
    // generator.intensityProfile<two_meson::jpacAmplitude<double_regge>>("etapi.cfg", "intensities.pdf");
    // generator.generatePhysics<two_meson::jpacAmplitude<double_regge>>(1E6, "etapi.cfg", outfile + ".root", 0.015); 

    // ------------------------------------------------------------------------------
    // Make histograms

    ThreeBodyReader reader({"eta", "pi", "p", outfile + ".root"});

    plotter plotter;
    histogram_1D ht      = plotter.new_histogram_1D("#it{t}  [GeV^{2}]", {-6, 1});
    histogram_1D hs12    = plotter.new_histogram_1D("#it{s}_{#eta#pi}  [GeV^{2}]", {0, 14});
    histogram_1D hs1     = plotter.new_histogram_1D("#it{s}_{#eta#it{p}}  [GeV^{2}]", {0, 14});
    histogram_1D hs2     = plotter.new_histogram_1D("#it{s}_{#pi#it{p}}  [GeV^{2}]", {0, 14});
    histogram_1D ht1     = plotter.new_histogram_1D("#it{t}_{#eta}  [GeV^{2}]", {-10, 1});
    histogram_1D ht2     = plotter.new_histogram_1D("#it{t}_{#pi}  [GeV^{2}]", {-10, 1});
    histogram_1D hcos    = plotter.new_histogram_1D("cos#theta_{GJ}", {-1.,1});
    histogram_1D hphi    = plotter.new_histogram_1D("#phi_{GJ}", {-PI, PI});
    histogram_1D hEbeam  = plotter.new_histogram_1D("#it{E}_{#gamma}", {7.5, 9.5});

    histogram_2D phi_s12 = plotter.new_histogram_2D("#it{s}_{#eta#pi}  [GeV^{2}]",    "#phi_{GJ}");
    histogram_2D phi_s1  = plotter.new_histogram_2D("#it{s}_{#eta#it{p}}  [GeV^{2}]", "#phi_{GJ}");
    histogram_2D phi_s2  = plotter.new_histogram_2D("#it{s}_{#pi#it{p}}  [GeV^{2}]",  "#phi_{GJ}");

    for (int i = 0; i < reader.numEvents(); i++)
    {
        auto kin = reader.getEvent();
        double w = kin->weight();
        
        double setapi = reader.getExtra(R::kS12);
        double spip   = reader.getExtra(R::kS2);

        if (setapi < 4.0 || spip < 4.0) continue;
        ht.fill(   reader.getExtra(R::kT)    , w);
        hs12.fill( reader.getExtra(R::kS12)  , w);
        hs1.fill(  reader.getExtra(R::kS1)   , w);
        hs2.fill(  reader.getExtra(R::kS2)   , w);
        ht1.fill(  reader.getExtra(R::kT1)   , w);
        ht2.fill(  reader.getExtra(R::kT2)   , w);
        hcos.fill( reader.getExtra(R::kCosGJ), w);
        hphi.fill( reader.getExtra(R::kPhiGJ), w);
        hEbeam.fill( kin->particle(0).E(), w);

        phi_s12.fill( reader.getExtra(R::kS12), reader.getExtra(R::kPhiGJ), w);
        phi_s1.fill(  reader.getExtra(R::kS1),  reader.getExtra(R::kPhiGJ), w);
        phi_s2.fill(  reader.getExtra(R::kS2),  reader.getExtra(R::kPhiGJ), w);
    };
    plotter.combine({3,4}, {&hphi, &hcos, &hs12, &hs1, &hs2, &ht1, &ht2, &ht, &hEbeam, &phi_s1, &phi_s2, &phi_s12}, outfile + ".pdf");
};