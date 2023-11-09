#include "DataReader.hpp"
#include "plotter.hpp"
#include "histogram.hpp"

void plot()
{
    using namespace jpacPhoto;
    using jpacPhoto::DataReader;
    using R = ThreeBodyReader;
    
    ThreeBodyReader reader({"eta", "pi", "p", "unweighted.root"});

    plotter plotter;
    histogram_1D ht   = plotter.new_histogram_1D("#it{t}  [GeV^{2}]");
    histogram_1D hs12 = plotter.new_histogram_1D("#it{s}_{#eta#pi}  [GeV^{2}]");
    histogram_1D hs1  = plotter.new_histogram_1D("#it{s}_{#eta}  [GeV^{2}]");
    histogram_1D hs2  = plotter.new_histogram_1D("#it{s}_{#pi}  [GeV^{2}]");
    histogram_1D ht1  = plotter.new_histogram_1D("#it{t}_{#eta}  [GeV^{2}]");

    histogram_1D ht2  = plotter.new_histogram_1D("#it{t}_{#pi}  [GeV^{2}]");
    histogram_1D hcos = plotter.new_histogram_1D("cos#theta_{GJ}");
    histogram_1D hphi = plotter.new_histogram_1D("#phi_{GJ}");

    histogram_2D phi_s12 = plotter.new_histogram_2D("s_{#eta#pi}  [GeV^{2}]", "#phi_{GJ}");
    histogram_2D phi_s1 = plotter.new_histogram_2D("s_{#eta}  [GeV^{2}]",     "#phi_{GJ}");
    histogram_2D phi_s2 = plotter.new_histogram_2D("s_{#pi}  [GeV^{2}]",      "#phi_{GJ}");

    for (int i = 0; i < reader.numEvents(); i++)
    {
        double w = reader.getEvent()->weight();
        
        ht.fill(   reader.getExtra(R::kT)    , w);
        hs12.fill( reader.getExtra(R::kS12)  , w);
        hs1.fill(  reader.getExtra(R::kS1)   , w);
        hs2.fill(  reader.getExtra(R::kS2)   , w);
        ht1.fill(  reader.getExtra(R::kT1)   , w);
        ht2.fill(  reader.getExtra(R::kT2)   , w);
        hcos.fill( reader.getExtra(R::kCosGJ), w);
        hphi.fill( reader.getExtra(R::kPhiGJ), w);

        phi_s12.fill( reader.getExtra(R::kS12), reader.getExtra(R::kPhiGJ), w);
        phi_s1.fill(  reader.getExtra(R::kS1),  reader.getExtra(R::kPhiGJ), w);
        phi_s2.fill(  reader.getExtra(R::kS2),  reader.getExtra(R::kPhiGJ), w);
    };
    plotter.combine({3,4}, {&ht, &ht1, &ht2, &hs12, &hs1, &hs2, &hcos, &hphi, nullptr, &phi_s12, &phi_s1, &phi_s2}, "etapi_plots.pdf");
};