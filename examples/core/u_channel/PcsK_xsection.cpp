#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "dirac_exchange.hpp"
#include "vector_exchange.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

void PcsK_xsection()
{
    // Set up kinematics of the Pcs K final state
    double M_PCS    = 4.5880;
    double gNucleon = -6.4469; 
    double cutoff   = 0.9;

    // Couplings 
    double g1  = 4.45E-3;
    double g3  = 3.01E-3;
    double g3p = 6.29E-3;

    // Spin half 
    reaction_kinematics kspinhalf (M_KAON, M_PCS);
    kspinhalf.set_meson_JP (0, -1);
    kspinhalf.set_baryon_JP(1, -1);

    dirac_exchange spinhalf(&kspinhalf, M_LAMBDA, "Spin-1/2");
    spinhalf.set_params({g1, gNucleon});
    spinhalf.set_formfactor(3, cutoff);

    // Spin three-halves
    reaction_kinematics kspinthreehalf (M_KAON, M_PCS);
    kspinthreehalf.set_meson_JP (0, -1);
    kspinthreehalf.set_baryon_JP(3, -1);

    dirac_exchange spinthreehalf(&kspinthreehalf, M_LAMBDA, "Spin-3/2");
    spinthreehalf.set_params({g3, gNucleon});
    spinthreehalf.set_formfactor(3, cutoff);

    // // ---------------------------------------------------------------------------

    int   N = 50;
    bool PRINT_TO_COMMANDLINE  = true;
    std::string filename  = "PcsK.pdf";

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    std::vector<amplitude*> amps;
    amps.push_back(&spinhalf);
    amps.push_back(&spinthreehalf);

    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    double  xmin = 5.;
    double  xmax = 12.;

    double  ymin = 2.E-3;
    double  ymax = 5.;

    std::string ylabel    = "#sigma(#gamma #it{p} #rightarrow #it{P}_{cs} #it{K})  [nb]";
    std::string xlabel    = "#it{W}_{#gamma #it{p}}  [GeV]";

    for (int i = 0; i < amps.size(); i++)
    {
        auto F = [&](double x)
        {
            return amps[i]->integrated_xsection(x*x);
        };

        plotter->AddEntry(N, F, {xmin, xmax}, amps[i]->get_id());
    }

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel,  ymin, ymax);
    plotter->SetYlogscale(1);
    plotter->RemoveLogo();

    plotter->SetLegend(0.7, 0.75);
    plotter->SetLegendOffset(0.3, 0.1);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
};