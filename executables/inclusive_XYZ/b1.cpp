#include "regge_trajectory.hpp"
#include "inclusive/triple_regge.hpp"
#include "inclusive/sigma_tot.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>
#include <memory>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // For all Z's we only have pion exchange
    double alphaPrime = 0.7, alpha0 = -M2_PION * alphaPrime;
    auto alphaPi = new saturated_trajectory(+1, alpha0, alphaPrime, "#pi trajectory");
    alphaPi->set_minimum_spin(0);

    double LamPi = .9;  // 900 MeV cutoff for formfactor

    // ---------------------------------------------------------------------------
    // b1 meson
    // ---------------------------------------------------------------------------

    auto   b1  = new triple_regge(M_B1, "b_{1}(1235)");
    double gb1 = 0.24; // Gamma - Pi - b1 coupling

    auto beta_b1 = [&](double t)
    {
        return sqrt(2.) * sqrt(alphaPrime) * (gb1 / M_B1) * sqrt(Kallen(t, 0., M2_B1));
    };

    b1->add_term(alphaPi, beta_b1, &sigma_tot_pipp);
    b1->set_formfactor(true, LamPi);

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<triple_regge*> amps;
    amps.push_back(b1);

    // Plotter object
    auto plotter = new jpacGraph1D();

    int N = 100;
    double xmin = 20.;
    double xmax = 100.;

    double ymin = 0.;
    double ymax = 1.;

    std::string filename = "sigma_b1.pdf";
    std::string xlabel   = "W    [GeV]";
    std::string ylabel   = "#sigma(#gamma p #rightarrow b_{1} + anything)   [#mub]";
    bool PRINT_TO_CMDLINE = true;

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    for (int n = 0; n < amps.size(); n++)
    {
        std::cout << "\nPrinting amplitude: " << amps[n]->_identifier << "\n";

        auto F = [&](double x)
        {
            return amps[n]->integrated_xsection(x*x) * 1.E-3;
        };

        std::array<std::vector<double>, 2> x_fx; 
        x_fx = vec_fill(N, F, xmin, xmax, PRINT_TO_CMDLINE);

        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->_identifier);
    }

    plotter->SetXaxis(xlabel, xmin, std::ceil(xmax));
    plotter->SetYaxis(ylabel, ymin, ymax);
    
    plotter->SetLegend(0.24, 0.28);
    plotter->SetLegendOffset(0.5, 0.08);

    // Output to file
    plotter->Plot(filename);

    return 1.;
};