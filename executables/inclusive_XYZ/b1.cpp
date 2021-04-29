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
    auto alphaPi = new linear_trajectory(+1, alpha0, alphaPrime, "#pi trajectory");
    alphaPi->set_minimum_spin(0);

    double LamPi = .9;  // 900 MeV cutoff for formfactor

    // ---------------------------------------------------------------------------
    // b1 meson
    // ---------------------------------------------------------------------------

    auto   b1  = new triple_regge(M_B1, "#Gamma(-#alpha(t))");
    double gb1 = 0.24; // Gamma - Pi - b1 coupling

    auto beta_b1 = [&](double t)
    {
        return sqrt(2.) * sqrt(alphaPrime) * (gb1 / M_B1) * sqrt(Kallen(t, 0., M2_B1));
    };

    b1->add_term(alphaPi, beta_b1, &sigma_tot_pipp);
    b1->set_formfactor(true, LamPi);

    auto   b11  = new triple_regge(M_B1, "1");
    b11->add_term(alphaPi, beta_b1, &sigma_tot_pipp);
    b11->set_formfactor(true, LamPi);
    b11->set_debug(1);

    auto   b12  = new triple_regge(M_B1, "2");
    b12->add_term(alphaPi, beta_b1, &sigma_tot_pipp);
    b12->set_formfactor(true, LamPi);
    b12->set_debug(2);

    auto   b13  = new triple_regge(M_B1, "#Gamma(-#alpha(t'))");
    b13->add_term(alphaPi, beta_b1, &sigma_tot_pipp);
    b13->set_formfactor(true, LamPi);
    b13->set_debug(3);

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<triple_regge*> amps;
    amps.push_back(b1);
    amps.push_back(b11);
    amps.push_back(b12);
    // amps.push_back(b13);

    // Plotter object
    auto plotter = new jpacGraph1D();

    int N = 100;
    double xmin = 0.6;
    double xmax = 0.9999;

    double ymin = 0.;
    double ymax = 10.;

    double W = 30.;
    double s = W*W;

    std::string filename = "sigma_b1.pdf";
    std::string xlabel   = "x    [GeV^{2}]";
    std::string ylabel   = "d#sigma / dx    [#mub]";
    bool PRINT_TO_CMDLINE = true;

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    for (int n = 0; n < amps.size(); n++)
    {
        std::cout << "\nPrinting amplitude: " << amps[n]->_identifier << "\n";

        auto F = [&](double x)
        {
            return amps[n]->dsigma_dx(s, x) * 1.E-3;
        };

        std::array<std::vector<double>, 2> x_fx; 
        x_fx = vec_fill(N, F, xmin, xmax, PRINT_TO_CMDLINE);

        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->_identifier);
    }

    plotter->SetXaxis(xlabel, xmin, std::ceil(xmax));
    plotter->SetYaxis(ylabel, ymin, ymax);
    
    plotter->SetLegend(0.54, 0.28);
    plotter->SetLegendOffset(0.5, 0.08);

    // Output to file
    plotter->Plot(filename);

    return 1.;
};