
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
    // Zc(3900)+
    // ---------------------------------------------------------------------------

    auto   Zc  = new triple_regge(M_ZC3900, "Z_{c}(3900)^{+}");
    double gZc = 5.17E-2; // Gamma - Pi - Zc coupling

    auto beta_Zc = [&](double t)
    {
        return sqrt(2.) * sqrt(alphaPrime) * (gZc / M_ZC3900) * sqrt(Kallen(t, 0., M_ZC3900*M_ZC3900));
    };

    Zc->add_term(alphaPi, beta_Zc, &sigma_tot_pipp);
    Zc->set_formfactor(true, LamPi);

    // ---------------------------------------------------------------------------
    // Zb(10610)+
    // ---------------------------------------------------------------------------

    auto   Zb  = new triple_regge(M_ZB10610, "Z_{b}(10610)^{+}");
    double gZb = 5.8E-2; // Gamma - Pi - Zb coupling

    auto beta_Zb = [&](double t)
    {
        return sqrt(2.) * sqrt(alphaPrime) * (gZb / M_ZB10610) * sqrt(Kallen(t, 0., M_ZB10610*M_ZB10610));
    };

    Zb->add_term(alphaPi, beta_Zb, &sigma_tot_pipp);
    Zb->set_formfactor(true, LamPi);

    // ---------------------------------------------------------------------------
    // Zb(10650)+
    // ---------------------------------------------------------------------------

    auto   Zbp  = new triple_regge(M_ZB10650, "Z_{b}(10650)^{+}");
    double gZbp = 2.9E-2; // Gamma - Pi - Zb' coupling

    auto beta_Zbp = [&](double t)
    {
        return sqrt(2.) * sqrt(alphaPrime) * (gZbp / M_ZB10650) * sqrt(Kallen(t, 0., M_ZB10650*M_ZB10650));
    };

    Zbp->add_term(alphaPi, beta_Zbp, &sigma_tot_pipp);
    Zbp->set_formfactor(true, LamPi);

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<triple_regge*> amps;
    amps.push_back(Zc);
    amps.push_back(Zb);
    amps.push_back(Zbp);

    int N = 100;

    // Plotter object
    auto plotter = new jpacGraph1D();


    // ---------------------------------------------------------------------------
    // Integrate in both t and M2
    // ---------------------------------------------------------------------------
    double xmin = 5.;
    double xmax = 100.;

    double ymax = 2.E0;
    double ymin = 3.E-3;

    std::string filename = "sigma_Z.pdf";
    std::string xlabel   = "W   [GeV]";
    std::string ylabel   = "#sigma(#gamma p #rightarrow Z + anything)   [nb]";

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    for (int n = 0; n < amps.size(); n++)
    {
        std::cout << "\nPrinting amplitude: " << amps[n]->_identifier << "\n";

        auto F = [&](double W)
        {
            return amps[n]->integrated_xsection(W*W);
        };

        std::array<std::vector<double>, 2> x_fx; 
        x_fx = vec_fill(N, F, xmin, xmax, true);
        
        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->_identifier);
    }

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    // plotter->SetYlogscale(true);
    
    plotter->SetLegend(0.24, 0.18);
    plotter->SetLegendOffset(0.5, 0.14);

    // Output to file
    plotter->Plot(filename);

    return 1.;
};