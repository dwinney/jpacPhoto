
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
    alphaPi->set_minimum_spin(0);`

    // ---------------------------------------------------------------------------
    // Zc(3900)+
    // ---------------------------------------------------------------------------

    auto Zc  = new triple_regge(M_ZC3900, "Z_{c}(3900)^{+}");
    double gZc = 5.17E-2; // Gamma - Pi - Zc coupling

    // We have two contributions from no-flip and flip contributions
    auto beta_nonflip_Zc = [&](double t)
    {
        return (sqrt(alphaPrime) / 2.) * (gZc / M_ZC3900) * (M2_PION - t);
    };
    auto beta_flip_Zc    = [&](double t)
    {
        return (sqrt(alphaPrime) / 4.) * (gZc / M_ZC3900) * (M2_PION - t) * (sqrt(-t) / M_ZC3900);
    };

    Zc->add_term(alphaPi, beta_nonflip_Zc, &sigma_tot_pipp);
    Zc->add_term(alphaPi, beta_flip_Zc,    &sigma_tot_pipp);

    // ---------------------------------------------------------------------------
    // Zb(10610)+
    // ---------------------------------------------------------------------------

    auto Zb  = new triple_regge(M_ZB10610, "Z_{b}(10610)^{+}");
    double gZb = 5.8E-2; // Gamma - Pi - Zb coupling

    // We have two contributions from no-flip and flip contributions
    auto beta_nonflip_Zb = [&](double t)
    {
        return (sqrt(alphaPrime) / 2.) * (gZb / M_ZB10610) * (M2_PION - t);
    };
    auto beta_flip_Zb    = [&](double t)
    {
        return (sqrt(alphaPrime) / 4.) * (gZb / M_ZB10610) * (M2_PION - t) * (sqrt(-t) / M_ZB10610);
    };

    Zb->add_term(alphaPi, beta_nonflip_Zb, &sigma_tot_pipp);
    Zb->add_term(alphaPi, beta_flip_Zb,    &sigma_tot_pipp);

    // ---------------------------------------------------------------------------
    // Zb(10650)+
    // ---------------------------------------------------------------------------

    auto Zbp  = new triple_regge(M_ZB10650, "Z_{b}(10650)^{+}");
    double gZbp = 2.9E-2; // Gamma - Pi - Zb' coupling

    // We have two contributions from no-flip and flip contributions
    auto beta_nonflip_Zbp = [&](double t)
    {
        return (sqrt(alphaPrime) / 2.) * (gZbp / M_ZB10650) * (M2_PION - t);
    };
    auto beta_flip_Zbp    = [&](double t)
    {
        return (sqrt(alphaPrime) / 4.) * (gZbp / M_ZB10650) * (M2_PION - t) * (sqrt(-t) / M_ZB10650);
    };

    Zbp->add_term(alphaPi, beta_nonflip_Zbp, &sigma_tot_pipp);
    Zbp->add_term(alphaPi, beta_flip_Zbp,    &sigma_tot_pipp);

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<triple_regge*> amps;
    amps.push_back(Zc);
    // amps.push_back(Zb);
    // amps.push_back(Zbp);

    int N = 100;

    double W = 40.;

    // ---------------------------------------------------------------------------
    // Integrate in M2
    // ---------------------------------------------------------------------------
    double  xmax = 30.;

    double  ymin = 2.E-2;
    double  ymax = 2.E10;

    std::string filename = "triple_Z.pdf";
    std::string xlabel   = "#it{-t} [GeV^{2}]";
    std::string ylabel   = "#frac{d#sigma}{dt}  [nb GeV^{-2}]";


    // Plotter object
    auto plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    for (int n = 0; n < amps.size(); n++)
    {
        std::cout << "\nPrinting amplitude: " << amps[n]->_identifier << "\n";
        double s = W*W;
        double xmin = -(amps[n]->_kinematics->t_bounds(-1, s) + 10. * EPS);
        

        auto F = [&](double mt)
        {
            return amps[n]->dsigma_dt(s, -mt);
        };

        std::array<std::vector<double>, 2> x_fx; 
        x_fx = vec_fill(N, F, xmin, xmax, true);

        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->_identifier);
    }

    plotter->SetXaxis(xlabel, 0., xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(true);
    
    std::ostringstream streamObj;
    streamObj << std::setprecision(4) << "W = " << W << " GeV";
    std::string header = streamObj.str();
    plotter->SetLegend(0.6, 0.55, header);
    plotter->SetLegendOffset(0.5, 0.18);

    // Output to file
    plotter->Plot(filename);

    return 1.;
};