
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

    // ---------------------------------------------------------------------------
    // Zc(3900)+
    // ---------------------------------------------------------------------------

    auto Zc  = new triple_regge(M_ZC3900, "Z_{c}(3900)");
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
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<triple_regge*> amps;
    amps.push_back(Zc);

    int N = 100;

    double  ymin = 0.;
    double  ymax = 300.;

    std::string filename = "intZ.pdf";
    std::string xlabel   = "-t      [GeV^{2}]";
    std::string ylabel   = "M^{2}   [GeV^{2}] ";

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    auto plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    double ws[3] = {21., 20., 19.};

    for (int i = 0; i < 3; i++)
    {
        double W = ws[i];
        double s = W*W;

        double xmin = -Zc->_kinematics->t_bounds(+1, s);   
        double xmax = -Zc->_kinematics->t_bounds(-1, s);

        auto F = [&](double t)
        {
            return Zc->_kinematics->M2_bounds(-1, s, t);
        };

        std::array<std::vector<double>, 2> x_fx; 
        x_fx = vec_fill(N, F, xmin, xmax, true);

        auto sF = [&](double t)
        {
            return Zc->_kinematics->M2_bounds(+1, s, -t);
        };

        std::array<std::vector<double>, 2> x_fx2; 
        x_fx2 = vec_fill(N, sF, xmin, xmax, true);

        x_fx[0].insert( x_fx[0].end(), x_fx2[0].begin(), x_fx2[0].end());
        x_fx[1].insert( x_fx[1].end(), x_fx2[1].begin(), x_fx2[1].end());

        std::ostringstream streamObj;
        streamObj << std::setprecision(4) << "W = " << W << " GeV";
        std::string label = streamObj.str();
        plotter->AddEntry(x_fx[0], x_fx[1], label);
    }

    double xMin = -Zc->_kinematics->t_bounds(+1, ws[2]*ws[2]);   
    double xMax = -Zc->_kinematics->t_bounds(-1, ws[2]*ws[2]);
    plotter->SetXaxis(xlabel, xMin, xMax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    

    plotter->SetLegend(0.2, 0.3);
    plotter->SetLegendOffset(0.5, 0.12);

    // Output to file
    plotter->Plot(filename);

    return 1.;
};