
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

    auto Zc  = new triple_regge(M_ZC3900,  "Z_{c}(3900)^{+}");
    auto Zb  = new triple_regge(M_ZB10610, "Z_{b}(10610)^{+}");
    auto Zbp = new triple_regge(M_ZB10650, "Z_{b}(10650)^{+}");

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<triple_regge*> amps;
    amps.push_back(Zc);
    amps.push_back(Zb);
    amps.push_back(Zbp);

    int N = 100;

    double W = 20.;
    double s = W*W;

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
    // Print the phase-space for each kinematics

    for (int i = 0; i < amps.size(); i++)
    {
        double xmin = -amps[i]->_kinematics->t_bounds(+1, s);   
        double xmax = -amps[i]->_kinematics->t_bounds(-1, s);

        auto F = [&](double t)
        {
            return amps[i]->_kinematics->M2_bounds(-1, s, t);
        };

        std::array<std::vector<double>, 2> x_fx; 
        x_fx = vec_fill(N, F, xmin, xmax, true);

        auto sF = [&](double t)
        {
            return amps[i]->_kinematics->M2_bounds(+1, s, -t);
        };

        std::array<std::vector<double>, 2> x_fx2; 
        x_fx2 = vec_fill(N, sF, xmin, xmax, true);

        x_fx[0].insert( x_fx[0].end(), x_fx2[0].begin(), x_fx2[0].end());
        x_fx[1].insert( x_fx[1].end(), x_fx2[1].begin(), x_fx2[1].end());

        plotter->AddEntry(x_fx[0], x_fx[1], amps[i]->_identifier);
    }

    double xMin = -Zc->_kinematics->t_bounds(+1, s);   
    double xMax = -Zc->_kinematics->t_bounds(-1, s);
    plotter->SetXaxis(xlabel, xMin, xMax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    

    std::ostringstream streamObj;
    streamObj << std::setprecision(4) << "W = " << W << " GeV";
    std::string header = streamObj.str();
    plotter->SetLegend(0.6, 0.6, header);
    plotter->SetLegendOffset(0.5, 0.12);

    // Output to file
    plotter->Plot(filename);

    return 1.;
};