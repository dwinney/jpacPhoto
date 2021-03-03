
#include "constants.hpp"
#include "regge_trajectory.hpp"
#include "inclusive/triple_regge.hpp"
#include "inclusive/field_fox_couplings.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>
#include <memory>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // trajectories
    auto alphaPom = new linear_trajectory(+1,  1., 0.37, "Pomeron");
    auto alphaReg = new linear_trajectory(+1, 0.5,   1., "Reggeon");

    // ---------------------------------------------------------------------------
    // X(3872)
    // ---------------------------------------------------------------------------

    auto X  = new triple_regge(M_X3872, "X(3872)");

    double gX = sqrt(3.6*3.6 + 8.2*8.2) * 1.E-3;
    double betaRpp = sqrt(56.08 / 0.389 / pow(0.9, 1. - 0.4625));

    // Reggeon - Reggeon - Pomeron coupling
    auto gX_RRP = [&](double t)
    {
        double betaX = gX*gX * (1./8. - t / (4. * M2_X3872));
        return betaX / (betaRpp * betaRpp) * G_RRP(t);
    };

    // Reggeon - Reggeon - Reggeon coupling
    auto gX_RRR = [&](double t)
    {
        double betaX = gX*gX * (1./8. - t / (4. * M2_X3872));
        return betaX / (betaRpp * betaRpp) * G_RRR(t);
    };

    X->add_term({alphaReg, alphaReg, alphaPom}, gX_RRP);
    X->add_term({alphaReg, alphaReg, alphaReg}, gX_RRR);

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    // which amps to plot
    std::vector<triple_regge*> amps;
    amps.push_back(X);

    int N = 100;

    double W = 50.;
    double x = 0.9;

    double  xmin = 0.;
    double  xmax = 1.;

    double  ymin = 1.E-3;
    double  ymax = 1.E1;

    std::string filename = "triple_X.pdf";
    std::string xlabel   = "#it{-t} [GeV^{2}]";
    std::string ylabel   = "E d#sigma / dt dM^{2}  [nb GeV^{-2}]";

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    auto plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude
    for (int n = 0; n < amps.size(); n++)
    {
        auto F = [&](double mt)
        {
            double s = W*W;
            double M2 = s * (1. - x);
            return amps[n]->invariant_xsection(s, -mt, M2);
        };

        std::array<std::vector<double>, 2> x_fx; 
        x_fx = vec_fill(N, F, xmin, xmax, true);

        plotter->AddEntry(x_fx[0], x_fx[1], amps[n]->_identifier);
    }

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(true);
    
    std::ostringstream streamObj;
    streamObj << std::setprecision(4) << "x = " << x << ",  W = " << W << " GeV";
    std::string header = streamObj.str();
    plotter->SetLegend(0.6, 0.6, header);

    // Output to file
    plotter->Plot(filename);

    return 1.;
};