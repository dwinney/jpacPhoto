#include "inclusive/total_xsection.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>
#include <memory>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------

    int N = 500;

    double  ymin = 7.;
    double  ymax = 4.E2;

    double xmin = (M_PION + M_PROTON) + 0.01;   
    double xmax = 3;

    std::string filename = "sigma.pdf";
    std::string xlabel   = "W  [GeV]";
    std::string ylabel   = "#sigma_{tot}^{#pip}   [mb] ";

    std::unique_ptr<total_xsection> PDG_pimp(  new PDG_parameterization(M_PION, M_PROTON, {-1., 1., 9.56, 1.767, 18.75}));
    std::unique_ptr<total_xsection> PDG_pipp(  new PDG_parameterization(M_PION, M_PROTON, {+1., 1., 9.56, 1.767, 18.75}));

    std::unique_ptr<total_xsection> JPAC_pimp_NR( new JPAC_parameterization(-1, false) );
    std::unique_ptr<total_xsection> JPAC_pipp_NR( new JPAC_parameterization(+1, false) );

    std::unique_ptr<total_xsection> JPAC_pimp_R( new JPAC_parameterization(-1, true) );
    std::unique_ptr<total_xsection> JPAC_pipp_R( new JPAC_parameterization(+1, true) );

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    std::unique_ptr<jpacGraph1D> plotter( new jpacGraph1D() );

    // ---------------------------------------------------------------------------
    // Print the phase-space for each kinematics

    std::array<std::vector<double>, 2> x_fx; 

    auto F = [&](double W)
    {
        double s = W*W;
        return JPAC_pimp_R->eval(s);
    };
    x_fx = vec_fill(N, F, xmin, xmax);
    plotter->AddEntry(x_fx[0], x_fx[1], "#pi^{-} p");

    x_fx[0].clear(); x_fx[1].clear();
    auto G = [&](double W)
    {
        double s = W*W;
        return JPAC_pipp_R->eval(s);
    };
    x_fx = vec_fill(N, G, xmin, xmax);
    plotter->AddEntry(x_fx[0], x_fx[1], "#pi^{+} p");
    
    x_fx[0].clear(); x_fx[1].clear();
    auto FF = [&](double W)
    {
        double s = W*W;
        return PDG_pimp->eval(s);
    };
    x_fx = vec_fill(N, FF, xmin, xmax);
    plotter->AddDashedEntry(x_fx[0], x_fx[1]);

    x_fx[0].clear(); x_fx[1].clear();

    auto GG = [&](double W)
    {
        double s = W*W;
        return PDG_pipp->eval(s);
    };
    x_fx = vec_fill(N, GG, xmin, xmax);
    plotter->AddDashedEntry(x_fx[0], x_fx[1]);

    // Axes and legend options
    plotter->SetXaxis(xlabel, 1., xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(true);

    plotter->SetLegend(0.7, 0.2);
    plotter->SetLegendOffset(0.3, 0.075);

    // Output to file
    plotter->Plot(filename);

    return 1.;
};