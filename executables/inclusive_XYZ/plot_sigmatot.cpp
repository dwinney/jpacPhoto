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
    std::array<double, 2> bounds = {xmin, xmax};
    
    std::string filename = "sigma.pdf";
    std::string xlabel   = "W  [GeV]";
    std::string ylabel   = "#sigma_{tot}^{#pip}   [mb] ";

    // PDG parameterizations
    std::unique_ptr<total_xsection> PDG_pimp(  new PDG_parameterization(M_PION, M_PROTON, {-1., 1., 9.56, 1.767, 18.75}));
    std::unique_ptr<total_xsection> PDG_pipp(  new PDG_parameterization(M_PION, M_PROTON, {+1., 1., 9.56, 1.767, 18.75}));

    // JPAC parameterizations
    std::unique_ptr<total_xsection> JPAC_pimp( new JPAC_parameterization(-1, true) );
    std::unique_ptr<total_xsection> JPAC_pipp( new JPAC_parameterization(+1, true) );

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    std::unique_ptr<jpacGraph1D> plotter( new jpacGraph1D() );

    // ---------------------------------------------------------------------------
    // Print the phase-space for each kinematics

    double piPlus, useJPAC;
    auto F = [&](double W)
    {
        double s = W*W;
        return piPlus * ( useJPAC * JPAC_pipp->eval(s) + !useJPAC * PDG_pipp->eval(s) ) 
            + !piPlus * ( useJPAC * JPAC_pimp->eval(s) + !useJPAC * PDG_pimp->eval(s) );
    };

    piPlus = true; useJPAC = true;
    plotter->AddEntry(N, F, {xmin, xmax}, "#pi^{+} p");

    piPlus = true; useJPAC = false;
    plotter->AddDashedEntry(N, F, {xmin, xmax});

    piPlus = false; useJPAC = true;
    plotter->AddEntry(N, F, {xmin, xmax}, "#pi^{-} p");

    piPlus = false; useJPAC = false;
    plotter->AddDashedEntry(N, F, {xmin, xmax});

    // Axes and legend options
    plotter->SetXaxis(xlabel, 1., xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(true);

    plotter->SetLegend(0.7, 0.2);
    plotter->SetLegendOffset(0.3, 0.075);

    // Output to file
    plotter->Plot(filename);

    return 0;
};