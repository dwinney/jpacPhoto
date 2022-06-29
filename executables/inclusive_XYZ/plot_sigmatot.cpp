#include "inclusive/total_xsection_options.hpp"

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

    double xmin = pow((M_PROTON + M_PION) + 0.03, 2.);   
    double xmax = 6.;
    std::array<double, 2> bounds = {xmin, xmax};
    
    std::string filename = "sigma.pdf";
    std::string xlabel   = "M^{2} [GeV^{2}]";
    std::string ylabel   = "#sigma_{tot}^{#pip}   [mb] ";

    // PDG parameterizations
    total_xsection * PDG_pipp = get_total_xsection( PDG_pipp_onlyRegge );
    total_xsection * PDG_pimp = get_total_xsection( PDG_pimp_onlyRegge );

    // JPAC - SAID parameterizations
    total_xsection * JPAC_pimp = get_total_xsection( JPAC_pimp_withResonances );
    total_xsection * JPAC_pipp = get_total_xsection( JPAC_pipp_withResonances );

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    std::unique_ptr<jpacGraph1D> plotter( new jpacGraph1D() );

    // ---------------------------------------------------------------------------
    // Print the phase-space for each kinematics

    double piPlus, useJPAC;
    auto F = [&](double s)
    {
        return piPlus * ( useJPAC * JPAC_pipp->eval(s, M2_PION) + !useJPAC * PDG_pipp->eval(s, M2_PION) ) 
            + !piPlus * ( useJPAC * JPAC_pimp->eval(s, M2_PION) + !useJPAC * PDG_pimp->eval(s, M2_PION) );
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