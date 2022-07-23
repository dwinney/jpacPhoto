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

    double  ymin = 6.;
    double  ymax = 4.E2;

    double xmin = pow((M_PROTON + M_PION) + 0.03, 2.);   
    double xmax = 6.;
    std::array<double, 2> bounds = {xmin, xmax};
    
    std::string filename = "sigma.pdf";
    std::string xlabel   = "#it{M}^{2} [GeV^{2}]";
    std::string ylabel   = "#sigma_{tot}^{#pi^{*+} p}   [mb] ";

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

    double q2 = M2_PION;
    
    total_xsection * sigma;
    auto F = [&](double s)
    {
        return sigma->eval(s, q2);
    };

    sigma = JPAC_pipp;
    plotter->AddEntry(N, F, {xmin, xmax}, "#pi^{#plus} p", 1);

    sigma = PDG_pipp;
    plotter->AddDashedEntry(N, F, {xmin, xmax});

    sigma = JPAC_pimp;
    plotter->AddEntry(N, F, {xmin, xmax}, "#pi^{#minus} p");

    sigma = PDG_pimp;
    plotter->AddDashedEntry(N, F, {xmin, xmax});

    // Axes and legend options
    plotter->SetXaxis(xlabel, 1., xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(true);

    plotter->SetLegend(0.7, 0.2);
    plotter->SetLegendOffset(0.3, 0.075);

    // Axes and legend options
    plotter->SetXaxis(xlabel, 1., xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);
    plotter->SetYlogscale(true);

    plotter->SetLegend(0.7, 0.2);
    plotter->SetLegendOffset(0.3, 0.1);

    // Output to file
    plotter->Plot(filename);

    return 0;
};