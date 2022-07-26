#include "inclusive_kinematics.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

void chew_low()
{
    // Kinematics object
    inclusive_kinematics kb1 (M_B1);
    kb1._s = 9.;

    // ---------------------------------------------------------------------------
    // Plotting options
    // ---------------------------------------------------------------------------
    int N = 200;

    double  ymin = 1.0;
    double  ymax = 3.3;

    double Wmin = M_PROTON + M_PION;
    double xmin = - kb1.TMINfromM2(Wmin*Wmin);   
    double xmax = - kb1.TMAXfromM2(Wmin*Wmin);  
    std::array<double, 2> bounds = {xmin, xmax};
    
    std::string filename = "phase-space.pdf";
    std::string xlabel   = "-#it{t} [GeV^{2}]";
    std::string ylabel   = "#it{M}^{2}  [GeV^{2}]";
    bool PRINT = true;

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();

    // ---------------------------------------------------------------------------
    // Print the phase-space for each kinematics

    auto M2max = [&](double mt)
    {
        return kb1.M2MAXfromT(-mt);
    };
    auto M2min = [&](double mt)
    {
        return kb1.M2MINfromT(-mt);
    };
    plotter->AddEntry(N, M2max, bounds, "", PRINT);
    plotter->AddEntry(N, M2min, bounds, "", PRINT);

    // Axes and legend options
    plotter->SetXaxis(xlabel, 0., xmax);
    plotter->SetYaxis(ylabel, ymin, ymax);

    plotter->SetLegend(false);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
};