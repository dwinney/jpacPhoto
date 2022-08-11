#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "dirac_exchange.hpp"
#include "vector_exchange.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

void compare_PcsK()
{
    // Set up kinematics of the Pcs K final state
    double M_PCS = 4.5880;

    // Couplings 
    double g1  = 4.45E-3;
    double g3  = 3.01E-3;
    double gNucleon = -6.4469; 
    double cutoff   = 0.9;

    //Spin-1/2
    reaction_kinematics kPcsK1 (M_KAON, M_PCS);
    kPcsK1.set_meson_JP (0, -1);
    kPcsK1.set_baryon_JP(1, -1);

    dirac_exchange LamEx1(&kPcsK1, M_LAMBDA, "#Lambda exchange");
    LamEx1.set_params({g1, 0., gNucleon});
    LamEx1.set_formfactor(3, cutoff);

    //Spin-3/2
    reaction_kinematics kPcsK3 (M_KAON, M_PCS);
    kPcsK3.set_meson_JP (0, -1);
    kPcsK3.set_baryon_JP(3, -1);

    dirac_exchange LamEx3(&kPcsK3, M_LAMBDA, "#Lambda exchange");
    LamEx3.set_params({g3, 0., gNucleon});
    LamEx3.set_formfactor(3, cutoff);

    // // ---------------------------------------------------------------------------

    int   N = 50;
    bool PRINT_TO_COMMANDLINE  = true;
    std::string filename  = "PcsK.pdf";

    // ---------------------------------------------------------------------------
    // You shouldnt need to change anything below this line
    // ---------------------------------------------------------------------------

    // Plotter object
    jpacGraph1D * plotter = new jpacGraph1D();

    std::string ylabel    = "";
    std::string xlabel    = "#theta_{s}  [deg.]";
    // ---------------------------------------------------------------------------
    // Print the desired observable for each amplitude

    double xmin = 0.;
    double xmax = 180;
    double ymin = -1.;
    double ymax = +1.;
    double w = 5.5; 
    double s = w*w;
   
    int o;
    auto F = [&](double theta)
    {
        theta *= DEG2RAD;
        double t = kPcsK1.t_man(s, theta);
        switch (o)
        {
            case 1: return LamEx1.A_LL(s, t);
            case 2: return LamEx1.K_LL(s, t);
            case 3: return LamEx1.beam_asymmetry_4pi(s, t);
            default: return 0.;
        };
        return 0.;
    };

    auto G = [&](double theta)
    {
        theta *= DEG2RAD;
        double t = kPcsK3.t_man(s, theta);
        switch (o)
        {
            case 1: return LamEx3.A_LL(s, t);
            case 2: return LamEx3.K_LL(s, t);
            case 3: return LamEx3.beam_asymmetry_4pi(s, t);
            default: return 0.;
        };
        return 0.;
    };

    o = 1;
    plotter->AddEntry(N, F, {xmin,xmax}, "A_{LL}");
    plotter->AddDashedEntry(N, G, {xmin,xmax});

    o = 2;
    plotter->AddEntry(N, F, {xmin,xmax}, "K_{LL}");
    plotter->AddDashedEntry(N, G, {xmin,xmax});

    o = 3;
    plotter->AddEntry(N, F, {xmin,xmax}, "#Sigma_{4#pi}");
    plotter->AddDashedEntry(N, G, {xmin,xmax});


    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel,  ymin, ymax);
    plotter->SetLegendOffset(0.3, 0.15);
    plotter->RemoveLogo();

    // Add a header to legend to specify the fixed Q2
    std::ostringstream streamObj;
    streamObj << std::setprecision(4) << "W = " << w << " GeV";
    std::string header = streamObj.str();
    plotter->SetLegend(0.7, 0.25, header);

    // Output to file
    plotter->Plot(filename);

    delete plotter;
};