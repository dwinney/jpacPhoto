#include "constants.hpp"
#include "reaction_kinematics.hpp"
#include "amplitudes/dirac_exchange.hpp"
#include "amplitudes/vector_exchange.hpp"

#include "jpacGraph1D.hpp"
#include "jpacUtils.hpp"

#include <cstring>

using namespace jpacPhoto;

int main( int argc, char** argv )
{
    // Set up kinematics of the Pcs K final state
    double M_PCS = 4.5880;
    reaction_kinematics kPcsK (M_KAON, M_PCS);
    kPcsK.set_meson_JP(0, -1);
    kPcsK.set_baryon_JP(1, -1);

    // Couplings 
    double gPhoton  =  0.0412501;
    double gNucleon = -6.4469; 
    double cutoff   = 0.9;

    // Create the amplitude
    dirac_exchange LamEx(&kPcsK, M_LAMBDA, "#Lambda exchange");
    LamEx.set_params({gPhoton, gNucleon});
    LamEx.set_formfactor(3, cutoff);

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
        double t = kPcsK.t_man(s, theta);
        switch (o)
        {
            case 1: return LamEx.A_LL(s, t);
            case 2: return LamEx.K_LL(s, t);
            case 3: return LamEx.beam_asymmetry_4pi(s, t);
            default: return 0.;
        };
        return 0.;
    };

    // ---------------------------------------------------------------------------
    // double  xmin = 5.;
    // double  xmax = 12.;

    // double  ymin = 0.;
    // double  ymax = 0.;

    // std::string ylabel    = "#sigma(#gamma p #rightarrow P K)  [nb]";
    // std::string xlabel    = "W_{#gammap}  [GeV]";

    // auto F = [&](double x)
    // {
    //     return LamEx.integrated_xsection(x*x);
    // };

    o = 1;
    plotter->AddEntry(N, F, {xmin,xmax}, "A_{LL}");
    o = 2;
    plotter->AddEntry(N, F, {xmin,xmax}, "K_{LL}");
    o = 3;
    plotter->AddEntry(N, F, {xmin,xmax}, "#Sigma_{4#pi}");

    plotter->SetXaxis(xlabel, xmin, xmax);
    plotter->SetYaxis(ylabel,  ymin, ymax);
    plotter->SetLegendOffset(0.3, 0.1);
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