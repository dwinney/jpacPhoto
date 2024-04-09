// Model for high energy photoproduction of the pi Delta system via natural
// and unnatural exchanges.
//
// Adapted from models originally written by J. Nys and V. Mathieu in [1]
// are fit to new GlueX data
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Univeritat Bonn (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// References:
// [1] - https://arxiv.org/abs/1710.09394
// ------------------------------------------------------------------------------


#include "kinematics.hpp"
#include "amplitude.hpp"
#include "regge/natural_exchange.hpp"
#include "regge/unnatural_exchange.hpp"
#include "fitter.hpp"
#include "piDelta/plots.hpp"
#include "piDelta/data.hpp"

#include "plotter.hpp"

using namespace jpacPhoto;

struct piDelta_analysis
{
    static std::string data_type(int i)
    {
        switch (i)
        {
            case 0:  return "dsig/dt";
            case 1:  return "rho_H (0, 1,  1)";
            case 2:  return "rho_H (0, 3,  1)";
            case 3:  return "rho_H (0, 3, -1)";
            case 4:  return "rho_H (1, 1,  1)";
            case 5:  return "rho_H (1, 3,  3)";
            case 6:  return "rho_H (1, 3,  1)";
            case 7:  return "rho_H (1, 3, -1)";
            case 8:  return "rho_H (2, 3,  1)";
            case 9:  return "rho_H (2, 3, -1)";
            case 10: return "Sigma_4pi";
            default: return "ERROR";
        };
    };  

    static double fcn(const std::vector<data_set> &sets, amplitude amp)
    {
        double chi2 = 0;
        for (auto datum : sets)
        {
            if      (datum._type == 0)  chi2 += chi2_dif(datum, amp); 
            else if (datum._type == 10) chi2 += chi2_BA( datum, amp);
            else                        chi2 += chi2_SDME(datum, amp);
        };
        return chi2;
    };

    static double chi2_dif(data_set &data, amplitude amp)
    {
        double chi2 = 0;
        for (int i = 0; i < data._N; i++)
        {
            double s = s_cm(data._extras[0]);
            double t = - data._x[i];

            double sig_th = amp->differential_xsection(s, t) * 1E-3;
            double sig_ex = data._z[i];
            double error  = data._zerr[0][i];
            chi2 += pow((sig_th - sig_ex) / error, 2);
        };
        return chi2;
    };

    static double chi2_BA(data_set &data, amplitude amp)
    {
        double chi2 = 0;
        for (int i = 0; i < data._N; i++)
        {
            double s = s_cm(data._extras[0]);
            double t = - data._x[i];

            double sig_th = amp->beam_asymmetry_4pi(s, t);
            double sig_ex = data._z[i];
            double error  = data._zerr[0][i];
            chi2 += pow((sig_th - sig_ex) / error, 2);
        };
        return chi2;
    };

    static double chi2_SDME(data_set &data, amplitude amp)
    {
        double chi2 = 0;
        std::array<int, 3> ids = piDelta::SDME_indices(data._type); 
        for (int i = 0; i < data._N; i++)
        {
            double s = s_cm(data._extras[0]);
            double t = - data._x[i];

            double sig_th = (ids[0] == 2) ? std::imag(amp->bSDME_H(ids[0], ids[1], ids[2], s, t))
                                          : std::real(amp->bSDME_H(ids[0], ids[1], ids[2], s, t));
            double sig_ex = data._z[i];
            double error  = data._zerr[0][i];
            chi2 += pow((sig_th - sig_ex) / error, 2);
        };
        return chi2;
    };
};

void fit()
{
    using namespace jpacPhoto;

    // ---------------------------------------------------------------------------
    // Set up amplitude
    // ---------------------------------------------------------------------------

    int kPole = regge::natural_exchange::kPole;
    int kCut  = regge::natural_exchange::kCut;

    //--------------------------------------------
    // Pole model 

    // Degenerate natural trajectory
    double aPN = 0.9, a0N  = 1 - aPN*M_RHO*M_RHO; 
    double aPU = 0.7, a0U  =   - aPU*M_PION*M_PION; 
    
    // Top couplings
    double ga2 = 0.204268,  grho = 0.109677;
    double gpi = 0.428507,  gb1  = 0.0976007;

    // Bottom couplings
    double g0 = 19.16;                        // gpiDelta
    double g1 = -48.2, g2 = -52.4, g3 = 40.2; // grhoDelta

    // Suppression factors
    double bN = -0.4161, bU = 0.06043;

    // Pion absorption constant
    double cpi = 1.0587;
    
    //-----------------------------------------
    // Cut model 

    double bUc =  0.135,    bNc =  2.122,  cpic =  1.041;
    double g1c = -370.777, g2c = -242.416, g3c  = -138.9583;
    
    //--------------------------------------------------------------------------
    // Kinematics for pi Delta

    kinematics kpi = new_kinematics(M_PION, M_DELTA);
    kpi->set_meson_JP( PSEUDOSCALAR );
    kpi->set_baryon_JP( THREEPLUS );

    //--------------------------------------------------------------------------
    // Natural exchanges

    amplitude rho = new_amplitude<regge::natural_exchange>(kpi, -1, "#rho");
    rho->set_parameters({a0N, aPN,  grho, g1, g2, g3, bN});
    
    amplitude a2  = new_amplitude<regge::natural_exchange>(kpi, +1, "a_{2}");
    a2->set_parameters( {a0N, aPN, -ga2,  g1, g2, g3, bN});

    amplitude N = rho + a2;
    N->set_id("Natural");

    //--------------------------------------------------------------------------
    // Unnatural exchanges

    amplitude pi  = new_amplitude<regge::unnatural_exchange>(kpi, +1, "#pi");
    pi->set_parameters( {a0U, aPU, -gpi,  g0, bU, cpi});

    amplitude b1  = new_amplitude<regge::unnatural_exchange>(kpi, -1, "b_{1}");
    b1->set_parameters( {a0U, aPU,  gb1,  g0, bU, cpi});

    amplitude U = pi + b1;
    U->set_id("Unnatural");

    //--------------------------------------------------------------------------
    // Total sum

    amplitude total = N + U;
    total->set_id("Sum");

    // ---------------------------------------------------------------------------
    // Set up fitter
    // ---------------------------------------------------------------------------


    std::vector<std::string> labels = {
        "a0 (rho)", "aP (rho)", "gT (rho)", "gB1 (rho)", "gB2 (rho)", "gB3 (rho)", "bN (rho)",
        "a0 (a2)",  "aP (a2)",  "gT (a2)",  "gB1 (a2)",  "gB2 (a2)",  "gB3 (a2)",  "bN (a2)",
        "a0 (pi)", "aP (pi)", "gT (pi)", "gB (pi)", "bU (pi)", "c (pi)",
        "a0 (b1)", "aP (b1)", "gT (b1)", "gB (b1)", "bU (b1)", "c (b1)"
    };
    std::vector<double> pars = {
        a0N, aPN,  grho, g1, g2, g3, bN, a0N, aPN, -ga2,  g1, g2, g3, bN,
        a0U, aPU, -gpi,  g0, bU, cpi, a0U, aPU,  gb1,  g0, bU, cpi
    };

    fitter<piDelta_analysis> fitter(total, "Midgrad", 1.E-4);
    fitter.set_parameter_labels(labels);
    for (int i = 0; i < labels.size(); i++) fitter.fix_parameter(labels[i], pars[i]);
    
    fitter.free_parameter("bU (b1)");
    fitter.free_parameter("bU (pi)");
    fitter.sync_parameter("bU (b1)", "bU (pi)");
    fitter.free_parameter("gB (b1)");
    fitter.free_parameter("gB (pi)");
    fitter.sync_parameter("gB (b1)", "gB (pi)");
    fitter.free_parameter("a0 (b1)");
    fitter.sync_parameter("c (b1)", "c (pi)");
    
    fitter.add_data( piDelta::differential() );
    fitter.add_data( piDelta::beam_asymmetry() );
    fitter.add_data( piDelta::SDMEs() );
    fitter.do_fit();

    // ---------------------------------------------------------------------------
    // Plot the results
    // ---------------------------------------------------------------------------

    plotter plotter;

    plot pdif = piDelta::plot_differential(plotter);
    pdif.add_curve({1E-3, 1.}, [&](double mt){ return total->differential_xsection(s_cm(8.0), -mt) * 1E-3; }, "Total");

    plot pBA = piDelta::plot_beam_asymmetry(plotter);
    pBA.add_curve({1E-3, 1.}, [&](double mt){ return total->beam_asymmetry_4pi(s_cm(8.0), -mt); }, "Total");
    for (auto amp : extract_subamplitudes(total))
    {
        pdif.add_curve({1E-3, 1.}, [&](double mt){ return amp->differential_xsection(s_cm(8.0), -mt) * 1E-3; }, amp->id());
        pBA.add_curve({1E-3, 1.}, [&](double mt){ return amp->beam_asymmetry_4pi(s_cm(8.0), -mt); }, amp->id());
    };

    std::vector<plot> pSDMEs = piDelta::plot_SDMEs(plotter);
    for (int i = 1; i <= 9; i++)
    {
        std::array<int,3> ind = piDelta::SDME_indices(i);
        if (ind[0] == 2) pSDMEs[i-1].add_curve({1E-3, 1.}, [&](double mt){ return std::imag(total->bSDME_H(ind[0], ind[1], ind[2], s_cm(8.5), -mt)); }, "Pole model");
        else             pSDMEs[i-1].add_curve({1E-3, 1.}, [&](double mt){ return std::real(total->bSDME_H(ind[0], ind[1], ind[2], s_cm(8.5), -mt)); }, "Pole model");
    };

    plotter.combine({2,1}, {pdif, pBA}, "diff_and_BA.pdf");
    plotter.combine({3,3},  pSDMEs,     "SDMEs.pdf");
};