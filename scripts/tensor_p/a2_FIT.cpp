// Predictions for a2(1320) photoproduction with both vector and axialvector exchanges
// Reproduces fig. 3 in [1] with couplings fit to data
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Universitat Bonn (HISKP)
// Email:        daniel.winney@iu.alumni.edu
//               winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------
// References 
// [1] - arXiv:2005.01617 [hep-ph]
// ------------------------------------------------------------------------------


#include "plotter.hpp"
#include "constants.hpp"
#include "regge/reggeon_exchange.hpp"

void a2_FIT()
{
    using namespace jpacPhoto;
    using namespace regge;

    // ------------------------------------------------------------------------------
    // Data (theres so few points that i just hard code em)

    std::vector<double> mt5  = {0.3, 0.55, 0.9, 1.3, 1.75};
    std::vector<double> mte5 = {0.1, 0.15, 0.2, 0.2, 0.25};
    std::vector<double> sig5 = {0.2375, 0.15871, 0.24709, 0.13775, 0.0702 };
    std::vector<double> err5 = { 0.28392 - 0.2375,  0.18564 - 0.15871, 0.26844 - 0.24709, 0.15167 - 0.13775, 0.07940 - 0.0702 };
    
    std::vector<double> mt4  = {0.3, 0.55, 0.9, 1.3};
    std::vector<double> mte4 = {0.1, 0.15, 0.2, 0.2};
    std::vector<double> sig4 = { 0.4306, 0.27198, 0.31951, 0.23616 };
    std::vector<double> err4 = { 0.4798 - 0.4306, 0.30633 - 0.27198, 0.34643 - 0.31951, 0.26123 - 0.23616 };

    // ------------------------------------------------------------------------------
    /// Couplings and inputs

    // Trajectories
    double alpha0_V = 0.5,                alphaP_V = 0.9;
    double alpha0_A = -0.7*M_PION*M_PION, alphaP_A = 0.7;

    // Top couplings
    std::array<double,2> beta_rho = {0.251/3, 1.060/3}, beta_omega = {0.251,  1.060};
    std::array<double,2> beta_h1  = {0.821, 0.581};
    
    std::array<double,5> gT_MIN = {0,     1/2., -1/sqrt(6), 0.,   0};
    std::array<double,5> gT_TMD = {-1./2, 1.0,  -1/sqrt(6), 0.,   0};
    std::array<double,5> gT_M1  = {0.,    1/4., -1/sqrt(6), 1/4., 0};

    // Bottom couplings
    std::array<double,2> gB_rho   = {1.63, 13.01};
    std::array<double,2> gB_omega = {8.13, 1.860};
    std::array<double,2> gB_h1    = {0.0,  25.24};

    // ------------------------------------------------------------------------------
    /// Amplitudes

    kinematics kA2 = new_kinematics(M_A2);
    kA2->set_meson_JP(TENSOR);

    amplitude rho_MIN = new_amplitude<reggeon_exchange>(kA2, +1, -1, "#rho");
    rho_MIN->set_option(reggeon_exchange::kRemoveZero);
    rho_MIN->set_option(reggeon_exchange::k2020_Minimal);
    rho_MIN->set_parameters( {alpha0_V, alphaP_V, 0,
                          beta_rho[0]*gT_MIN[0], beta_rho[0]*gT_MIN[1], beta_rho[0]*gT_MIN[2], beta_rho[0]*gT_MIN[3], beta_rho[0]*gT_MIN[4],
                          gB_rho[0], gB_rho[1]} );

    amplitude omega_MIN = new_amplitude<reggeon_exchange>(kA2, +1, -1, "#omega");
    omega_MIN->set_option(reggeon_exchange::k2020_Minimal);
    omega_MIN->set_parameters( {alpha0_V, alphaP_V, 0,
                            beta_omega[0]*gT_MIN[0], beta_omega[0]*gT_MIN[1], beta_omega[0]*gT_MIN[2], beta_omega[0]*gT_MIN[3], beta_omega[0]*gT_MIN[4], 
                            gB_omega[0], gB_omega[1]} );

    amplitude h1_MIN = new_amplitude<reggeon_exchange>(kA2, -1, -1, "#it{h}_{1}");
    h1_MIN->set_option(reggeon_exchange::k2020_Minimal);
    h1_MIN->set_parameters( {alpha0_A, alphaP_A, 0,
                            beta_h1[0]*gT_M1[0], beta_h1[0]*gT_M1[1], beta_h1[0]*gT_M1[2], beta_h1[0]*gT_M1[3], beta_h1[0]*gT_M1[4], 
                            gB_h1[0], gB_h1[1]} );

    amplitude rho_TMD = new_amplitude<reggeon_exchange>(kA2, +1, -1, "#rho");
    rho_TMD->set_option(reggeon_exchange::kRemoveZero);
    rho_TMD->set_option(reggeon_exchange::k2020_TMD);
    rho_TMD->set_parameters( {alpha0_V, alphaP_V, 0,
                          beta_rho[1]*gT_TMD[0], beta_rho[1]*gT_TMD[1], beta_rho[1]*gT_TMD[2], beta_rho[1]*gT_TMD[3], beta_rho[1]*gT_TMD[4],
                          gB_rho[0], gB_rho[1]} );

    amplitude omega_TMD = new_amplitude<reggeon_exchange>(kA2, +1, -1, "#omega");
    omega_TMD->set_option(reggeon_exchange::k2020_TMD);
    omega_TMD->set_parameters( {alpha0_V, alphaP_V, 0,
                            beta_omega[1]*gT_TMD[0], beta_omega[1]*gT_TMD[1], beta_omega[1]*gT_TMD[2], beta_omega[1]*gT_TMD[3], beta_omega[1]*gT_TMD[4], 
                            gB_omega[0], gB_omega[1]} );

    amplitude h1_TMD = new_amplitude<reggeon_exchange>(kA2, -1, -1, "#it{h}_{1}");
    h1_TMD->set_option(reggeon_exchange::k2020_Minimal);
    h1_TMD->set_parameters( {alpha0_A, alphaP_A, 0,
                            beta_h1[1]*gT_M1[0], beta_h1[1]*gT_M1[1], beta_h1[1]*gT_M1[2], beta_h1[1]*gT_M1[3], beta_h1[1]*gT_M1[4], 
                            gB_h1[0], gB_h1[1]} );

    amplitude MIN = rho_MIN + omega_MIN + h1_MIN;
    amplitude TMD = rho_TMD + omega_TMD + h1_TMD;

    // ------------------------------------------------------------------------------
    /// Plot results

    plotter plotter;
    double s;

    plot p1 = plotter.new_plot();
    p1.set_labels("#minus#it{t}  [GeV^{2}]", "d#sigma/d#it{t}  [#mub/GeV^{2}]");
    p1.set_ranges({0, 2.0}, {0.0, 0.6});
    p1.add_header("Minimal");
    p1.set_legend(0.65, 0.65);
    p1.set_legend_spacing(0.03);
    p1.add_data({mt4, sig4}, {mte4, err4}, jpacColor::Blue);
    p1.add_data({mt5, sig5}, {mte5, err5}, jpacColor::Red);
    s = s_cm(4);
    p1.add_curve( {0.07, 2}, [&](double mt){ return MIN->differential_xsection(s, -mt)    * 1E-3; }, "#it{E}_{#gamma} = 4 GeV");
    p1.add_dashed({0.07, 2}, [&](double mt){ return h1_MIN->differential_xsection(s, -mt) * 1E-3; });
    s = s_cm(5);
    p1.add_curve( {0.07, 2}, [&](double mt){ return MIN->differential_xsection(s, -mt)    * 1E-3; }, "#it{E}_{#gamma} = 5 GeV");
    p1.add_dashed({0.07, 2}, [&](double mt){ return h1_MIN->differential_xsection(s, -mt) * 1E-3; });

    plot p2 = plotter.new_plot();
    p2.set_labels("#minus#it{t}  [GeV^{2}]", "d#sigma/d#it{t}  [#mub/GeV^{2}]");
    p2.set_ranges({0, 2.0}, {0, 0.6});
    p2.add_header("TMD");
    p2.set_legend(0.65, 0.65);
    p2.add_data({mt4, sig4}, {mte4, err4}, jpacColor::Blue);
    p2.add_data({mt5, sig5}, {mte5, err5}, jpacColor::Red);
    s = s_cm(4);
    p2.add_curve( {0.07, 2}, [&](double mt){ return TMD->differential_xsection(s, -mt)    * 1E-3; });
    p2.add_dashed({0.07, 2}, [&](double mt){ return h1_TMD->differential_xsection(s, -mt) * 1E-3; });
    s = s_cm(5);
    p2.add_curve( {0.07, 2}, [&](double mt){ return TMD->differential_xsection(s, -mt)    * 1E-3; });
    p2.add_dashed({0.07, 2}, [&](double mt){ return h1_TMD->differential_xsection(s, -mt) * 1E-3; });

    plotter.combine({2,1}, {p1,p2}, "fig3.pdf");
};  