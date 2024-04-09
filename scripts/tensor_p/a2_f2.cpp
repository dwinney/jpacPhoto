// Predictions for a2 and f2 photoproduction with both vector and axialvector exchanges
// Reproduces fig. 4 in [1] with couplings fit to data
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

void a2_f2()
{
    using namespace jpacPhoto;
    using namespace regge;

    // ------------------------------------------------------------------------------
    // Data (theres so few points that i just hard code em)

    // A2 data
    std::vector<double> mt5  = {0.3, 0.55, 0.9, 1.3, 1.75};
    std::vector<double> mte5 = {0.1, 0.15, 0.2, 0.2, 0.25};
    std::vector<double> sig5 = {0.2375, 0.15871, 0.24709, 0.13775, 0.0702 };
    std::vector<double> err5 = { 0.28392 - 0.2375,  0.18564 - 0.15871, 0.26844 - 0.24709, 0.15167 - 0.13775, 0.07940 - 0.0702 };
    
    std::vector<double> mt4  = {0.3, 0.55, 0.9, 1.3};
    std::vector<double> mte4 = {0.1, 0.15, 0.2, 0.2};
    std::vector<double> sig4 = { 0.4306, 0.27198, 0.31951, 0.23616 };
    std::vector<double> err4 = { 0.4798 - 0.4306, 0.30633 - 0.27198, 0.34643 - 0.31951, 0.26123 - 0.23616 };

    // F2 data
    std::vector<double> mt3  = {0.45,     0.55,     0.65,     0.75,     0.85,     0.95};
    std::vector<double> mte3 = {0.05,     0.05,     0.05,     0.05,     0.05,     0.05};
    std::vector<double> sig3 = {0.307332, 0.283866, 0.244473, 0.266191, 0.245707, 0.224500};
    std::vector<double> err3 = {0.122933, 0.113546, 0.097789, 0.106476, 0.098283, 0.089800};

    // ------------------------------------------------------------------------------
    /// Couplings and inputs

    // Trajectories
    double alpha0_V = 0.5,                alphaP_V = 0.9;
    double alpha0_A = -0.7*M_PION*M_PION, alphaP_A = 0.7;

    // Top couplings
    std::array<double,2> abeta_rho = {0.251/3, 1.060/3}, abeta_omega = {0.251,  1.060};
    std::array<double,2> abeta_h1  = {0.821, 0.581};

    std::array<double,2> fbeta_rho = {0.251, 1.060}, fbeta_omega = {0.251/3,  1.060/3};
    std::array<double,2> fbeta_h1  = {0.821/3, 0.581/3};
    
    std::array<double,5> gT_MIN = {0,     1/2., -1/sqrt(6), 0.,   0};
    std::array<double,5> gT_TMD = {-1./2, 1.0,  -1/sqrt(6), 0.,   0};
    std::array<double,5> gT_M1  = {0.,    1/4., -1/sqrt(6), 1/4., 0};

    // Bottom couplings
    std::array<double,2> gB_rho   = {1.63, 13.01};
    std::array<double,2> gB_omega = {8.13, 1.860};
    std::array<double,2> gB_h1    = {0.0,  25.24};

    // ------------------------------------------------------------------------------
    /// a2 Amplitudes

    kinematics kA2 = new_kinematics(M_A2);
    kA2->set_meson_JP(TENSOR);

    amplitude a2rho_MIN = new_amplitude<reggeon_exchange>(kA2, +1, -1, "#rho");
    a2rho_MIN->set_option(reggeon_exchange::kRemoveZero);
    a2rho_MIN->set_option(reggeon_exchange::k2020_Minimal);
    a2rho_MIN->set_parameters( {alpha0_V, alphaP_V, 0,
                          abeta_rho[0]*gT_MIN[0], abeta_rho[0]*gT_MIN[1], abeta_rho[0]*gT_MIN[2], abeta_rho[0]*gT_MIN[3], abeta_rho[0]*gT_MIN[4],
                          gB_rho[0], gB_rho[1]} );

    amplitude a2omega_MIN = new_amplitude<reggeon_exchange>(kA2, +1, -1, "#omega");
    a2omega_MIN->set_option(reggeon_exchange::k2020_Minimal);
    a2omega_MIN->set_parameters( {alpha0_V, alphaP_V, 0,
                            abeta_omega[0]*gT_MIN[0], abeta_omega[0]*gT_MIN[1], abeta_omega[0]*gT_MIN[2], abeta_omega[0]*gT_MIN[3], abeta_omega[0]*gT_MIN[4], 
                            gB_omega[0], gB_omega[1]} );

    amplitude a2h1_MIN = new_amplitude<reggeon_exchange>(kA2, -1, -1, "#it{h}_{1}");
    a2h1_MIN->set_option(reggeon_exchange::k2020_Minimal);
    a2h1_MIN->set_parameters( {alpha0_A, alphaP_A, 0,
                            abeta_h1[0]*gT_M1[0], abeta_h1[0]*gT_M1[1], abeta_h1[0]*gT_M1[2], abeta_h1[0]*gT_M1[3], abeta_h1[0]*gT_M1[4], 
                            gB_h1[0], gB_h1[1]} );

    amplitude a2rho_TMD = new_amplitude<reggeon_exchange>(kA2, +1, -1, "#rho");
    a2rho_TMD->set_option(reggeon_exchange::kRemoveZero);
    a2rho_TMD->set_option(reggeon_exchange::k2020_TMD);
    a2rho_TMD->set_parameters( {alpha0_V, alphaP_V, 0,
                          abeta_rho[1]*gT_TMD[0], abeta_rho[1]*gT_TMD[1], abeta_rho[1]*gT_TMD[2], abeta_rho[1]*gT_TMD[3], abeta_rho[1]*gT_TMD[4],
                          gB_rho[0], gB_rho[1]} );

    amplitude a2omega_TMD = new_amplitude<reggeon_exchange>(kA2, +1, -1, "#omega");
    a2omega_TMD->set_option(reggeon_exchange::k2020_TMD);
    a2omega_TMD->set_parameters( {alpha0_V, alphaP_V, 0,
                            abeta_omega[1]*gT_TMD[0], abeta_omega[1]*gT_TMD[1], abeta_omega[1]*gT_TMD[2], abeta_omega[1]*gT_TMD[3], abeta_omega[1]*gT_TMD[4], 
                            gB_omega[0], gB_omega[1]} );

    amplitude a2h1_TMD = new_amplitude<reggeon_exchange>(kA2, -1, -1, "#it{h}_{1}");
    a2h1_TMD->set_option(reggeon_exchange::k2020_Minimal);
    a2h1_TMD->set_parameters( {alpha0_A, alphaP_A, 0,
                            abeta_h1[1]*gT_M1[0], abeta_h1[1]*gT_M1[1], abeta_h1[1]*gT_M1[2], abeta_h1[1]*gT_M1[3], abeta_h1[1]*gT_M1[4], 
                            gB_h1[0], gB_h1[1]} );

    amplitude a2_MIN = a2rho_MIN + a2omega_MIN + a2h1_MIN;
    amplitude a2_TMD = a2rho_TMD + a2omega_TMD + a2h1_TMD;

    // ------------------------------------------------------------------------------
    /// f2 Amplitudes

    kinematics kF2 = new_kinematics(M_F2);
    kF2->set_meson_JP(TENSOR);

    amplitude f2rho_MIN = new_amplitude<reggeon_exchange>(kF2, +1, -1, "#rho");
    f2rho_MIN->set_option(reggeon_exchange::kRemoveZero);
    f2rho_MIN->set_option(reggeon_exchange::k2020_Minimal);
    f2rho_MIN->set_parameters( {alpha0_V, alphaP_V, 0,
                          fbeta_rho[0]*gT_MIN[0], fbeta_rho[0]*gT_MIN[1], fbeta_rho[0]*gT_MIN[2], fbeta_rho[0]*gT_MIN[3], fbeta_rho[0]*gT_MIN[4],
                          gB_rho[0], gB_rho[1]} );

    amplitude f2omega_MIN = new_amplitude<reggeon_exchange>(kF2, +1, -1, "#omega");
    f2omega_MIN->set_option(reggeon_exchange::k2020_Minimal);
    f2omega_MIN->set_parameters( {alpha0_V, alphaP_V, 0,
                            fbeta_omega[0]*gT_MIN[0], fbeta_omega[0]*gT_MIN[1], fbeta_omega[0]*gT_MIN[2], fbeta_omega[0]*gT_MIN[3], fbeta_omega[0]*gT_MIN[4], 
                            gB_omega[0], gB_omega[1]} );

    amplitude f2h1_MIN = new_amplitude<reggeon_exchange>(kF2, -1, -1, "#it{h}_{1}");
    f2h1_MIN->set_option(reggeon_exchange::k2020_Minimal);
    f2h1_MIN->set_parameters( {alpha0_A, alphaP_A, 0,
                            fbeta_h1[0]*gT_M1[0], fbeta_h1[0]*gT_M1[1], fbeta_h1[0]*gT_M1[2], fbeta_h1[0]*gT_M1[3], fbeta_h1[0]*gT_M1[4], 
                            gB_h1[0], gB_h1[1]} );

    amplitude f2rho_TMD = new_amplitude<reggeon_exchange>(kF2, +1, -1, "#rho");
    f2rho_TMD->set_option(reggeon_exchange::kRemoveZero);
    f2rho_TMD->set_option(reggeon_exchange::k2020_TMD);
    f2rho_TMD->set_parameters( {alpha0_V, alphaP_V, 0,
                          fbeta_rho[1]*gT_TMD[0], fbeta_rho[1]*gT_TMD[1], fbeta_rho[1]*gT_TMD[2], fbeta_rho[1]*gT_TMD[3], fbeta_rho[1]*gT_TMD[4],
                          gB_rho[0], gB_rho[1]} );

    amplitude f2omega_TMD = new_amplitude<reggeon_exchange>(kF2, +1, -1, "#omega");
    f2omega_TMD->set_option(reggeon_exchange::k2020_TMD);
    f2omega_TMD->set_parameters( {alpha0_V, alphaP_V, 0,
                            fbeta_omega[1]*gT_TMD[0], fbeta_omega[1]*gT_TMD[1], fbeta_omega[1]*gT_TMD[2], fbeta_omega[1]*gT_TMD[3], fbeta_omega[1]*gT_TMD[4], 
                            gB_omega[0], gB_omega[1]} );

    amplitude f2h1_TMD = new_amplitude<reggeon_exchange>(kF2, -1, -1, "#it{h}_{1}");
    f2h1_TMD->set_option(reggeon_exchange::k2020_Minimal);
    f2h1_TMD->set_parameters( {alpha0_A, alphaP_A, 0,
                            fbeta_h1[1]*gT_M1[0], fbeta_h1[1]*gT_M1[1], fbeta_h1[1]*gT_M1[2], fbeta_h1[1]*gT_M1[3], fbeta_h1[1]*gT_M1[4], 
                            gB_h1[0], gB_h1[1]} );

    amplitude f2_MIN = f2rho_MIN + f2omega_MIN + f2h1_MIN;
    amplitude f2_TMD = f2rho_TMD + f2omega_TMD + f2h1_TMD;

    // ------------------------------------------------------------------------------
    /// Plot results

    plotter plotter;

    auto print_dxs = [](plot& p, double Egam, amplitude amp, bool legend = false)
    {
        double s = s_cm(Egam);
        auto dsig = [amp,s](double mt)
        {
            return amp->differential_xsection(s, -mt) * 1E-3; // in mub
        };
        if (!legend) p.add_curve({-amp->get_kinematics()->t_min(s), 2}, dsig);
        else         p.add_curve({-amp->get_kinematics()->t_min(s), 2}, dsig, var_def("#it{E}_{#gamma}", Egam, "GeV"));
    };

    double s;

    plot p1 = plotter.new_plot();
    p1.set_labels("#minus#it{t}  [GeV^{2}]", "d#sigma/d#it{t}  [#mub/GeV^{2}]");
    p1.set_ranges({0, 2.0}, {0.0, 0.6});
    p1.add_data({mt4, sig4}, {mte4, err4}, jpacColor::Blue);
    p1.add_data({mt5, sig5}, {mte5, err5}, jpacColor::Red);
    p1.set_legend_spacing(0.025);
    p1.set_legend(0.7, 0.55);
    print_dxs(p1, 4,   a2_MIN);
    print_dxs(p1, 5,   a2_MIN);
    print_dxs(p1, 8.5, a2_MIN);


    plot p2 = plotter.new_plot();
    p2.set_labels("#minus#it{t}  [GeV^{2}]", "d#sigma/d#it{t}  [#mub/GeV^{2}]");
    p2.set_ranges({0, 2.0}, {0, 0.6});
    p2.add_data({mt4, sig4}, {mte4, err4}, jpacColor::Blue);
    p2.add_data({mt5, sig5}, {mte5, err5}, jpacColor::Red);
    print_dxs(p2, 4,   a2_TMD);
    print_dxs(p2, 5,   a2_TMD);
    print_dxs(p2, 8.5, a2_TMD);


    plot p3 = plotter.new_plot();
    p3.set_labels("#minus#it{t}  [GeV^{2}]", "d#sigma/d#it{t}  [#mub/GeV^{2}]");
    p3.set_ranges({0, 2.0}, {0.0, 1.5});
    p3.add_data({mt3, sig3}, {mte3, err3}, jpacColor::Purple);
    p3.set_legend_spacing(0.025);
    p3.set_legend(0.7, 0.55);
    print_dxs(p3, 4,   f2_MIN);
    print_dxs(p3, 5,   f2_MIN);
    print_dxs(p3, 8.5, f2_MIN);
    p3.color_offset(1);
    print_dxs(p3, 3.4, f2_MIN);

    plot p4 = plotter.new_plot();
    p4.set_labels("#minus#it{t}  [GeV^{2}]", "d#sigma/d#it{t}  [#mub/GeV^{2}]");
    p4.set_ranges({0, 2.0}, {0, 1.5});
    p4.add_data({mt3, sig3}, {mte3, err3}, jpacColor::Purple);
    print_dxs(p4, 4,   f2_TMD);
    print_dxs(p4, 5,   f2_TMD);
    print_dxs(p4, 8.5, f2_TMD);
    p4.color_offset(1);
    print_dxs(p4, 3.4, f2_TMD);

    plotter.combine({2,2}, {p1,p2,p3,p4}, "fig4.pdf");
};  
