// Methods for quickly visualizing data and plotting theoretical curves
// compared to data points
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Instituto de Ciencias Nucleares (ICN),
//               Universidad Nacional Autonoma de México (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#ifndef JPSI_PLOTS_HPP
#define JPSI_PLOTS_HPP

#include "utilities.hpp"
#include "plotter.hpp"

#include "jpsi/data/clas/data.hpp"
#include "jpsi/data/gluex/data.hpp"
#include "jpsi/data/jpsi007/data.hpp"

namespace jpacPhoto
{
    namespace jpsi
    {
        // ------------------------------------------------------------------------------
        // INTEGRATED PLOT

        // Take in a plotter (just to prove colors have been initialized)
        // and produce the pot for the integrated cross_section
        inline plot plot_integrated(plotter& p, amplitude amp = nullptr)
        {
            data_set integrated_gluex = gluex::integrated(); 
            data_set integrated_clas  = clas::integrated(); 
    
            plot pint = p.new_plot();

            if (amp != nullptr)
            {
                pint.add_curve({8, 12}, [&](double Eg){ return amp->integrated_xsection(s_cm(Eg)); }, solid(jpacColor::DarkGrey, amp->id()));
            };
            pint.add_data( integrated_gluex, dot(jpacColor::Blue, "GlueX (2023)"));
            pint.add_data( integrated_clas, square(jpacColor::Red, "CLAS (2026)"));
            pint.set_logscale(false, true);
            pint.set_legend(0.2, 0.65);
            pint.set_ranges({8, 12}, {1E-2, 10});
            pint.set_labels("#it{E}_{#gamma}  [GeV]", 
                            "#sigma(#gamma #it{p} #rightarrow #it{J}/#psi #it{p})  [nb]");
    
            return pint;
        };
    
        // ------------------------------------------------------------------------------
        // GLUEX DIFFERENTIAL PLOTS

        // Plot the differential slices
        inline plot get_gluex_slice(plotter& p, int i, amplitude amp = nullptr)
        {
            data_set slice = gluex::slice(i); 
    
            double avg_E = slice._extras["Eavg"];
            std::array<double,2> yrange;
            switch (i)
            {
                case 0: yrange = {7E-3,  1}; break;
                case 1: yrange = {2E-3,  2}; break;
                case 2: yrange = {4E-4,  5}; break;
            };
            
            plot pdif = p.new_plot();
            pdif.add_data( slice, dot(jpacColor::Blue, "GlueX (2023)"));
            if (amp != nullptr)
            {
                double s = s_cm(avg_E);
                double tmin = -amp->get_kinematics()->t_min(s);
                double tmax = -amp->get_kinematics()->t_max(s);
                pdif.add_curve({tmin, tmax}, [&](double mt){ return amp->differential_xsection(s, -mt); }, solid(jpacColor::DarkGrey, amp->id()));
            };
            pdif.set_logscale(false, true);
            pdif.set_ranges({0, 10.5}, yrange);
            pdif.set_legend(0.6, 0.6);
            pdif.add_header("#LT#it{E}_{#gamma}#GT =", avg_E, "GeV");
            pdif.set_labels("#minus#it{t}  [GeV^{2}]", "#it{d}#sigma/#it{dt}  [nb / GeV^{2}]");
    
            return pdif;
        };

        inline void plot_gluex_differential(plotter& p, amplitude amp = nullptr, std::string filename = "gluex_differential.pdf")
        {
            std::vector<plot> ps;
            for (int i = 0; i<3; i++)
            {
                ps.push_back(get_gluex_slice(p,i, amp));
                if (i != 0) ps.back().set_legend(false);
            };
            p.stack(ps, filename);
        };

        // ------------------------------------------------------------------------------
        // CLAS DIFFERENTIAL PLOTS

        // Plot a specific slice 
        inline plot get_clas_slice(plotter& p, int i, amplitude amp = nullptr)
        {
            data_set slice = clas::slice(i);
    
            double avg_E = slice._extras["Eavg"];
            std::array<double,2> yrange;
            switch (i)
            {
                case 0: yrange = {7E-3,  1}; break;
                case 1: yrange = {2E-3,  2}; break;
                case 2: yrange = {4E-4,  5}; break;
            };

            plot pdif = p.new_plot();
            pdif.add_data( slice, square(jpacColor::Red, "CLAS (2026)"));
            if (amp != nullptr)
            {
                double s = s_cm(avg_E);
                double tmin = -amp->get_kinematics()->t_min(s);
                double tmax = -amp->get_kinematics()->t_max(s);
                pdif.add_curve({tmin, tmax}, [&](double mt){ return amp->differential_xsection(s, -mt); }, solid(jpacColor::DarkGrey, amp->id()));
            };
            pdif.set_logscale(false, true);
            pdif.set_ranges({0, 6.5}, yrange);
            pdif.set_legend(0.6, 0.6);
            pdif.add_header("#LT#it{E}_{#gamma}#GT =", avg_E, "GeV");
            pdif.set_labels("#minus#it{t}  [GeV^{2}]", "#it{d}#sigma/#it{dt}  [nb / GeV^{2}]");
    
            return pdif;
        };


        inline void plot_clas_differential(plotter& p, amplitude amp = nullptr, std::string filename = "clas_differential.pdf")
        {
            std::vector<plot> ps;
            for (int i = 0; i<3; i++)
            {
                ps.push_back(get_clas_slice(p,i, amp));
                if (i != 0) ps.back().set_legend(false);
            };
            p.stack(ps, filename);
        };

        // ------------------------------------------------------------------------------
        // JPSI007 DIFFERENTIAL PLOTS

        inline plot get_jpsi007_slice(plotter& p, int i, amplitude amp = nullptr)
        {
            data_set slice = jpsi007::slice(i); 
            double avg_E = slice._extras["Eavg"];

            // Grab the position of upper edge of last t-bin
            int max_i = std::distance(slice._x.begin(), max_element(slice._x.begin(), slice._x.end()));
            double tmax = slice._x[max_i] + slice._xerr[1][max_i];

            plot pdif = p.new_plot();
            pdif.add_data(slice, triangle(jpacColor::Green, "#it{J}/#psi-007"));
            if (amp != nullptr)
            {
                double s = s_cm(avg_E);
                double tmin = amp->get_kinematics()->t_min(s);
                double tmax = amp->get_kinematics()->t_max(s);
                pdif.add_curve({0, -tmax+tmin}, [&](double mt){ return amp->differential_xsection(s, -mt+tmin); }, solid(jpacColor::DarkGrey, amp->id()));
            };
            pdif.set_logscale(false, true);
            pdif.set_legend(0.6, 0.8);
            pdif.set_ranges({0, tmax + 0.2}, {1E-4, 6});
            pdif.add_header("#LT#it{E}_{#gamma}#GT =", slice._extras["Eavg"], "GeV");
            pdif.set_labels("|#it{t} - #it{t}_{min}|  [GeV^{2}]", "#it{d}#sigma/#it{dt}  [nb / GeV^{2}]");

            return pdif;
        };

        inline void plot_jpsi007_differential(plotter& p, amplitude amp = nullptr, std::string filename = "jpsi007_differential.pdf")
        {
            std::vector<plot> ps;
            for (int i = 1; i<=12; i++)
            {
                ps.push_back(get_jpsi007_slice(p,i, amp));
                if (i != 0) ps.back().set_legend(false);
            };
            p.combine({3,4}, ps, filename);
        };
    };
};

#endif