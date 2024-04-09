// Interface functions for the Boyarski and GlueX data of pi delta photoproduction
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Universitat Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------

#ifndef PIDELTA_PLOTS_HPP
#define PIDELTA_PLOTS_HPP

#include "data_set.hpp"
#include "plotter.hpp"
#include "piDelta/data.hpp"
#include "constants.hpp"

namespace jpacPhoto 
{ 
    namespace piDelta 
    {
        inline plot plot_differential(plotter& plotr)
        {
            plot p = plotr.new_plot();

            p.add_data(piDelta::differential());
            p.set_labels("#minus#it{t} [GeV]", "d#sigma/d#it{t}  [#mub]");
            p.set_logscale(false, true);
            p.add_header("#it{E}_{#gamma} = 8 GeV");
            p.set_legend(0.6, 0.5);
            p.set_ranges({1E-3, 1.08}, {5E-2, 7});
            return p;
        };

        inline plot plot_SDME(plotter& plotr, int a, int m, int mp)
        {
            plot p = plotr.new_plot();
            p.add_data(piDelta::SDME(a, m, mp));
            p.add_header("#it{E}_{#gamma} = 8.5 GeV");
            p.set_legend(0.7, 0.2);
            p.set_ranges({0, 1}, {-0.5, 0.5});

            std::string label = "#rho^{" + std::to_string(a) + "}_{ " + std::to_string(m) + std::to_string(mp) + "}";
            if (a == 2) label = "Im " + label;
            p.set_labels("#minus #it{t}  [GeV^{2}]", label);
            return p;
        };

        inline std::vector<plot> plot_SDMEs(plotter & plotr)
        {
            std::vector<plot> ps;
            for (int i = 1; i <= 9; i++) 
            {
                std::array<int,3> ids = SDME_indices(i);
                ps.push_back(plot_SDME(plotr, ids[0], ids[1], ids[2]));
            };
            return ps;
        };

        inline plot plot_beam_asymmetry(plotter& plotr)
        {
            plot p = plotr.new_plot();
            p.add_data(piDelta::beam_asymmetry());
            p.add_header("#it{E}_{#gamma} = 8.5 GeV");
            p.set_legend(0.7, 0.2);
            p.set_ranges({0, 1.2}, {-1, 1});
            p.set_labels("#minus #it{t}  [GeV^{2}]", 
                            "#Sigma_{4#pi}");
            return p;
        };

    }; 
};

#endif