// Methods for quickly visualizing J/psi-007 data and plotting theoretical curves
// compared to data points
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef JPSI007_PLOTS_HPP
#define JPSI007_PLOTS_HPP

#include "data_set.hpp"
#include "plotter.hpp"
#include "data.hpp"

#include <algorithm>

namespace jpacPhoto
{
    namespace jpsi007
    {
        inline plot plot_slice(plotter& p, int i)
        {
            data_set slice = jpsi007::slice(i); 
            std::string header = slice._id;
            slice._id = "#it{J}/#psi-007"; 

            // Grab the position of upper edge of last t-bin
            int max_i = std::distance(slice._x.begin(), max_element(slice._x.begin(), slice._x.end()));
            double tmax = slice._x[max_i] + slice._xerr[1][max_i];

            plot pdif = p.new_plot();
            pdif.add_data(slice);
            pdif.set_logscale(false, true);
            pdif.set_legend(0.5, 0.65);
            pdif.set_ranges({0, tmax + 0.2}, {1E-3, 6});
            pdif.add_header(header);
            pdif.set_labels("|#it{t} - #it{t}_{min}|  [GeV^{2}]", "#it{d}#sigma/#it{dt}  [nb / GeV^{2}]");

            return pdif;
        };
    };
};

#endif