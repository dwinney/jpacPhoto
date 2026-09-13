// Class and methods for handling data sets used for fitting
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef DATA_SET_HPP
#define DATA_SET_HPP

#include "utilities.hpp"

#include <fstream>
#include <sstream>

namespace jpacPhoto
{
    struct data_set
    {     
        // Number of data points
        int _N = 0;

        std::string _id = "data_set";

        // Each data set should specify what kind of data it is
        // This should match whatever is expected for a particular fitter
        int _type;
        
        // Save up to three data members for each "point"
        // These can include s, t, dsig/dt for example
        std::vector<double> _x, _y, _z;

        // Other possible vectors to store things like bin sizes, etc
        std::array<std::vector<double>, 2> _xerr, _yerr, _zerr;

        // In additon, save any number of extra parameters that may be needed to 
        // identify the data set
        std::vector<double> _extras; 

        // If we want a data entry in the legend when plotting
        bool _add_to_legend = false;
    };

    // For plotters we'll always plot x as the independent variable and z as the dependent one
    // This method swaps x and y
    inline data_set swap_dependent_variable(data_set d)
    {
        data_set new_d;
        new_d._N      = d._N;
        new_d._id     = d._id;
        new_d._type   = d._type;
        new_d._extras = d._extras;
        new_d._x    = d._y;    new_d._xerr = d._yerr;
        new_d._y    = d._x;    new_d._yerr = d._xerr;
        new_d._z    = d._z;    new_d._zerr = d._zerr;
        new_d._add_to_legend = d._add_to_legend;

        return new_d;
    };
};

#endif