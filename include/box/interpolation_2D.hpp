// Utility class to create interpolations of 2D data arrays
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef TWOD_INTERP
#define TWOD_INTERP

#include <functional>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <array>
#include <vector>
#include "constants.hpp"
#include <Math/Interpolator.h>

namespace jpacPhoto
{
    class interpolation_2D
    {
        // ---------------------------------------------------------------------------
        public: 

        // Constructor
        interpolation_2D(bool verbose = false)
        : _verbose(verbose)
        {};

        ~interpolation_2D()
        {
            clear_grid();
        };
        
        // Set limits for each variable 
        void set_limits(std::array<double,2> x, std::array<double,2> y)
        {
            _xmin = x[0]; _xmax = x[1];
            _ymin = y[0]; _ymax = y[1];
            _limits_set = true;
        };

        void set_grid_size(int nx,  int ny)
        {
            _xN = nx; _yN = ny;
        }
        
        // Delete the current grid
        void clear_grid()
        {
            _x_values.clear();
            _y_values.clear();
            _f_values.clear();

            for (int i=0; i < _y_slices.size(); i++) delete _y_slices[i];
            _y_slices.clear(); 
        };

        // Take in a function and generate a grid of points based in set grid parameters
        void generate_grid(std::function<double(double, double)> f, bool skip_interp = false);

        // Output the function at a given value
        double eval(double x, double y);

        // Export grid to file 
        void export_grid(std::string filename);

        // or import a grid from a file
        void import_grid(std::string filename);

        // Access function value saved at a given index location of the grid
        double get_data_point(int i, int j){ return _f_values[i][j]; };

        // Whether to show messages like file locations etc
        void set_verbose(bool x){ _verbose = x; };

        // ---------------------------------------------------------------------------
        private: 
        
        // Whether to show command line messages
        bool _verbose = false;

        // Grid parameters
        bool _limits_set = false;
        double _xmin, _xmax, _ymin, _ymax; // Limits of the two variables
        int _xN = 100, _yN = 100;          // Number of points in each variable

        std::vector<double>              _x_values;
        std::vector<std::vector<double>> _y_values;
        std::vector<std::vector<double>> _f_values;
        
        // 1D interpolations of f as function of y at fixed x
        void set_up_slices();
        bool _slices_made = false;
        std::vector<ROOT::Math::Interpolator*> _y_slices;
    };
};

#endif