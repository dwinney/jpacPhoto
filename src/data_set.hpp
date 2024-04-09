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

#include "constants.hpp"
#include "elementwise.hpp"

#include <fstream>
#include <sstream>

namespace jpacPhoto
{

    // Importing data sets we'll need to be able to find the /data/ directory from the 
    // top level one. Thus we need to be able to access the JPACPHOTO environment variable
    inline std::string jpacPhoto_dir()
    {
       // Find the correct data file using the top level repo directory
        std::string top_dir;
        char const * env = std::getenv("JPACPHOTO");
        if ( env == NULL || std::string(env) == "" )
        {
            return error("import_data", "Cannot find environment variable JPACPHOTO!", "");
        }
        return std::string(env);  
    };

    // Import a set of data with N columns with relative path 
    // and full path jpacPhoto_dir/ + rel_path
    template<int N> 
    inline std::array<std::vector<double>,N> import_data(std::string rel_path)
    {
        // Check if rel_path starts with a / or not
        // if not we add one
        if (rel_path.front() != '/') rel_path = "/" + rel_path;

        // Add the top level dir path to get full file path
        std::array<std::vector<double>, N> result;
        std::string file_path = jpacPhoto_dir() + rel_path;
        std::ifstream infile(file_path);

        if (!infile.is_open())
        {
            return error("import_data", "Cannot open file " + file_path + "!", result);
        };

        // Import data!
        std::string line;
        while (std::getline(infile, line))
        {   
            if (line.empty()) continue;        // skips empty lines
            if (line.front() == '#') continue; // Skip comment lines 
            std::istringstream is(line);   

            for (int i = 0; i < N; i++)
            {
                double x;
                is >> x;
                result[i].push_back(x);
            };
        };
            
        return result;
    };

    // Similar to above except that the data is transposed, i.e. rows are the "categories"
    // and the columns are data points. We specify the number of rows in this case
    template<int N>
    inline std::array<std::vector<double>, N> import_transposed(std::string rel_path)
    {
        // Check if rel_path starts with a / or not
        // if not we add one
        if (rel_path.front() != '/') rel_path = "/" + rel_path;

        // Add the top level dir path to get full file path
        std::array<std::vector<double>, N> result;
        std::string file_path = jpacPhoto_dir() + rel_path;
        std::ifstream infile(file_path);

        if (!infile.is_open())
        {
            return error("import_data", "Cannot open file " + file_path + "!", result);
        };

        // Import data!
        for (int i = 0; i < N; i++)
        {   
            std::string line;
            std::getline(infile, line);
            if (line.empty()) continue;        // skips empty lines
            if (line.front() == '#') continue; // Skip comment lines 
            std::istringstream is(line);   

            double x;
            while(is >> x)
            {
                result[i].push_back(x);
            };
        };
            
        return result;
    };

    // If data file has more columns than are actually needed,
    // import with import_data and use this to throw out all but the desired columns
    template<int Nin, int Nout> 
    inline std::array<std::vector<double>,Nout> reshape_data(std::array<std::vector<double>,Nin> data, std::array<int,Nout> to_keep)
    {
        std::array<std::vector<double>, Nout> result;

        for (int i = 0; i < Nout; i++)
        {
            result[i] = data[ to_keep[i] ];
        };
        return result;
    };

    // Make sure all the vectors are the correct size
    template<int S>
    inline int check(std::array<std::vector<double>,S> data, std::string id)
    {
        // Grab the size of the first entry
        int N = data[0].size();
        
        // And compare to the rest
        for (auto column : data)
        {
            if (column.size() != N)
            {
                warning("data_set", "Input vectors of " + id + " have mismatching sizes!");
                return 0;
            };
        };

        return N;
    };

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