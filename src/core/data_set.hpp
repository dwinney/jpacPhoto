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

    // ---------------------------------------------------------------------------
    // Element-wise operations on data vectors

    inline std::vector<double>operator+(std::vector<double> lhs, std::vector<double> rhs)
    {
        if (lhs.size() != rhs.size()) return error("Attempted to add two vectors of different sizes!", std::vector<double>());

        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i] + rhs[i] );
        };
        return result;
    };

    inline std::vector<double>operator-(std::vector<double> lhs, std::vector<double> rhs)
    {
        if (lhs.size() != rhs.size()) return error("Attempted to add two vectors of different sizes!", std::vector<double>());

        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i] - rhs[i] );
        };
        return result;
    };

    // Given two vector<double>s of the same size, calculate the average element wise
    inline std::vector<double>operator*( std::vector<double> lhs, double c)
    {
        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i]*c );
        };
        return result;
    };

    inline std::vector<double>operator*(double c, std::vector<double> rhs)
    {
        std::vector<double> result;
        for (int i = 0; i < rhs.size(); i++)
        {
            result.push_back( c*rhs[i] );
        };
        return result;
    };

    inline std::vector<double>operator/( std::vector<double> lhs, double c)
    {
        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i]/c );
        };
        return result;
    };

    class data_set
    {
        public:

        // Empty constructor
        data_set()
        {};

        // Constructor with two columns for integrated data without error bars
        data_set(std::array<std::vector<double>,2> data, std::string id = "data_set")
        : _id(id),
          _w(data[0]), _sigma(data[1])
        {
            check<2>(data);
        };

        // Constructor with three columns for integrated data with error bars
        data_set(std::array<std::vector<double>,3> data, std::string id = "data_set")
        : _id(id),
          _w(data[0]), _sigma(data[1]), _error(data[2])
        {
            check<3>(data);
        };

        // Constructor with four columns for differential data
        data_set(std::array<std::vector<double>,4> data, std::string id = "data_set")
        : _id(id),
          _w(data[0]), _t(data[1]), _sigma(data[2]), _error(data[3])
        { 
            check<4>(data);
        };

        int _N = 0;
        std::string _id = "data_set";

        // Vectors to store energy and momentum transfer variables.
        // Observable and its error
        std::vector<double> _w, _t, _sigma, _error;

        // Other possible vectors to store things like bin sizes, etc
        std::vector<double> _werr, _terr;

        // Whether the values stored in _w correspond to invariant energy W = sqrt(s) (false)
        // or lab frame energy Egamma (true)
        bool _lab = false;

        // Whether the momentum transfer values stored in _t correspond to invariant t (false)
        // or t' = t - t_min (true)
        bool _tprime = false;

        // Whether saves values in _t are positive or negative t
        // i.e. -t (true) vs t (false)
        bool _negt   = false;

        // For a differential set it may be useful to have an average s 
        double _avg_s = 0;

        private: 

        // Make sure all the vectors are the correct size
        template<int S>
        inline void check(std::array<std::vector<double>,S> data)
        {
            // Grab the size of the first entry
            int N = data[0].size();
            
            // And compare to the rest
            for (auto column : data)
            {
                if (column.size() != N)
                {
                    warning("data_set", "Input vectors of " + _id + " have mismatching sizes!");
                    clear();
                    return;
                };
            };

            // If they all pass, save everything
            _N = N;
        };

        // Delete all data stored
        inline void clear()
        {
            _N = 0; 
            _w.clear(); _t.clear(), _sigma.clear(), _error.clear();
        };
    };
};

#endif