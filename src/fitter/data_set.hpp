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

namespace jpacPhoto
{
    class data_set
    {
        public:

        data_set(std::array<std::vector<double>,3> data, std::string id = "data_set")
        : _id(id),
          _w(data[0]), _sigma(data[1]), _error(data[2])
        {
            check<3>(data);
        };

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

        // Whether the values stored in _w correspond to invariant energy W = sqrt(s) (false)
        // or lab frame energy Egamma (true)
        bool _lab = false;

        // Whether the momentum transfer values stored in _t correspond to invariant t (false)
        // or t' = t - t_min (true)
        bool _tprime = false;

        // Whether saves values in _t are positive or negative t
        // i.e. -t (true) vs t (false)
        bool _negt   = false;

        private: 

        // Make sure all the vectors are the correct size
        template<int S>
        inline void check(std::array<std::vector<double>,S> data)
        {
            // Grab the size of the first entry
            int N = data[0].size();

            // And compare to the rest
            for (int i = 1; i < S; i++)
            {
                if (data[i].size() != N)
                warning("data_set", "Input vectors of " + _id + " have mismatching sizes!");
                clear();
                return;
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