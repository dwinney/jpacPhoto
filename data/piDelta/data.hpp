// Interface functions for the Boyarski and GlueX data of pi delta photoproduction
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Universitat Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------

#ifndef PIDELTA_DATA_HPP
#define PIDELTA_DATA_HPP

#include "data_set.hpp"
#include "constants.hpp"

namespace jpacPhoto 
{ 
    namespace piDelta 
    {
        // Old differential data
        inline data_set differential()
        {
            std::string id = "SLAC (1968)";

            auto raw = import_data<3>("/data/piDelta/dcs_pimDelta_Boyarski.txt");
            int  N   = check<3>(raw, id);

            data_set wrapped;
            wrapped._id   = id;
            wrapped._type = 0;
            wrapped._N    = N;
            wrapped._add_to_legend = true;

            wrapped._x       = raw[0]; // -t [GeV2]
            wrapped._z       = raw[1]; // dsig/dt [mub]
            wrapped._zerr[0] = raw[2];
            wrapped._zerr[1] = raw[2];

            // save Egam in extras
            wrapped._extras.push_back(8.0);

            return wrapped;
        };

        // translate the indices of SDME to the index in the data file
        inline int SDME_index(int a, int m, int mp)
        {
            int x = a*100 + m*10 + mp;  
            switch (x)
            {
                case  11: return 1;
                case  31: return 2;
                case  29: return 3;
                case 111: return 4;
                case 133: return 5;
                case 131: return 6;
                case 129: return 7;
                case 231: return 8;
                case 229: return 9;
                default: return -1;
            };
        };

        // translate back
        inline std::array<int,3> SDME_indices(int i)
        {
            switch (i)
            {
                case  1: return {0, 1,  1};
                case  2: return {0, 3,  1};
                case  3: return {0, 3, -1};
                case  4: return {1, 1,  1};
                case  5: return {1, 3,  3};
                case  6: return {1, 3,  1};
                case  7: return {1, 3, -1};
                case  8: return {2, 3,  1};
                case  9: return {2, 3, -1};
                default: return {-1, -1, -1};
            };
        };

        // New SDME data
        inline data_set SDME(int a, int m, int mp)
        {
            std::string id = "GlueX (2024)";

            auto raw = import_data<1+2*9>("/data/piDelta/data_sdmes_helicityframe.txt");
            int    N = check<1+2*9>(raw, id);

            int index = SDME_index(a, m, mp);
            if (index < 0) 
            {
                return error("piDelta::SDME", "Invalid SDME indices passed!", data_set());
            };

            data_set wrapped;
            wrapped._id   = id;
            wrapped._type = index;
            wrapped._N    = N;
            wrapped._add_to_legend = true;

            wrapped._x       = raw[0];             // -t [GeV2]
            wrapped._z       = raw[1+2*(index-1)]; // dsig/dt [mub]
            wrapped._zerr[0] = raw[2+2*(index-1)];
            wrapped._zerr[1] = raw[2+2*(index-1)];

            // save Egam in extras
            wrapped._extras.push_back(8.5);

            return wrapped;
        };

        inline std::vector<data_set> SDMEs()
        {
            std::vector<data_set> ds;
            for (int i = 1; i <= 9; i++) 
            {
                std::array<int,3> ids = SDME_indices(i);
                ds.push_back(SDME(ids[0], ids[1], ids[2]));
            };
            return ds;
        };

        inline data_set beam_asymmetry()
        {
            std::string id = "GlueX (2024)";
            auto raw = import_data<4>("/data/piDelta/data_BSA.txt");
            int    N = check<4>(raw, id);

            data_set wrapped;
            wrapped._id = id;
            wrapped._type = 10;
            wrapped._N = N;
            wrapped._add_to_legend = true;

            wrapped._x = raw[0];       // -t [GeV2]
            wrapped._xerr[0] = raw[1];
            wrapped._xerr[1] = raw[1];
            wrapped._z = raw[2];       // Sigma_4pi
            wrapped._zerr[0] = raw[3];
            wrapped._zerr[1] = raw[3];

            // save Egam in extras
            wrapped._extras.push_back(8.5);
            return wrapped;
        };
    }; 
};

#endif