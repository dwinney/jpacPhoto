// Methods for quickly parsing Gluex 2022 data files into the data_set format
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef GLUEX_DATA_HPP
#define GLUEX_DATA_HPP

#include "data_set.hpp"
#include "constants.hpp"
#include "kinematics.hpp"

namespace jpacPhoto
{
    namespace gluex
    {
        // Integrated cross-section data as a function of Egamma
        inline data_set integrated()
        {
            std::string id = "GlueX (2023)";

            auto raw  = import_data<6>("/data/jpsip/gluex/gluex2022_int.tsv");
            int  N    = check<6>(raw, id);
        
            data_set wrapped;
            wrapped._id   = id;
            wrapped._type = 0;
            wrapped._N    = N;
            wrapped._add_to_legend = true;

            // energy values are lab frame energies 
            wrapped._x       = raw[1];
            wrapped._xerr[0] = raw[1] - raw[4];
            wrapped._xerr[1] = raw[5] - raw[1];

            // cross section is column 2
            // column 3 is the error
            // compute symmetric upper and lower error bars
            wrapped._z       = raw[2];
            wrapped._zerr[0] = raw[3]; 
            wrapped._zerr[1] = raw[3]; 

            return wrapped;
        };

        // Parses the three slices available for differential cross-section
        inline data_set slice(int sliceid)
        {
            
            std::string id;
            std::string filename;
            double Eavg;

            switch (sliceid)
            {
                case 0: 
                {
                    id = "GlueX (E = 8.93 GeV)";
                    filename = "/data/jpsip/gluex/gluex2022_diff_E0893.tsv";
                    Eavg = 8.92877;
                    break;
                }
                case 1:
                {
                    id = "GlueX (E = 9.85 GeV)";
                    filename = "/data/jpsip/gluex/gluex2022_diff_E0985.tsv";
                    Eavg = 9.8583;
                    break;
                }
                case 2:
                {
                    id = "GlueX (E = 10.82 GeV)";
                    filename = "/data/jpsip/gluex/gluex2022_diff_E1082.tsv";
                    Eavg = 10.8205;
                    break;
                }
                default: return data_set();
            };
            

            auto raw  = import_data<7>(filename);
            int  N    = check<7>(raw, id);

            data_set wrapped;
            wrapped._id    = id;
            wrapped._N     = N;
            wrapped._type = 1;
            wrapped._add_to_legend = true;

            wrapped._x         = raw[1]; // t
            wrapped._xerr[0]   = raw[1] - raw[4];
            wrapped._xerr[1]   = raw[5] - raw[1];
            wrapped._y         = raw[6]; // avg E values
            wrapped._z         = raw[2]; // dsig/dt
            wrapped._zerr[0]   = raw[3];
            wrapped._zerr[1]   = raw[3];
            wrapped._extras.push_back(Eavg);

            return wrapped;
        };

        // Produce a vector with all the 2023 data from GlueX
        // Indices correspond to:
        // [0] slice-0 differential
        // [1] slice-1 differential
        // [2] slice-2 differential
        inline std::vector<data_set> differential()
        {
            std::vector<data_set> data;
            for (int i = 0; i < 3; i++) data.push_back( gluex::slice(i) );
            return data;
        };

        // Produce a vector with all the 2023 data from GlueX
        // Indices correspond to:
        // [0] slice-0 differential
        // [1] slice-1 differential
        // [2] slice-2 differential
        // [3] integrated xsection
        inline std::vector<data_set> all()
        {
            std::vector<data_set> data = differential();
            data.push_back( gluex::integrated() );
            return data;
        };
    };
};

#endif