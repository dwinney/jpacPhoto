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

namespace jpacPhoto
{
    inline data_set gluex_integrated2022()
    {
        std::string id = "GlueX 2022";

        auto raw  = import_data<6>("/data/gluex/gluex2022_int.tsv");
        int  N    = check<6>(raw, id);
    
        data_set wrapped;
        wrapped._id   = id;
        wrapped._N    = N;
        wrapped._lab  = true;
        wrapped._type = integrated_data;

        // energy values are lab frame energies 
        wrapped._w       = raw[1];
        wrapped._werr[0] = raw[1] - raw[4];
        wrapped._werr[1] = raw[5] - raw[1];

        // cross section is column 2
        // column 3 is the error
        // compute symmetric upper and lower error bars
        wrapped._obs       = raw[2];
        wrapped._obserr[0] = raw[3]/2; 
        wrapped._obserr[1] = raw[3]/2; 

        return wrapped;
    };

    inline data_set gluex_slice0()
    {
        std::string id = "GlueX (E = 8.93 GeV)";

        auto raw  = import_data<7>("/data/gluex/gluex2022_diff_E0893.tsv");
        int  N    = check<7>(raw, id);

        data_set wrapped;
        wrapped._id    = id;
        wrapped._N     = N;
        wrapped._type  = differential_data;
        wrapped._lab   = true;
        wrapped._negt  = true;
        wrapped._avg_w = 8.92877;

        // Bin widths 
        wrapped._w         = raw[6];
        wrapped._t         = raw[1];
        wrapped._terr[0]   = raw[1] - raw[4];
        wrapped._terr[1]   = raw[5] - raw[1];
        wrapped._obs       = raw[2];
        wrapped._obserr[0] = raw[3]/2;
        wrapped._obserr[1] = raw[3]/2;

        return wrapped;
    };

    inline data_set gluex_slice1()
    {
        std::string id = "GlueX (E = 9.85 GeV)";

        auto raw  = import_data<7>("/data/gluex/gluex2022_diff_E0985.tsv");
        int  N    = check<7>(raw, id);

        data_set wrapped;
        wrapped._id    = id;
        wrapped._N     = N;
        wrapped._type  = differential_data;
        wrapped._lab   = true;
        wrapped._negt  = true;
        wrapped._avg_w = 9.8583;

        // Bin widths 
        wrapped._w         = raw[6];
        wrapped._t         = raw[1];
        wrapped._terr[0]   = raw[1] - raw[4];
        wrapped._terr[1]   = raw[5] - raw[1];
        wrapped._obs       = raw[2];
        wrapped._obserr[0] = raw[3]/2;
        wrapped._obserr[1] = raw[3]/2;

        return wrapped;
    };

    inline data_set gluex_slice2()
    {
        std::string id = "GlueX (E = 10.82 GeV)";

        auto raw  = import_data<7>("/data/gluex/gluex2022_diff_E1082.tsv");
        int  N    = check<7>(raw, id);

        data_set wrapped;
        wrapped._id    = id;
        wrapped._N     = N;
        wrapped._type  = differential_data;
        wrapped._lab   = true;
        wrapped._negt  = true;
        wrapped._avg_w = 10.8205;

        // Bin widths 
        wrapped._w         = raw[6];
        wrapped._t         = raw[1];
        wrapped._terr[0]   = raw[1] - raw[4];
        wrapped._terr[1]   = raw[5] - raw[1];
        wrapped._obs       = raw[2];
        wrapped._obserr[0] = raw[3]/2;
        wrapped._obserr[1] = raw[3]/2;

        return wrapped;
    };

};

#endif