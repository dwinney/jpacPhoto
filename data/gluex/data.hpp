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
        auto raw      = import_data<6>("/data/gluex/gluex2022_int.tsv");
        auto reduced  = reshape_data<6,3>(raw, {1,2,3});

        data_set wrapped(reduced, "GlueX 2022");

        // energy values are lab frame energies 
        wrapped._lab  = true;

        // combine the bin-sizes to calculate the error bars
        wrapped._werr = (raw[5] - raw[4]) / 2;

        return wrapped;
    };

    inline data_set gluex_slice0()
    {
        auto raw     = import_data<7>("/data/gluex/gluex2022_diff_E0893.tsv");
        auto reduced = reshape_data<7,4>(raw, {6,1,2,3});

        data_set wrapped(reduced, "GlueX (E = 8.93 GeV)");

        wrapped._lab   = true;
        wrapped._negt  = true;
        wrapped._avg_s = 8.92877;

        // Bin widths 
        wrapped._terr = (raw[5] - raw[4]) / 2;

        return wrapped;
    };

    inline data_set gluex_slice1()
    {
        auto raw     = import_data<7>("/data/gluex/gluex2022_diff_E0985.tsv");
        auto reduced = reshape_data<7,4>(raw, {6,1,2,3});

        data_set wrapped(reduced, "Gluex (E = 9.85 GeV)");

        wrapped._lab   = true;
        wrapped._negt  = true;
        wrapped._avg_s = 9.8583;

        return wrapped;
    };

    inline data_set gluex_slice2()
    {
        auto raw     = import_data<7>("/data/gluex/gluex2022_diff_E1082.tsv");
        auto reduced = reshape_data<7,4>(raw, {6,1,2,3});

        data_set wrapped(reduced, "GlueX (E = 10.82 GeV)");

        wrapped._lab   = true;
        wrapped._negt  = true;
        wrapped._avg_s = 10.82;

        return wrapped;
    };
};

#endif