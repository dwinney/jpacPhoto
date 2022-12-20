// Methods for quickly parsing J/psi-007 2022 data files into the data_set format
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef JPSI007_HPP
#define JPSI007_HPP

#include "data_set.hpp"

namespace jpacPhoto
{
    inline data_set jpsi007_all()
    {
        auto raw     = import_data<15>("data/jpsi007/jpsi007_all.tsv");
        auto reduced = reshape_data<15,4>(raw, {8, 10, 11, 14});

        data_set wrapped(reduced, "J/psi-007");

        wrapped._lab  = true;
        wrapped._negt = true;

        return wrapped;
    };

    inline data_set jpsi007_setting1()
    {
        auto raw     = import_data<15>("data/jpsi007/jpsi007_setting1.tsv");
        auto reduced = reshape_data<15,4>(raw, {8, 10, 11, 14});

        data_set wrapped(reduced, "J/psi-007 (Setting 1)");

        wrapped._lab  = true;
        wrapped._negt = true;
        
        return wrapped;
    };

    inline data_set jpsi007_setting2()
    {
        auto raw     = import_data<15>("data/jpsi007/jpsi007_setting2.tsv");
        auto reduced = reshape_data<15,4>(raw, {8, 10, 11, 14});

        data_set wrapped(reduced, "J/psi-007 (Setting 2)");

        wrapped._lab  = true;
        wrapped._negt = true;
        
        return wrapped;
    };

    inline data_set jpsi007_setting3()
    {
        auto raw     = import_data<15>("data/jpsi007/jpsi007_setting3.tsv");
        auto reduced = reshape_data<15,4>(raw, {8, 10, 11, 14});

        data_set wrapped(reduced, "J/psi-007 (Setting 3)");

        wrapped._lab  = true;
        wrapped._negt = true;
        
        return wrapped;
    };

    inline data_set jpsi007_setting4()
    {
        auto raw     = import_data<15>("data/jpsi007/jpsi007_setting4.tsv");
        auto reduced = reshape_data<15,4>(raw, {8, 10, 11, 14});

        data_set wrapped(reduced, "J/psi-007 (Setting 4)");

        wrapped._lab  = true;
        wrapped._negt = true;
        
        return wrapped;
    };
};

#endif