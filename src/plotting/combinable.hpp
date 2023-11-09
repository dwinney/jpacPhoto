// Abstract class to interface with the plotter::combine method
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef COMBINABLE_HPP
#define COMBINABLE_HPP

namespace jpacPhoto
{
    class combinable
    {
        public:

        combinable(){};

        virtual void combine_draw(double scale) = 0;
    };
};

#endif
