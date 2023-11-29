// Abstract class to represent a function of two variables which generalizes 
// a bottom vertex to inclusive process.
// This is usually a total hadronic cross section or structure function
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#ifndef INCLUSIVE_FUNCTION_HPP
#define INCLUSIVE_FUNCTION_HPP

#include <memory>
#include "constants.hpp"

namespace jpacPhoto
{
    class raw_inclusive_function;
    using inclusive_function = std::shared_ptr<raw_inclusive_function>;

    template<class T>
    inline inclusive_function new_inclusive_function()
    {
        inclusive_function ptr = std::make_shared<T>();
        return ptr;
    };

    template<class T, class A>
    inline inclusive_function new_inclusive_function(A par)
    {
        inclusive_function ptr = std::make_shared<T>(par);
        return ptr;
    };

    template<class T, class A, class B>
    inline inclusive_function new_inclusive_function(A par1, B par2)
    {
        inclusive_function ptr = std::make_shared<T>(par1, par2);
        return ptr;
    };

    template<class T, class A, class B, class C>
    inline inclusive_function new_inclusive_function(A par1, B par2, C par3)
    {
        inclusive_function ptr = std::make_shared<T>(par1, par2, par3);
        return ptr;
    };

    class raw_inclusive_function
    {
        public: 
        
        raw_inclusive_function(){};

        raw_inclusive_function(std::array<double, 2> m)
        : _mB(m[0]), _mT(m[1])
        {};

        // Just need to define how to evaluate the cross-section
        // We may allow an additional variable q2 to consider virtuality effects 
        virtual double evaluate(double s, double q2) = 0;

        // With no specified q2 we default to the on-shell mass of the beam
        inline  double evaluate(double s){ return evaluate(s, _mB*_mB); };

        // Mass of beam and target
        double _mB = 0, _mT = 0;
    };
};

#endif