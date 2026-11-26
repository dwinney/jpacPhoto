// Header file with useful functions
//
// ---------------------------------------------------------------------------
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <ios>
#include "constants.hpp"
#include "jpacUtils/print.hpp"

namespace jpacPhoto
{
    using namespace jpacUtils;
    
    // ---------------------------------------------------------------------------
    // Frame conversion methods (specifically for photoproduction)

    // Photon lab energy
    double E_beam(double W);
    // Center of mass energy given beam energy
    double W_cm(double egam);
    // Center of mass energy given beam energy
    double s_cm(double egam);

    // Kinematic function for 2->2 scattering (see eq. 5.23 in Byckling & Kajantie)
    double G(double x, double y, double z, double u, double v, double w);

    // ---------------------------------------------------------------------------
    // Importing data sets we'll need to be able to find the main directory from the 
    // top level one. Thus we need to be able to access the appropriate environment variable
    std::string main_dir();

    // Same as above but looks for DESKTOP
    std::string desktop();

    // Sometimes we'll want to work on an analysis outside the main JPACPHOTO directory
    // this looks for the variable JPACPHOTO_WORKING_DIR
    std::string working_dir();

    // ---------------------------------------------------------------------------
    // Element-wise operations on data vectors
    
    // Given two vector<double>s of the same size, calculate the average element wise
    template<typename T>
    inline std::vector<T> operator*( std::vector<T> lhs, double c)
    {
        std::vector<T> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i]*c );
        };
        return result;
    };

    template<typename T>
    inline std::vector<T> operator*(double c, std::vector<T> rhs)
    {
        std::vector<T> result;
        for (int i = 0; i < rhs.size(); i++)
        {
            result.push_back( c*rhs[i] );
        };
        return result;
    };

    template<typename T>
    inline std::vector<T> operator/( std::vector<T> lhs, double c)
    {
        std::vector<T> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i]/c );
        };
        return result;
    };

    template<typename T>
    inline std::vector<T> operator-(const std::vector<T> & x)
    {
        return -1. * x;
    };

    template<typename T>
    inline std::vector<T> operator/=(const std::vector<T> & x, double c)
    {
        return x/c;
    };

    template<typename T>
    inline std::vector<T> operator*=(const std::vector<T> & x, double c)
    {
        return x*c;
    };
    
    template<typename T>
    inline std::vector<T> operator+(std::vector<T> lhs, std::vector<T> rhs)
    {
        if (lhs.size() != rhs.size()) return error("Attempted to add two vectors of different sizes!", std::vector<double>());

        std::vector<T> result;
        for (uint i = 0; i < lhs.size(); i++) result.push_back( lhs[i] + rhs[i] );
        return result;
    };

    template<typename T>
    inline std::vector<T> operator-(std::vector<T> lhs, std::vector<T> rhs)
    {
        if (lhs.size() != rhs.size()) return error("Attempted to add two vectors of different sizes!", std::vector<double>());

        std::vector<double> result;
        for (uint i = 0; i < lhs.size(); i++) result.push_back( lhs[i] - rhs[i] );
        return result;
    };

    // Add a constant to all elements of a vector
    template<typename T>
    inline std::vector<T> operator+(double lhs, std::vector<T> rhs)
    {
        std::vector<T> result;
        for (uint i = 0; i < rhs.size(); i++) result.push_back( lhs + rhs[i] );
        return result;
    };   
    
    template<typename T>
    inline std::vector<T> operator-(double lhs, std::vector<T> rhs) 
    {
        return lhs + (-rhs);
    };

    template<typename T>
    inline std::vector<T> operator+(std::vector<T> lhs, double rhs)
    {
        std::vector<T> result;
        for (uint i = 0; i < lhs.size(); i++) result.push_back( rhs + lhs[i] );
        return result;
    };

    template<typename T>
    inline std::vector<T> operator-(std::vector<T> lhs, double rhs) 
    {
        return lhs + (-rhs);
    };

    inline std::vector<double> multiply_elementwise(std::vector<double> in1, std::vector<double> in2)
    {
        if (in1.size() != in2.size()) warning("multiply_elementwise()", "Input vectors not the same size!");

        std::vector<double> out;
        for (uint i = 0; i < in1.size(); i++)
        {
            out.push_back(in1[i]*in2[i]);
        }
        return out;
    };
};
// ---------------------------------------------------------------------------

#endif
