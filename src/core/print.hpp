// Functions for strings and printing things to the commandline
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef PRINT_HPP
#define PRINT_HPP

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <complex> 

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------   
    // Output an empty line to the terminal
    inline void line()
    {
        std::cout << std::endl;
    };

    // Print out a horizontal line
    inline void divider()
    {
        std::cout << "--------------------------------------------------------------" << std::endl;
    };

    // Default spacing value
    const int PRINT_SPACING = 15;

    template<typename T>
    inline void print(T x)
    {
        std::cout << std::boolalpha << std::left;  
        std::cout << std::setw(PRINT_SPACING) << x << std::endl;
    };

    template <typename First, typename... Rest>
    inline void print(First first, Rest... rest)
    {
        std::cout << std::boolalpha << std::left; 
        std::cout << std::setw(PRINT_SPACING) << first;
        print(rest...);
    } 

    template<typename T>
    inline void print(std::vector<T> v)
    {
        std::cout << std::boolalpha; 
        for (auto vi : v)
        {
            std::cout << std::left << std::setw(PRINT_SPACING) << vi;
        };
        std::cout << std::endl;
    };

    // ---------------------------------------------------------------------------
    // String operations

    const int STRING_PRECISION = 3;

    // Produce a string with the format "name = value units"

    template <typename T>
    inline std::string var_def(std::string name, T value, std::string units = "")
    {
        std::stringstream ss;
        ss << std::setprecision(STRING_PRECISION) << name + " = " << value << " " + units;
        return ss.str();
    };
};

#endif