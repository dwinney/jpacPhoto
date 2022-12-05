// Useful methods for debugging and error handling
//
// ---------------------------------------------------------------------------
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef DEBUG_HPP
#define DEBUG_HPP

#include <iostream>
#include <iomanip>

namespace jpacPhoto
{
    // -----------------------------------------------------------------------

    // Default spacing value
    const int DEBUG_SPACING = 15;

    // Functions for printing to screen instead of having to copy this line all the time
    template<typename T>
    inline void debug(T x)
    {
        std::cout << std::boolalpha; 
        std::cout << x << std::endl;
    };

    template<typename T, typename F>
    inline void debug(T x, F y)
    {
        std::cout << std::boolalpha; 
        std::cout << std::left << std::setw(DEBUG_SPACING) << x;
        std::cout << std::left << std::setw(DEBUG_SPACING) << y << std::endl;
    };

    template<typename T, typename F, typename G>
    inline void debug(T x, F y, G z)
    {
        std::cout << std::boolalpha; 
        std::cout << std::left << std::setw(DEBUG_SPACING) << x;
        std::cout << std::left << std::setw(DEBUG_SPACING) << y;
        std::cout << std::left << std::setw(DEBUG_SPACING) << z << std::endl;
    };

    template<typename T, typename F, typename G, typename H>
    inline void debug(T x, F y, G z, H a)
    {
        std::cout << std::boolalpha; 
        std::cout << std::left << std::setw(DEBUG_SPACING) << x;
        std::cout << std::left << std::setw(DEBUG_SPACING) << y;
        std::cout << std::left << std::setw(DEBUG_SPACING) << z;
        std::cout << std::left << std::setw(DEBUG_SPACING) << a << std::endl;
    };

    template<typename T, typename F, typename G, typename H, typename I>
    inline void debug(T x, F y, G z, H a, I b)
    {
        std::cout << std::boolalpha; 
        std::cout << std::left << std::setw(DEBUG_SPACING) << x;
        std::cout << std::left << std::setw(DEBUG_SPACING) << y;
        std::cout << std::left << std::setw(DEBUG_SPACING) << z;
        std::cout << std::left << std::setw(DEBUG_SPACING) << a;
        std::cout << std::left << std::setw(DEBUG_SPACING) << b << std::endl;
    };

    // ---------------------------------------------------------------------------
    // ERROR Messages

    // Throw an error message then quits code 
    inline void fatal()
    {
        std::cout << std::left << "FATAL ERROR! Quiting..." << std::endl;
        exit( EXIT_FAILURE );
    };

    // Error message with location and reason messages too
    inline void fatal(std::string location, std::string reason = "")
    {
        std::cout << std::left << "FATAL ERROR! " + location + ": " + reason << std::endl;
        std::cout << std::left << "Quiting..." << std::endl;

        exit( EXIT_FAILURE );
    };

    // Warning message does not exit code or returns simply throws a message up
    inline void warning(std::string message)
    {
        std::cout << std::left << "WARNING! " + message << std::endl;
    };

    // Warning message with additional location
    inline void warning(std::string location, std::string message = "")
    {
        std::cout << std::left << "WARNING! " + location + ": " + message << std::endl;
    };

    // Throw an error message and return a value
    template<typename T> 
    inline T error(std::string location, std::string message, T return_value )
    {
        warning(location, message);
        return return_value;
    };

    // Alternatively without a return value, this simply returns void type
    inline void error(std::string location, std::string message)
    {
        warning(location, message);
        return;
    };
};

#endif
