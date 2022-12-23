// Useful methods for debugging and error handling
//
// ---------------------------------------------------------------------------
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef DEBUG_HPP
#define DEBUG_HPP

#include <ios>
#include <iostream>
#include <iomanip>
#include <vector>

namespace jpacPhoto
{
    // -----------------------------------------------------------------------
    
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
    inline void warning(std::string location, std::string message)
    {
        std::cout << std::left << "WARNING! " + location + ": " + message << std::endl;
    };

        // Throw an error message without location and return a value
    template<typename T> 
    inline T error(std::string message, T return_value )
    {
        warning(message);
        return return_value;
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
