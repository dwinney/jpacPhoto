// Header file with useful functions
//
// ---------------------------------------------------------------------------
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cmath>
#include <complex>
#include <limits>
#include <ios>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <vector>
#include <functional>
#include <fstream>

#include "constants.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // the complex type is a liTtle dim in c++ so we need to define int & bool multiplication

    complex operator*(const int& c, const complex& rhs);
    complex operator*(const complex& lhs, const int& c);
    complex operator*(const bool& c, const complex& rhs);
    complex operator*(const complex& lhs, const bool& c);
    complex operator/(const complex&c, const int& i);
    complex operator/(const int& i, const complex&c);
    complex operator+(const complex&c, const int& i);
    complex operator+(const int& i, const complex & c);
    complex operator-(const complex&c, const int& i);
    complex operator-(const int& i, const complex & c);

    // This makes it so we always default to complex regardless of whether the input is an int or double
    template<typename T>
    complex csqrt(T x){ return sqrt(x * XR); };

    inline unsigned int factorial(unsigned int n) 
    {
        if (n == 0)
        return 1;
        return n * factorial(n - 1);
    };

    // ---------------------------------------------------------------------------
    // Frame conversion methods (specifically for photoproduction)

    // Photon lab energy
    double E_beam(double W);
    // Center of mass energy given beam energy
    double W_cm(double egam);
    // Center of mass energy given beam energy
    double s_cm(double egam);

    // ---------------------------------------------------------------------------
    // Kallen Triangle function

    // Only way to get a double or int Kallen is if all inputs are double/int
    template<typename T>
    inline T Kallen(T x, T y, T z)
    {
        return x*x + y*y + z*z - 2. * (x*y + x*z + y*z);
    };

    // If any of them are complex, return complex
    complex Kallen(complex z, double a, double b);
    complex Kallen(double a, complex z, double b);
    complex Kallen(double a, double b, complex z);

    // Kinematic function for 2->2 scattering (see eq. 5.23 in Byckling & Kajantie)
    double G(double x, double y, double z, double u, double v, double w);

    // ---------------------------------------------------------------------------
    // Function for easier comparison of doubles using the EPS value defined above
    // be careful when using this in general purposes since its a fixed-tolerance comparision and not always appropriate

    bool are_equal(double a,  double b,  double tol = EPS);
    bool are_equal(complex a, complex b, double tol = EPS);

    // Aliases for special cases of the above
    bool is_zero(double a, double tol = EPS);
    bool is_zero(complex a, double tol = EPS);

    // ---------------------------------------------------------------------------
    // ERROR Messages
    
    // Error message with location and reason messages too
    void fatal(std::string reason = "");

    // Warning message does not exit code or returns simply throws a message up
    void warning(std::string message);

    // Throw an error message without location and return a value
    template<typename T> 
    inline T error(std::string message, T return_value )
    {
        warning(message);
        return return_value;
    };

    // ---------------------------------------------------------------------------   
    // Methods that make printing messages easier 

    // Output an empty line to the terminal
    void line();

    // Print out a horizontal line
    void divider();
    void divider(int n);
    void dashed_divider();

    template<typename T>
    inline void print(T x)
    {
        std::cout << std::boolalpha << std::left << std::setprecision(9);  
        std::cout << std::setw(PRINT_SPACING) << x << std::endl;
    };

    template <typename First, typename... Rest>
    inline void print(First first, Rest... rest)
    {
        std::cout << std::boolalpha << std::left << std::setprecision(9);  
        std::cout << std::setw(PRINT_SPACING) << first;
        print(rest...);
    } 

    template<typename T>
    inline void print(std::vector<T> v)
    {
        std::cout << std::boolalpha << std::setprecision(9);  
        for (auto vi : v)
        {
            std::cout << std::left << std::setw(PRINT_SPACING) << vi << std::endl;
        };
        std::cout << std::endl;
    };

    // ---------------------------------------------------------------------------
    // Importing data sets we'll need to be able to find the main directory from the 
    // top level one. Thus we need to be able to access the appropriate environment variable
    inline std::string main_dir()
    {
       // Find the correct data file using the top level repo directory
        std::string top_dir;
        char const * env = std::getenv("JPACPHOTO");
        if ( env == NULL || std::string(env) == "" )
        {
            return error("main_dir(): Cannot find environment variable JPACPHOTO!", "");
        }
        return std::string(env);  
    };

    // Same as above but looks for DESKTOP
    inline std::string desktop()
    {
       // Find the correct data file using the top level repo directory
        std::string top_dir;
        char const * env = std::getenv("DESKTOP");
        if ( env == NULL || std::string(env) == "" )
        {
            return error("desktop(): Cannot find environment variable DESKTOP!", "");
        }
        return std::string(env);  
    };

    // ---------------------------------------------------------------------------
    // String operations

    // Print a string centered on the terminal 
    inline void centered(std::string words)
    {
        int x = words.length();
        int gap_width = (TEXT_WIDTH - x)/2;
        std::cout << std::left << std::setw(gap_width) << "" << std::setw(x) << words << std::setw(gap_width) << "" << std::endl;
    };

    // ---------------------------------------------------------------------------
    // Print functions evaluated on a grid to an ascii file

    // Take in a function and print an ascii file of values on a grid
    inline void print_to_file(std::array<double, 2> bounds, std::function<double(double)> F, std::string file)
    {
        std::ofstream out;
        out.open(file);
        
        for (int i = 0; i < PRINT_POINTS; i++)
        {
            double xi = bounds[0] + double(i) * (bounds[1] - bounds[0]) / double(PRINT_POINTS-1);
            out << std::left << std::setw(PRINT_SPACING) << xi << std::setw(PRINT_SPACING) << F(xi) << "\n";
        };
        out.close();
    }; 

    inline void print_to_file(std::array<double, 2> bounds, std::vector<std::function<double(double)>> Fs, std::string file)
    {
        std::ofstream out;
        out.open(file);
        
        for (int i = 0; i < PRINT_POINTS; i++)
        {
            double xi = bounds[0] + double(i) * (bounds[1] - bounds[0]) / double(PRINT_POINTS-1);
            out << std::left << std::setw(PRINT_SPACING) << xi; 
            for (auto F : Fs) out << std::setw(PRINT_SPACING) << F(xi);
            out << "\n";
        };
        out.close();
    }; 

    inline void print_to_file(std::array<double,2> boundsx, std::array<double,2> boundsy, std::function<double(double,double)> F, std::string file)
    {
        std::ofstream out;
        out.open(file);
        
        for (int i = 0; i < PRINT_POINTS; i++)
        {
            double xi = boundsx[0] + double(i) * (boundsx[1] - boundsx[0]) / double(PRINT_POINTS-1);
            for (int j = 0; j < PRINT_POINTS; j++)
            {
                double xj = boundsy[0] + double(j) * (boundsy[1] - boundsy[0]) / double(PRINT_POINTS-1);
                out << std::left;
                out << std::setw(PRINT_SPACING) << xi;
                out << std::setw(PRINT_SPACING) << xj;
                out << std::setw(PRINT_SPACING) << F(xi, xj);
                out << "\n";
            };
        };
        out.close();
    }; 
};
// ---------------------------------------------------------------------------

#endif
