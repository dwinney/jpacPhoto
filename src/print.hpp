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
#include <functional>
#include <fstream>

namespace jpacPhoto
{
    // Default values
    const int TEXT_WIDTH       = 62;
    const int PRINT_SPACING    = 18;
    const int PRINT_PRECISION  = 9;    
    const int STRING_PRECISION = 3;
    const int PRINT_POINTS     = 100;
    const std::string UNIT_DIV = std::string(PRINT_SPACING, '-');

    // ---------------------------------------------------------------------------   
    // Output an empty line to the terminal
    inline void line()
    {
        std::cout << std::endl;
    };

    // Print out a horizontal line
    inline void divider()
    {
        std::cout << std::string(TEXT_WIDTH, '-') << std::endl;
    };

    inline void divider(int n)
    {
        std::string div;
        for (int i = 0; i < n; i++)
        {
            div = div + UNIT_DIV;
        }
        std::cout << div << std::endl;
    };
    
    inline void dashed_divider()
    {
        std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
    };

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
            std::cout << std::left << std::setw(PRINT_SPACING) << vi;
        };
        std::cout << std::endl;
    };

    template<typename T>
    inline void print_quit(T x)
    {
        std::cout << std::boolalpha << std::left << std::setprecision(9);  
        std::cout << std::setw(PRINT_SPACING) << x << std::endl;
        exit(1);
    };
    template <typename First, typename... Rest>
    inline void print_quit(First first, Rest... rest)
    {
        std::cout << std::boolalpha << std::left << std::setprecision(9);  
        std::cout << std::setw(PRINT_SPACING) << first;
        print_quit(rest...);
    } 

    // ---------------------------------------------------------------------------
    // String operations

    // Produce a string with the format "name = value units"

    template <typename T>
    inline std::string var_def(std::string name, T value, std::string units = "")
    {
        std::stringstream ss;
        ss << std::setprecision(STRING_PRECISION) << name + " = " << value << " " + units;
        return ss.str();
    };

    // Print a string centered on the terminal 
    inline void centered(std::string words)
    {
        int x = words.length();
        int gap_width = (TEXT_WIDTH - x)/2;
        std::cout << std::left << std::setw(gap_width) << "" << std::setw(x) << words << std::setw(gap_width) << "" << std::endl;
    };

    // ---------------------------------------------------------------------------
    // Print things to file

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

    inline void print_to_file(std::array<double,2> boundsx, std::array<double,2> boundsy, std::vector<std::function<double(double,double)>> Fs, std::string file)
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
                for (auto F : Fs) out << std::setw(PRINT_SPACING) << F(xi, xj);
                out << "\n";
            };
        };
        out.close();
    }; 
};

#endif