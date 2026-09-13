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
    // Useful math functions
    
    inline double degrees(double radians){ return radians*DEG2RAD; };
    inline double radians(double degrees){ return degrees/DEG2RAD; };

    // Gamma function of complex argument 'z'
    // optional parameter OPT = 1 can be used to return the log of Gamma(z)
    std::complex<double> cgamma(std::complex<double> z, int OPT = 0);

    // This makes it so we always default to complex regardless of whether the input is an int or double
    template<typename T>
    complex csqrt(T x){ return sqrt(x * XR); };

    unsigned int factorial(unsigned int n);

    // ---------------------------------------------------------------------------
    // Wigner d-func coefficient of leading power
    double wigner_leading_coeff(int j, int lam1, int lam2);

    // Wigner d-function for half-integer spin
    double wigner_d_half(int j, int lam1, int lam2, double theta);

    // Wigner d-function for integer spin
    double wigner_d_int(int j, int lam1, int lam2, double theta);

    // Wigner d-function for integer spin in terms of the cosine of theta not theta
    complex wigner_d_int_cos(int j, int lam1, int lam2, double cos);

    // Legendre function in terms of cosine theta
    double legendre(int l, double z);
    
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


    bool is_real(complex a, double tol = EPS);
    bool is_imaginary(complex a, double tol = EPS);

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

    // Methods that wrap objects in std::cout's and std::endl's to print em to screen
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
    std::string main_dir();

    // Same as above but looks for DESKTOP
    std::string desktop();

    // ---------------------------------------------------------------------------
    // String operations

    // to_string but for doubles
    std::string to_string(double d, uint precision = 8);

    // Print a string centered on the terminal 
    void centered(std::string words);

    // ---------------------------------------------------------------------------
    // Print functions evaluated on a grid to an ascii file

    // Take in a function and print an ascii file of values on a grid
    void print_to_file(std::array<double, 2> bounds, std::function<double(double)> F, std::string file);
    void print_to_file(std::array<double, 2> bounds, std::vector<std::function<double(double)>> Fs, std::string file);
    void print_to_file(std::array<double, 2> boundsx, std::array<double,2> boundsy, std::function<double(double,double)> F, std::string file);

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

    std::vector<double> multiply_elementwise(std::vector<double> in1, std::vector<double> in2);
    std::vector<double> square_elementwise(std::vector<double> in);
    std::vector<double> real(std::vector<complex> vx);
    std::vector<double> imag(std::vector<complex> vx);
   
    // ---------------------------------------------------------------------------
    // Read in data into vectors

    // Import a set of data with N columns with relative path 
    // and full path jpacPhoto_dir/ + rel_path
    template<int N> 
    inline std::array<std::vector<double>,N> import_data(std::string rel_path)
    {
        // Check if rel_path starts with a / or not
        // if not we add one
        if (rel_path.front() != '/') rel_path = "/" + rel_path;

        // Add the top level dir path to get full file path
        std::array<std::vector<double>, N> result;
        std::string file_path = main_dir() + rel_path;
        std::ifstream infile(file_path);

        if (!infile.is_open())
        {
            return error("import_data - Cannot open file " + file_path + "!", result);
        };

        // Import data!
        std::string line;
        while (std::getline(infile, line))
        {   
            if (line.empty()) continue;        // skips empty lines
            if (line.front() == '#') continue; // Skip comment lines 
            std::istringstream is(line);   

            for (int i = 0; i < N; i++)
            {
                double x;
                is >> x;
                result[i].push_back(x);
            };
        };
            
        return result;
    };

    // Similar to above except that the data is transposed, i.e. rows are the "categories"
    // and the columns are data points. We specify the number of rows in this case
    template<int N>
    inline std::array<std::vector<double>, N> import_transposed(std::string rel_path)
    {
        // Check if rel_path starts with a / or not
        // if not we add one
        if (rel_path.front() != '/') rel_path = "/" + rel_path;

        // Add the top level dir path to get full file path
        std::array<std::vector<double>, N> result;
        std::string file_path = main_dir() + rel_path;
        std::ifstream infile(file_path);

        if (!infile.is_open())
        {
            return error("import_data - Cannot open file " + file_path + "!", result);
        };

        // Import data!
        for (int i = 0; i < N; i++)
        {   
            std::string line;
            std::getline(infile, line);
            if (line.empty()) continue;        // skips empty lines
            if (line.front() == '#') continue; // Skip comment lines 
            std::istringstream is(line);   

            double x;
            while(is >> x)
            {
                result[i].push_back(x);
            };
        };
            
        return result;
    };

    // If data file has more columns than are actually needed,
    // import with import_data and use this to throw out all but the desired columns
    template<int Nin, int Nout> 
    inline std::array<std::vector<double>,Nout> reshape_data(std::array<std::vector<double>,Nin> data, std::array<int,Nout> to_keep)
    {
        std::array<std::vector<double>, Nout> result;

        for (int i = 0; i < Nout; i++)
        {
            result[i] = data[ to_keep[i] ];
        };
        return result;
    };

    // Make sure all the vectors are the correct size
    template<int S>
    inline int check(std::array<std::vector<double>,S> data, std::string id)
    {
        // Grab the size of the first entry
        int N = data[0].size();
        
        // And compare to the rest
        for (auto column : data)
        {
            if (column.size() != N)
            {
                warning("data_set - Input vectors of " + id + " have mismatching sizes!");
                return 0;
            };
        };

        return N;
    };

    // ---------------------------------------------------------------------------
    // Simple functions to calculate numerical derivatives up to error ~ h^4
    // Coefficients taken from https://en.wikipedia.org/wiki/Finite_difference_coefficient
    
    // Central difference 
    template<typename T> 
    inline T central_difference_derivative(uint n, std::function<T(double)> F, double x, double h = 1E-3)
    {
        T f3m, f2m, fm, f, fp, f2p, f3p;
        f   = F(x);
        fp  = F(x+h),    fm = F(x-h);
        f2p = F(x+2*h), f2m = F(x-2*h);
        if (n >= 3) { f3p = F(x+3*h); f3m = F(x-3*h); };
        
        T num;
        switch (n)
        {
            case 0  : return F(x);
            case 1  : num =          +f2m/12  -2*fm/3          +2*fp/3   -f2p/12;          break;
            case 2  : num =          -f2m/12  +4*fm/3  -5*f/2  +4*fp/3   -f2p/12;          break;
            case 3  : num = +f3m/8   -f2m    +13*fm/8         -13*fp/8   +f2p     -f3p/8;  break;
            case 4  : num = -f3m/6 +2*f2m    -13*fm/2 +28*f/3 -13*fp/2 +2*f2p     -f3p/6;  break;
            default : 
            {
                warning("central_difference_derivative - Order n = "+to_string(n)+" derivatives not implemented!"); 
                return NaN<T>();
            };
        };
        return num / pow(h, n);
    };

    // Forward difference
    template<typename T> 
    inline T forward_difference_derivative(uint n, std::function<T(double)> F, double x, double h = 1E-3)
    {
        if (n == 0) return F(x);

        std::vector<double> c;
        switch (n)
        {
            // O(h^4)
            case 1  : c = {-25./12, 4., -3., 4./3, -1./4}; break;
            case 2  : c = {15./4, -77./6, 107./6, -13., 61./12, -5./6}; break;
            case 3  : c = {-49./8, 29., -461./8, 62., -307./8, 13., -15./8}; break;
            case 4  : c = {28./3, -111./2, 142., -1219/6., 176., -185./2, 82./3 ,-7./2}; break;
            default : 
            {
                warning("forward_difference_derivative - Order n = "+to_string(n)+" derivatives not implemented!"); 
                return NaN<T>();
            };
        };
        T num = 0;
        for (int i = 0; i < c.size(); i++) num += c[i] * F(x+i*h);

        return num / pow(h, n);
    };
    
    // and finally backward finite difference
    template<typename T> 
    inline T backward_difference_derivative(uint n, std::function<T(double)> F, double x, double h = 1E-3)
    {
        if (n == 0) return F(x);

        std::vector<double> c;
        switch (n)
        {
            // These are O(h^4)
            case 1  : c = {-25./12, 4., -3., 4./3, -1./4}; break;
            case 2  : c = {15./4, -77./6, 107./6, -13., 61./12, -5./6}; break;
            case 3  : c = {-49./8, 29., -461./8, 62., -307./8, 13., -15./8}; break;
            case 4  : c = {28./3, -111./2, 142., -1219/6., 176., -185./2, 82./3 ,-7./2}; break;
            default : 
            {
                warning("forward_difference_derivative - Order n = "+to_string(n)+" derivatives not implemented!"); 
                return NaN<T>();
            };
        };
        T num = 0;
        for (int i = 0; i < c.size(); i++) num += c[i] * F(x-i*h);

        return pow(-1, n) * num / pow(h, n);
    };

    // Mixed central derivative of a function of 2 variables d2F(x,y)/dxdy
    template<typename T>
    inline T mixed_central_derivatives(std::function<T(double,double)> F, std::array<double,2> xs, double e)
    {
        double x = xs[0], y = xs[1];
        T f2p = central_difference_derivative<T>(1, [&](double s){ return F(s, y+2*e); },  x, e);
        T fp  = central_difference_derivative<T>(1, [&](double s){ return F(s, y+e);   },  x, e);
        T fm  = central_difference_derivative<T>(1, [&](double s){ return F(s, y-e);   },  x, e);
        T f2m = central_difference_derivative<T>(1, [&](double s){ return F(s, y-2*e); },  x, e);

        return (+f2m/12.-2.*fm/3.+2.*fp/3.-f2p/12.)/e;
    };

    // Mixed forward derivative of a function of 2 variables d2F(x,y)/dxdy
    template<typename T>
    inline T mixed_forward_derivatives(std::function<T(double,double)> F, std::array<double,2> xs, double e)
    {
        double x = xs[0], y = xs[1];
        std::array<double,5> c = {-25./12, 4., -3., 4./3, -1./4};
        T sum = 0.;
        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                sum += c[i]*c[j]*F(x+i*e, y+j*e);
            };
        };
        return sum/e/e;
    };

    template<typename T>
    inline T mixed_backward_derivatives(std::function<T(double,double)> F, std::array<double,2> xs, double e)
    {
        double x = xs[0], y = xs[1];
        std::array<double,5> c = {-25./12, 4., -3., 4./3, -1./4};
        T sum = 0.;
        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                sum += c[i]*c[j]*F(x-i*e, y-j*e);
            };
        };
        return sum/e/e;
    };
};
// ---------------------------------------------------------------------------

#endif
