// Header file with useful functions
//
// ---------------------------------------------------------------------------
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "utilities.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // the complex type is a liTtle dim in c++ so we need to define int & bool multiplication

    complex operator*(const int& c, const complex& rhs)
    {
        return complex(c*rhs.real(), c*rhs.imag());
    };

    complex operator*(const complex& lhs, const int& c)
    {
        return complex(c*lhs.real(), c*lhs.imag());
    };

    complex operator*(const bool& c, const complex& rhs)
    {
        return (c) ? rhs : 0.;
    };

    complex operator*(const complex& lhs, const bool& c)
    {
        return (c) ? lhs : 0.;
    };

    complex operator/(const complex&c, const int& i)
    {
        return (1./i)*c;
    };

    complex operator/(const int& i, const complex&c)
    {
        return (1./c)*i;
    };

    complex operator+(const complex&c, const int& i)
    {
        return c + XR*i;
    };

    complex operator+(const int& i, const complex & c)
    {
        return XR*i + c;
    };

    complex operator-(const complex&c, const int& i)
    {
        return c - XR*i;
    };

    complex operator-(const int& i, const complex & c)
    {
        return XR*i - c;
    };

    // ---------------------------------------------------------------------------
    // Frame conversion methods (specifically for photoproduction)

    // Photon lab energy
    double E_beam(double W)
    {
        return (W*W / M_PROTON - M_PROTON) / 2.;
    };

    // Center of mass energy given beam energy
    double W_cm(double egam)
    {
        return sqrt(M_PROTON * (2. * egam + M_PROTON));
    };

    // Center of mass energy given beam energy
    double s_cm(double egam)
    {
        return M_PROTON * (2. * egam + M_PROTON);
    };

    // ---------------------------------------------------------------------------
    // Kallen Triangle function

    // If any of them are complex, return complex
    complex Kallen(complex z, double a, double b) { return Kallen<complex>(z, XR*a, XR*b); };
    complex Kallen(double a, complex z, double b) { return Kallen<complex>(XR*a, z, XR*b); };
    complex Kallen(double a, double b, complex z) { return Kallen<complex>(XR*a, XR*b, z); };

    // Kinematic function for 2->2 scattering (see eq. 5.23 in Byckling & Kajantie)
    double G(double x, double y, double z, double u, double v, double w)
    {
        return  x*x*y + x*y*y + z*z*u + z*u*u + v*v*w + v*w*w 
              + x*z*w + x*u*v + y*z*w + y*u*w + y*z*v - y*z*w
              - x*y*(z + u + v + w) - z*u*(x + y + v + w) - v*w*(x + y + z + u);
    };

     // ---------------------------------------------------------------------------
    // Function for easier comparison of doubles using the EPS value defined above
    // be careful when using this in general purposes since its a fixed-tolerance comparision and not always appropriate

    bool are_equal(double a, double b, double tol)
    {
        return ( std::abs(a - b) < tol );
    }

    // Same thing for comparing complex doubles
    bool are_equal(complex a, complex b, double tol)
    {
        return (are_equal(real(a), real(b), tol) && are_equal(imag(a), imag(b), tol));
    };

    // Aliases for special cases of the above
    bool is_zero(double a, double tol)
    {
        return (std::abs(a) < tol);
    };

    bool is_zero(complex a, double tol)
    {
        return (std::abs(a) < tol);
    };

    // ---------------------------------------------------------------------------
    // ERROR Messages
    
    // Error message with location and reason messages too
    void fatal(std::string reason)
    {
        std::cout << std::left << "FATAL ERROR! " + reason << std::endl;
        std::cout << std::left << "Quiting..." << std::endl;
        exit( EXIT_FAILURE );
    };

    // Warning message does not exit code or returns simply throws a message up
    void warning(std::string message)
    {
        std::cout << std::left << "WARNING! " + message << std::endl;
    };

    // ---------------------------------------------------------------------------   
    // Methods that make printing messages easier 

    // Output an empty line to the terminal
    void line()
    {
        std::cout << std::endl;
    };

    // Print out a horizontal line
    void divider()
    {
        std::cout << std::string(TEXT_WIDTH, '-') << std::endl;
    };

    void divider(int n)
    {
        std::string div;
        for (int i = 0; i < n; i++)
        {
            div = div + UNIT_DIV;
        }
        std::cout << div << std::endl;
    };
    
    void dashed_divider()
    {
        std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
    };
};