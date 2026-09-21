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
    // Useful math functions

    std::complex<double> cgamma(std::complex<double> z,int OPT)
    // OPT = 0 for Gamma ; OPT = 1 for log(Gamma)
    {
        std::complex<double> I(0,1);
        std::complex<double> g, infini= 1e308+ 0.*I; // z0,z1
        double x0,q1,q2,x,y,th,th1,th2,g0,gr,gi,gr1,gi1;
        double na,t,x1,y1,sr,si;
        int j,k;
        x1=9e9;
        na=9e9;

        static double a[] = {
            8.333333333333333e-02,
            -2.777777777777778e-03,
            7.936507936507937e-04,
            -5.952380952380952e-04,
            8.417508417508418e-04,
            -1.917526917526918e-03,
            6.410256410256410e-03,
            -2.955065359477124e-02,
            1.796443723688307e-01,
            -1.39243221690590};

        x = real(z); x1 = x;
        y = imag(z); y1 = y;
        if (x > 171) return infini;
        if ((y == 0.0) && (x == (int)x) && (x <= 0.0))
            return infini;
        else if (x < 0.0) {
            x = -x;
            y = -y;
        }
        x0 = x;
        if (x <= 7.0) {
            na = (int)(7.0-x);
            x0 = x+na;
        }
        q1 = sqrt(x0*x0+y*y);
        th = atan(y/x0);
        gr = (x0-0.5)*log(q1)-th*y-x0+0.5*log(2.0*M_PI);
        gi = th*(x0-0.5)+y*log(q1)-y;
        for (k=0;k<10;k++){
            t = pow(q1,-1.0-2.0*k);
            gr += (a[k]*t*cos((2.0*k+1.0)*th));
            gi -= (a[k]*t*sin((2.0*k+1.0)*th));
        }
        if (x <= 7.0) {
            gr1 = 0.0;
            gi1 = 0.0;
            for (j=0;j<na;j++) {
            gr1 += (0.5*log((x+j)*(x+j)+y*y));
            gi1 += atan(y/(x+j));
            }
            gr -= gr1;
            gi -= gi1;
        }
        if (x1 <= 0.0) {
            q1 = sqrt(x*x+y*y);
            th1 = atan(y/x);
            sr = -sin(M_PI*x)*cosh(M_PI*y);
            si = -cos(M_PI*x)*sinh(M_PI*y);
            q2 = sqrt(sr*sr+si*si);
            th2 = atan(si/sr);
            if (sr < 0.0) th2 += M_PI;
            gr = log(M_PI/(q1*q2))-gr;
            gi = -th1-th2-gi;
            x = x1;
            y = y1;
        }
        if (OPT == 0) {
            g0 = exp(gr);
            gr = g0*cos(gi);
            gi = g0*sin(gi);
        }
        g = gr + I*gi;
        return g;
    };

    unsigned int factorial(unsigned int n) 
    {
        if (n == 0)
        return 1;
        return n * factorial(n - 1);
    };

      double legendre(int l, double z)
    {
        switch (l) 
        {
            case 0: return 1.;
            case 1: return z;
            case 2: return 0.5*(3.*z*z - 1.);
            case 3: return 0.5*z*(5.*z*z - 3.);
            case 4: return (35.*z*z*z*z - 30.*z*z + 3.)/8.;
            case 5: return z*(63.*z*z*z*z - 70.*z*z + 15.)/8.;
            default:
            {
                error("legendre - L value " + std::to_string(l) + " not available! Returning 0.", NaN<double>());
            }
        };

        return 0.;
    };

    // --------------------------------------------------------------------------
    // Angular function

    double wigner_leading_coeff(int j, int lam1, int lam2)
    {
        int M = std::max(std::abs(lam1), std::abs(lam2));
        int N = std::min(std::abs(lam1), std::abs(lam2));

        int lambda = std::abs(lam1 - lam2) + lam1 - lam2;

        double result = (double) factorial(2*j);
        result /= sqrt( (double) factorial(j-M));
        result /= sqrt( (double) factorial(j+M));
        result /= sqrt( (double) factorial(j-N));
        result /= sqrt( (double) factorial(j+N));
        result /= pow(2.,  double(j-M));
        result *= pow(-1., double(lambda)/2.);

        return result;
    };

    // USING WIKIPEDIA SIGN CONVENTION
    // theta is in radians
    // lam1 = 2 * lambda and lam2 = 2 * lambda^prime are integers
    double wigner_d_half(int j, int lam1, int lam2, double theta)
    {
        double phase = 1.;
        if ( j % 2 == 0 || (lam1 + lam2) % 2 != 0 )
        {
            error("wigner_d_half - Invalid arguments passed! Returning 0.", NaN<double>());
        };

        if (theta < 0)
        {
            phase *= pow(-1., double(lam1 - lam2)/2.);
            theta *= -1.;
        };

        // If first lam argument is smaller, switch them
        if (std::abs(lam1) < std::abs(lam2))
        {
            int temp = lam1;
            lam1 = lam2;
            lam2 = temp;
            phase *= pow(-1., double((lam1 - lam2)/2.));
        };

        // If first lam is negative, switch them
        if (lam1 < 0)
        {
            lam1 *= -1;
            lam2 *= -1;
            phase *= pow(-1., double((lam1 - lam2)/2.));
        }

        double result = 0.;

        int id = ((lam2 > 0) - (lam2 < 0)) * (j * 100 + lam1 * 10 + std::abs(lam2)); // negative sign refers to negative lam2
        switch (id)
        {
            // spin 1/2 
            case  111: 
            {
                result =  cos(theta / 2.); 
                break;
            };
            case -111: 
            {
                result = -sin(theta / 2.); 
                break;
            };

            // spin 3/2
            case  333:       
            {
                result = cos(theta / 2.) / 2.;
                result *= (1. + cos(theta));
                break;
            }
            case  331:
            {
                result = - sqrt(3.) / 2.;
                result *= sin(theta / 2.);
                result *= 1. + cos(theta);
                break;
            }
            case -331:
            {
                result = sqrt(3.) / 2.;
                result *= cos(theta / 2.);
                result *= 1. - cos(theta);
                break;
            }
            case -333:
            {
                result = - sin(theta / 2.) / 2.;
                result *= 1. - cos(theta);
                break;
            }
            case  311:
            {
                result = 1. / 2.;
                result *= 3. * cos(theta) - 1.;
                result *= cos(theta / 2.);
                break;
            }
            case -311:
            {
                result = -1. / 2.;
                result *= 3. * cos(theta) + 1.;
                result *= sin(theta / 2.);
                break;
            }

            // Spin- 5/2
            case  533:
            {
                result = -1. / 4.;
                result *= cos(theta / 2.);
                result *= (1. + cos(theta)) * (3. - 5. * cos(theta));
                break;
            }
            case  531:
            {
                result = sqrt(2.) / 4.;
                result *= sin(theta / 2.);
                result *= (1. + cos(theta)) * (1. - 5. * cos(theta));
                break;
            }
            case -531:
            {
                result =  sqrt(2.) / 4.;
                result *= cos(theta / 2.);
                result *= (1. - cos(theta)) * (1. + 5. * cos(theta));
                break;
            }
            case -533:
            {
                result = -1. / 4.;
                result *= sin(theta / 2.);
                result *= (1. - cos(theta)) * (3. + 5. * cos(theta));
                break;
            }
            case  511:
            {
                result = -1. / 2.;
                result *= cos(theta / 2.);
                result *= (1. + 2. * cos(theta) - 5. * cos(theta)*cos(theta));
                break;
            }
            case -511:
            {
                result = 1. / 2.;
                result *= sin(theta / 2.);
                result *= (1. - 2. * cos(theta) - 5. * cos(theta)*cos(theta));
                break;
            }

            default: return NaN<double>();
        };

        return phase * result;
    };

    double wigner_d_int(int j, int lam1, int lam2, double theta)
    {

        double phase = 1.;

        if (theta < 0)
        {
            phase *= pow(-1., double(lam1 - lam2));
            theta *= -1.;
        };

        // If first lam argument is smaller, switch them
        if (std::abs(lam1) < std::abs(lam2))
        {
            int temp = lam1;
            lam1 = lam2;
            lam2 = temp;
            phase *= pow(-1., double(lam1 - lam2));
        };

        // If first lam is negative, smitch them
        if (lam1 < 0)
        {
            lam1 *= -1;
            lam2 *= -1;
            phase *= pow(-1., double(lam1 - lam2));
        }

        // Output
        double result = 0.;
        int id = ((lam2 >= 0) - (lam2 < 0)) * (j * 100 + lam1 * 10 + std::abs(lam2)); // negative sign refers to negative lam2
        switch (id)
        {   
            // Spin 1
            case  111:
            {
                result = (1. + cos(theta)) / 2.;
                break;
            }
            case  110:
            {
                result = - sin(theta) / sqrt(2.);
                break;
            }
            case -111:
            {
                result = (1. - cos(theta)) / 2.;
                break;
            }
            case  100:
            {
                result = cos(theta);
                break;
            }

            default: return 0.;
        }

        return phase * result;
    };

    // Wigner functions butn ow in terms of the costheta, this allows 
    // an analytic continuation to complex angular polynomials
    complex wigner_d_int_cos(int j, int lam1, int lam2, complex cosine)
    {
        // Careful because this loses the +- phase of the sintheta. 
        complex sine = sqrt(XR - cosine * cosine);

        double phase = 1.;
        // If first lam argument is smaller, switch them
        if (std::abs(lam1) < std::abs(lam2))
        {
            int temp = lam1;
            lam1 = lam2;
            lam2 = temp;

            phase *= pow(-1., double(lam1 - lam2));
        };

        // If first lam is negative, smitch them
        if (lam1 < 0)
        {
            lam1 *= -1;
            lam2 *= -1;

            phase *= pow(-1., double(lam1 - lam2));
        }

        complex result = 0.;
        int id = ((lam2 >= 0) - (lam2 < 0)) * (j * 100 + lam1 * 10 + std::abs(lam2)); // negative sign refers to negative lam2
        switch (id)
        {   
            // Spin 1
            case  111:
            {
                result = (1. + cosine) / 2.;
                break;
            }
            case  110:
            {
                result = - sine / sqrt(2.);
                break;
            }
            case -111:
            {
                result = (1. - cosine) / 2.;
                break;
            }
            case  100:
            {
                result = cosine;
                break;
            }
            default: return NaN<complex>();
        }

        return phase * result;
    };

     complex wigner_d_half_cos(int j, int lam1, int lam2, complex cosine)
    {
        // Careful because this loses the +- phase of the sintheta. 
        complex sine = sqrt(XR - cosine * cosine);

        // Also need the half-angle factors
        complex sinhalf =  sqrt((XR - cosine) / 2.);
        complex coshalf =  sqrt((XR + cosine) / 2.);

        double phase = 1.;
        if ( j % 2 == 0 || (lam1 + lam2) % 2 != 0 )
        {
            error("wigner_d_half - Invalid arguments passed! Returning 0.", NaN<complex>());
        };

        // If first lam argument is smaller, switch them
        if (std::abs(lam1) < std::abs(lam2))
        {
            int temp = lam1;
            lam1 = lam2;
            lam2 = temp;
            phase *= pow(-1., double((lam1 - lam2)/2.));
        };

        // If first lam is negative, switch them
        if (lam1 < 0)
        {
            lam1 *= -1;
            lam2 *= -1;
            phase *= pow(-1., double((lam1 - lam2)/2.));
        }

        
        int id = ((lam2 > 0) - (lam2 < 0)) * (j * 100 + lam1 * 10 + std::abs(lam2)); // negative sign refers to negative lam2
        complex result = 0.;
        switch (id)
        {
            // spin 1/2 
            case  111: 
            {
                result =  coshalf; 
                break;
            };
            case -111: 
            {
                result = -sinhalf; 
                break;
            };

            // spin 3/2
            case  333:       
            {
                result = coshalf / 2.;
                result *= (1. + cosine);
                break;
            }
            case  331:
            {
                result = - sqrt(3.) / 2.;
                result *= coshalf;
                result *= 1. + cosine;
                break;
            }
            case -331:
            {
                result = sqrt(3.) / 2.;
                result *= coshalf;
                result *= 1. - cosine;
                break;
            }
            case -333:
            {
                result = - sinhalf / 2.;
                result *= 1. - cosine;
                break;
            }
            case  311:
            {
                result = 1. / 2.;
                result *= coshalf;
                result *= 3. * cosine - 1.;
                break;
            }
            case -311:
            {
                result = -1. / 2.;
                result *= sinhalf;
                result *= 3. * cosine + 1.;
                break;
            }

            // Spin- 5/2
            case  533:
            {
                result = -1. / 4.;
                result *= coshalf;
                result *= (1. + cosine) * (3. - 5. * cosine);
                break;
            }
            case  531:
            {
                result = sqrt(2.) / 4.;
                result *= sinhalf;
                result *= (1. + cosine) * (1. - 5. * cosine);
                break;
            }
            case -531:
            {
                result =  sqrt(2.) / 4.;
                result *= coshalf;
                result *= (1. - cosine) * (1. + 5. * cosine);
                break;
            }
            case -533:
            {
                result = -1. / 4.;
                result *= sinhalf;
                result *= (1. - cosine) * (3. + 5. * cosine);
                break;
            }
            case  511:
            {
                result = -1. / 2.;
                result *= coshalf;
                result *= (1. + 2. * cosine - 5. * cosine*cosine);
                break;
            }
            case -511:
            {
                result = 1. / 2.;
                result *= sinhalf;
                result *= (1. - 2. * cosine - 5. * cosine*cosine);
                break;
            }
            default: return NaN<complex>();
        };
        return phase * result;
    };

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
    // kallen Triangle function

    // If any of them are complex, return complex
    complex kallen(complex z, double a, double b) { return kallen<complex>(z, XR*a, XR*b); };
    complex kallen(double a, complex z, double b) { return kallen<complex>(XR*a, z, XR*b); };
    complex kallen(double a, double b, complex z) { return kallen<complex>(XR*a, XR*b, z); };

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


    bool is_real(complex a, double tol)
    {
        return is_zero(imag(a), tol);
    };

    bool is_imaginary(complex a, double tol )
    {
        return is_zero(real(a), tol);
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

    // ---------------------------------------------------------------------------
    // Importing data sets we'll need to be able to find the main directory from the 
    // top level one. Thus we need to be able to access the appropriate environment variable
    std::string main_dir()
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
    std::string desktop()
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

    std::string to_string(double d, uint precision)
    {
        std::stringstream ss;
        ss << std::setprecision(precision) << d;
        return ss.str();
    };

    // Print a string centered on the terminal 
    void centered(std::string words)
    {
        int x = words.length();
        int gap_width = (TEXT_WIDTH - x)/2;
        std::cout << std::left << std::setw(gap_width) << "" << std::setw(x) << words << std::setw(gap_width) << "" << std::endl;
    };

    // ---------------------------------------------------------------------------
    // Print functions evaluated on a grid to an ascii file

    // Take in a function and print an ascii file of values on a grid
    void print_to_file(std::array<double, 2> bounds, std::function<double(double)> F, std::string file)
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

    void print_to_file(std::array<double, 2> bounds, std::vector<std::function<double(double)>> Fs, std::string file)
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

    void print_to_file(std::array<double,2> boundsx, std::array<double,2> boundsy, std::function<double(double,double)> F, std::string file)
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

    // ---------------------------------------------------------------------------
    // Element-wise operations on data vectors

    std::vector<double> multiply_elementwise(std::vector<double> in1, std::vector<double> in2)
    {
        if (in1.size() != in2.size()) warning("multiply_elementwise() - Input vectors not the same size!");

        std::vector<double> out;
        for (uint i = 0; i < in1.size(); i++)
        {
            out.push_back(in1[i]*in2[i]);
        }
        return out;
    };

    std::vector<double> square_elementwise(std::vector<double> in){ return multiply_elementwise(in, in); };

    std::vector<double> real(std::vector<complex> vx)
    {
        std::vector<double> out;
        for (auto x : vx) out.push_back(real(x));
        return out;
    };
    
    std::vector<double> imag(std::vector<complex> vx)
    {
        std::vector<double> out;
        for (auto x : vx) out.push_back(imag(x));
        return out;
    };

};