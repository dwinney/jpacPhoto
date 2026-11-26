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

    // Kinematic function for 2->2 scattering (see eq. 5.23 in Byckling & Kajantie)
    double G(double x, double y, double z, double u, double v, double w)
    {
        return  x*x*y + x*y*y + z*z*u + z*u*u + v*v*w + v*w*w 
              + x*z*w + x*u*v + y*z*w + y*u*w + y*z*v - y*z*w
              - x*y*(z + u + v + w) - z*u*(x + y + v + w) - v*w*(x + y + z + u);
    };

    // ---------------------------------------------------------------------------
    // Importing data sets we'll need to be able to find the main directory from the 
    // top level one. Thus we need to be able to access the appropriate environment variable
    std::string main_dir()
    {
       // Find the correct data file using the top level repo directory
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
        char const * env = std::getenv("DESKTOP");
        if ( env == NULL || std::string(env) == "" )
        {
            return error("desktop(): Cannot find environment variable DESKTOP!", "");
        }
        return std::string(env);  
    };

    // Sometimes we'll want to work on an analysis outside the main JPACPHOTO directory
    // this looks for the variable JPACPHOTO_WORKING_DIR
    std::string working_dir()
    {
       // Find the correct data file using the top level repo directory
        char const * env = std::getenv("JPACPHOTO_WORKING_DIR");
        if ( env == NULL || std::string(env) == "" )
        {
            return error("desktop(): Cannot find environment variable JPACPHOTO_WORKING_DIR!", "");
        }
        return std::string(env);  
    };
};