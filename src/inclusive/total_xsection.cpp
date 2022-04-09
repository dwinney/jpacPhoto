// Phenomenological expressions for the total cross-sections.
// We use a generic class callable by double total_xsection(double) to select different
// parameterizations or reactions
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "inclusive/total_xsection.hpp"

// ---------------------------------------------------------------------------
// Actualy evaluate the cross-section

double jpacPhoto::total_xsection_PDG::eval(double s)
{
    double result = 0.;
    // if (s < _threshold + 10.*EPS) 
    // {
    //     return 0.;
    // } 
    // else if ((s < _cutoff) && (_sigma.size() > 0))
    // {
    //     result = interp.Eval( pLab(s) );
    // }
    // else
    // {   
        result = PDG_parameterization(s);
    // }
    return result; // output in mb!
};
// ---------------------------------------------------------------------------
// Open up .dat file, import available data and use set up an interpolation
void jpacPhoto::total_xsection_PDG::import_data(std::string datfile)
{
    // Find the correct data file using the top level repo directory
    std::string top_dir;
    char const * env = std::getenv("JPACPHOTO");
    if ( env == NULL || std::string(env) == "" )
    {
        std::cout << "Error! Cannot find top directory. \n";
        std::cout << "Make sure to set environment variable JPACPHOTO!\n";
        exit(0);
    }
    else
    {
        top_dir = std::string(env);
    }
    std::string full_path = top_dir + "/include/inclusive/total_xsection_data/" + datfile;

    std::ifstream infile(full_path);

    if (!infile.is_open())
    {
        std::cout << "ERROR! Cannot open file " << full_path << "!";
        exit(0);
    };

    // Add a zero at exactly threshold
    double plab_thresh = pLab(_threshold + EPS);
    _plab.push_back(plab_thresh);
    _sigma.push_back(0.);
    double old_plab = plab_thresh;
    
    // Import data!
    std::string line;
    while (std::getline(infile, line))
    {   
        if (line.empty()) continue;     // skips empty lines
        std::istringstream is(line);   

        int n;
        double plab, sigma;
        std::string trash; // columns in the file i dont care about

        is >> n;
        is >> plab;
        is >> trash >> trash;
        is >> sigma;

        if (std::abs(old_plab - plab) < 1.E-5 || old_plab > plab) 
        {   
            continue; //no duplicates
        };

        old_plab = plab;

        _plab.push_back(plab);
        _sigma.push_back(sigma);
    };

    interp.SetData(_plab, _sigma);
    return;
};