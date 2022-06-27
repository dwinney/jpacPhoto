// Phenomenological expressions for the total cross-sections from the PDG
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "inclusive/sigma_PDG.hpp"

// ---------------------------------------------------------------------------
// METHODS FOR THE PDG_parameterization
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Import available data and use set up an interpolation
// The .dat files for low energy in PDG format:
// POINT_NUMBER PLAB(GEV/C) PLAB_MIN PLAB_MAX SIG(MB) STA_ERR+ STA_ERR- SY_ER+(PCT) SY_ER-(PCT) REFERENCE FLAG
void jpacPhoto::PDG_parameterization::import_data(std::string datfile)
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
    double plab_thresh = pLab(_sth + EPS);
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

    _interp.SetData(_plab, _sigma);
    return;
};
