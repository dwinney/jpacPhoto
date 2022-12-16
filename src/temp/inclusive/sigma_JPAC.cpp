// Phenomenological expressions for the total cross-sections using SAID PW's at low energies
// and Regge poles at high energies
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "sigma_JPAC.hpp"

// ---------------------------------------------------------------------------
// METHODS FOR THE JPAC_parameterization
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// For a given isospin projection, assemble the 16 partial wave amplitudes
void jpacPhoto::JPAC_parameterization::initialize_PWAs()
{
    // Get 8 waves per parity and per isospin 
    for (int L = 0; L <= _Lmax ; L++)
    {
        _pw_1p.push_back( new SAID_PWA(L, 1, 2*L+1) );
        _pw_1m.push_back( new SAID_PWA(L, 1, 2*L-1) );
        _pw_3p.push_back( new SAID_PWA(L, 3, 2*L+1) );
        _pw_3m.push_back( new SAID_PWA(L, 3, 2*L-1) );
    };

    return;
};

// ---------------------------------------------------------------------------
void jpacPhoto::JPAC_parameterization::update_amplitudes(double s)
{
    _CpL.clear(); _CmL.clear();

    double W    = sqrt(s);
    double E    = (s + M2_PROTON - M2_PION) / (2.*W);
    double Elab = (s - M2_PROTON - M2_PION)/(2.*M_PROTON);
    double qcm  = sqrt(Kallen(s, M2_PROTON, M2_PION)) / (2.*W);

    // Calculate the s-channel isospin amplitudes
    for (int L = 0; L <= _Lmax; L++)
    {
        double f1 = 0., f3 = 0., g1 = 0., g3 = 0.;
        double fp, fm, gp, gm;
        double f1p, f1m, f2p, f2m;
        double Ap, Am, Bp, Bm, Cp, Cm;

        // Assemble amplitudes from PWAs
        f1 = ((L+1)*_pw_1p[L]->imaginary_part(s) + (L)*_pw_1m[L]->imaginary_part(s));
        f3 = ((L+1)*_pw_3p[L]->imaginary_part(s) + (L)*_pw_3m[L]->imaginary_part(s));
        g1 = coeffP(L) * (_pw_1p[L]->imaginary_part(s) - _pw_1m[L]->imaginary_part(s));
        g3 = coeffP(L) * (_pw_3p[L]->imaginary_part(s) - _pw_3m[L]->imaginary_part(s));

        f1 /= qcm; f3 /= qcm; g1 /= qcm; g3 /= qcm;

        // Use these to calculate the t-channel amplitudes
        fp = (f1 + 2.*f3) / 3.;
        gp = (g1 + 2.*g3) / 3.;
        fm = (f1 - f3)    / 3.;
        gm = (g1 - g3)    / 3.;


        // Amplitudes f1 and f2 with t-channel isospin
        f1p = fp + gp;
        f2p = -gp;
        f1m = fm + gm;
        f2m = -gm;

        // Invariant amplitudes
        Ap = (W+M_PROTON)/(E+M_PROTON)*f1p - (W-M_PROTON)/(E-M_PROTON)*f2p;
        Am = (W+M_PROTON)/(E+M_PROTON)*f1m - (W-M_PROTON)/(E-M_PROTON)*f2m;

        Bp = f1p/(E+M_PROTON) + f2p/(E-M_PROTON);
        Bm = f1m/(E+M_PROTON) + f2m/(E-M_PROTON);

        Cp = 4.*M_PI * (Ap + Elab * Bp);
        Cm = 4.*M_PI * (Am + Elab * Bm);

        _CpL.push_back(Cp);
        _CmL.push_back(Cm);
    };
};

// ---------------------------------------------------------------------------
// Print out the appropriate filename for the partial wave whose quantum numbers
// are stored
std::string jpacPhoto::JPAC_parameterization::SAID_PWA::get_filename()
{
    std::string filename = "/include/inclusive/total_xsection_data/SAID_PW/SAID_PiN_";
    switch (_L)
    {
        case 0: filename += "S"; break;
        case 1: filename += "P"; break;
        case 2: filename += "D"; break;
        case 3: filename += "F"; break;
        case 4: filename += "G"; break;
        case 5: filename += "H"; break;
        case 6: filename += "I"; break;
        case 7: filename += "J"; break;
    }

    filename += std::to_string(_I);
    filename += std::to_string(_J);
    filename += ".txt";

    return filename;
};

// ---------------------------------------------------------------------------
// Import available data and use set up an interpolation
void jpacPhoto::JPAC_parameterization::SAID_PWA::import_data()
{
    std::string filename = get_filename();

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
    std::string full_path = top_dir + filename;
    std::ifstream infile(full_path);

    if (!infile.is_open())
    {
        std::cout << "ERROR! Cannot open file " << full_path << "!";
        exit(0);
    };

    // Import data!
    std::string line;
    while (std::getline(infile, line))
    {   
        if (line.empty()) continue;     // skips empty lines
        std::istringstream is(line);   

        double plab, imA, reA;
        double trash; // columns in the file i dont care about

        is >> plab;
        is >> trash >> trash >> trash >> trash;
        is >> reA >> imA >> trash >> trash;

        _plab.push_back(plab * 1.E-3);
        _realPart.push_back(reA);
        _imagPart.push_back(imA);
    };

    _interpReal.SetData(_plab, _realPart);
    _interpImag.SetData(_plab, _imagPart);

    return;
};
