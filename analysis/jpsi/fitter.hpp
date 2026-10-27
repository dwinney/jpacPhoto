// Interface for amplitude fitting to J/psi photoproduction data
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Instituto de Ciencias Nucleares (ICN),
//               Universidad Nacional Autonoma de México (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#ifndef JPSI_FITTER_HPP
#define JPSI_FITTER_HPP

#include "amplitude.hpp"
#include "utilities.hpp"

namespace jpacPhoto
{
    namespace jpsi
    {
        // Return 
        std::string data_dir(){ return main_dir()+"/analysis/jpsi/data"; };
    
        // data_set stypes
        const uint kIntegrated      = 0; // sigma     vs Egam
        const uint kDifferential    = 1; // dsigma/dt vs -t
        const uint kDifferential_tp = 2; // dsigma/dt vs -t'

        // Actual implementation of the fitter
        struct unpolarized
        {
            // ------------------------------------------------------------------------------
            // Data handling 

            // String to display what type of data each set contains
            static std::string data_type(int i)
            {
                switch (i)
                {
                    case kIntegrated:      return "Integrated (Eg)";
                    case kDifferential:    return "Differential (Eg vs -t)";
                    case kDifferential_tp: return "Differential (Eg vs -t')";
                    default: return "ERROR";  
                };
                return "ERROR";
            };

            // ------------------------------------------------------------------------------
            // Chi2 calculation

            // Function being minimized, sum of chi2's of each data set
            static double fcn(const std::vector<data_set> &data_sets, amplitude amp)
            {
                double chi2 = 0;
                for (auto data : data_sets)
                {
                    switch (data._type)
                    {
                        case kIntegrated:      { chi2 += chi2_integrated(data, amp);          break; }
                        case kDifferential:    { chi2 += chi2_differential(data, amp, false); break; }
                        case kDifferential_tp: { chi2 += chi2_differential(data, amp, true);  break; }
                        default: { chi2 += NaN<double>(); };
                    };
                };
                return chi2;
            };

            // Calculate chi2 from an integrated cross section data set
            static double chi2_integrated(data_set &data, amplitude amp)
            {
                double chi2 = 0;
                for (int i = 0; i < data._N; i++)
                {
                    double s = s_cm(data._x[i]);
                    double sig_th = amp->integrated_xsection(s);
                    double sig_ex = data._z[i];
                    double error  = data._zerr[0][i];
                    chi2 += pow((sig_th - sig_ex) / error, 2);
                };
                return chi2;
            };

            // Calculate chi2 from a differential data set
            static double chi2_differential(data_set &data, amplitude amp, bool use_tp)
            {
                double chi2 = 0;
                for (int i = 0; i < data._N; i++)
                {
                    double s  = s_cm(data._y[i]);
                    double t  = (use_tp) ? -data._x[i]+amp->get_kinematics()->t_min(s) : -data._x[i];

                    double sig_th = amp->differential_xsection(s, t);
                    double sig_ex = data._z[i];
                    double error  = data._zerr[0][i];
                    chi2 += pow((sig_th - sig_ex) / error, 2);
                };
                return chi2;
            };

            // ------------------------------------------------------------------------------
            // Export 

            static void export_parameters(std::vector<double> pars, std::string filename, std::string description = "")
            {
                using std::endl, std::setw, std::left;

                std::ofstream out;
                out.open(filename);
                int  precision = 12, spacing = precision+10;
                out << std::setprecision(precision);
                out << std::showpos;

                // Preamble info
                if (!description.empty()) out << left << "# " + description << endl;
                out << left << "#"+std::string(3*spacing-10, '-') << endl;
                out << left << setw(12) << "# n" << setw(spacing) << "    f_n" << setw(spacing) << "    a_n" << endl;
                out << left << "#"+std::string(3*spacing-10, '-') << endl;
                out << left << "#                        J = 1/2" << endl;
                out << left << "#"+std::string(3*spacing-10, '-') << endl;

                for (int n = 0; n < 3; n++)
                {
                    out << left << setw(12) << "  "+to_string(n);
                    if (n != 2) out << left << setw(spacing) << pars[n];
                    else        out << left << setw(spacing) << "";
                    out << left << setw(spacing) << pars[2+n] << endl;
                };
                out << left << "#"+std::string(3*spacing-10, '-') << endl;
                out << left << "#                        J = 3/2" << endl;
                out << left << "#"+std::string(3*spacing-10, '-') << endl;

                for (int n = 0; n < 6; n++)
                {
                    out << left << setw(12) << "  "+to_string(n);
                    out << left << setw(spacing) << pars[5+n];
                    out << left << setw(spacing) << pars[11+n] << endl;
                };
                out.close();
            };

            static std::vector<double> import_parameters(std::string filename)
            {
                std::vector<double> fs, as;

                std::ifstream infile(filename);
                std::string line;

                bool already_skipped = false;
                while (std::getline(infile, line))
                {   
                    if (line.empty())        continue; // skips empty lines
                    if (line.front() == '#') continue; // Skip comment lines 
                    std::istringstream is(line);  
                    
                    double n, f, a;
                    // Dont care about first two columns
                    is >> n;
                    if (n == 2 && !already_skipped)
                    {
                      is >> a;
                      as.push_back(a);
                      already_skipped = true;
                    }
                    else
                    {
                      is >> f >> a;
                      fs.push_back(f); as.push_back(a);
                    } 
                    continue;
                };

                std::vector<double> pars;
                // J = 1/2
                for (int i = 0; i < 2; i++) pars.push_back(fs[i]);
                for (int i = 0; i < 3; i++) pars.push_back(as[i]);
                // J = 3/2
                for (int i = 0; i < 6; i++) pars.push_back(fs[2+i]);
                for (int i = 0; i < 6; i++) pars.push_back(as[3+i]);

                return pars;
            };
        };
    };
};

#endif