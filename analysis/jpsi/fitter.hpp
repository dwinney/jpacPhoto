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

#include "jpsi/fitter.hpp"
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
            // String to display what type of data each set contains
            static std::string data_type(int i)
            {
                switch (i)
                {
                    case kIntegrated:      return "Integrated (Eg)";
                    case kDifferential:    return "Differential (Eg vs -t)";
                    case kDifferential_tp: return "Differential (EG vs -t')";
                    default: return "ERROR";  
                };
                return "ERROR";
            };

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
        };
    };
};

#endif