// Class which takes in an amplitude and some data and produces a fit based on chi2 minimization with MINUIT
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitude_fitter.hpp"

// ---------------------------------------------------------------------------
// Calculate the chi2 by summing over data sets

// Cummulative chi2 fron integrated data supplied
double jpacPhoto::amplitude_fitter::chi2_integrated(const double *par)
{
    // Clear the last evaluation of chi2's
    _integrated_chi2s.clear();

    // convert double *par to vector
    std::vector<double> params = convert(par);

    // Iterate over each set data set
    double chi2 = 0.;
    for (int i = 0; i < _integrated_data.size(); i++)
    {
        // Save the contributions from each individual data set
        double chi2_i =  chi2_integrated_i(i, params);
        _integrated_chi2s.push_back(chi2_i);
        chi2 += chi2_i;
    };

    // Return the cummulative total
    return chi2;
};

// Calculate chi2 from the i-th data set added 
double jpacPhoto::amplitude_fitter::chi2_integrated_i(int i, std::vector<double> pars)
{   
    // Pass the parameters to the amplitude
    _amplitude->set_params(pars);

    // Individual quantities
    std::vector<double> vs     = _integrated_data[i]._s;
    std::vector<double> vsigma = _integrated_data[i]._sigma;
    std::vector<double> verror = _integrated_data[i]._error;

    // Calculate chi2
    double chi2 = 0;
    for (int n = 0; n < vs.size(); n++)
    {   
        double s;
        
        (_integrated_data[i]._useEgamma) ? (s = pow(W_cm(vs[n]), 2.)) : (s = vs[n]);

        double sigma_th = _amplitude->integrated_xsection(s);
        double sigma_ex = vsigma[n];
        double err      = verror[n];

        chi2 += pow( (sigma_th - sigma_ex) / err, 2.);
    }; 

    return chi2;
};

// Cummulative chi2 fron differential data supplied
double jpacPhoto::amplitude_fitter::chi2_differential(const double *par)
{
    // Clear the last evaluation of chi2's
    _differential_chi2s.clear();

    // convert double *par to vector
    std::vector<double> params = convert(par);


    // Iterate over each set data set
    double chi2 = 0.;
    for (int i = 0; i < _differential_data.size(); i++)
    {
        // Save the contributions from each individual data set
        double chi2_i =  chi2_differential_i(i, params);
        _differential_chi2s.push_back(chi2_i);

        chi2 += chi2_i;
    };

    // Return the cummulative total
    return chi2;
};

// Calculate chi2 from the i-th data set added 
double jpacPhoto::amplitude_fitter::chi2_differential_i(int i, std::vector<double> pars)
{   
    // Pass the parameters to the amplitude
    _amplitude->set_params(pars);

    // Grab the ith data set 
    data_set ith_data_set = _differential_data[i];

    // Individual quantities
    double x                  = ith_data_set._fixed_s;
    std::vector<double> vt    = ith_data_set._t;
    std::vector<double> sigma = ith_data_set._dsigma;
    std::vector<double> error = ith_data_set._derror;

    // Calculate chi2
    double chi2 = 0;
    for (int n = 0; n < vt.size(); n++)
    {   
        double s, t;

        (ith_data_set._useEgamma) ? (s = pow(W_cm(x), 2.)) : (s = x);

        // Check if momentum transfer is negative (i.e. is what was saves -t or t ?)
        double xt = vt[n]; 
        if (xt > 0) xt *= -1.;

        // Then convert from t' to t if necessary
        (ith_data_set._useTPrime) ? (t = xt + _amplitude->_kinematics->t_man(s, 0.)) : (t = xt);

        double sigma_th = _amplitude->differential_xsection(s, t);
        double sigma_ex = sigma[n];
        double err      = error[n];

        chi2 += pow( (sigma_th - sigma_ex) / err, 2.);
    }; 

    return chi2;
};

//-----------------------------------------------------------------------
// These are the actual fitting routines which call instances of MINUIT and do the fit

double jpacPhoto::amplitude_fitter::do_fit(std::vector<double> starting_guess)
{
    if (starting_guess.size() != _pars.size()) 
    {
        std::cout << "amplitude_fitter::fit_differential() : Error! Size of starting vector not the expected size of parameters (" << std::to_string(_pars.size()) <<"). \n";
        std::cout << "Returning without fit..." << std::endl;
        return -1.;
    };

    set_up(starting_guess);

    new_line(); data_info();
    new_line(); divider(); 
    new_line(); variable_info(starting_guess, 0);
    new_line(); divider(); new_line();

    std::cout << "Beginning fit..." << std::flush;    
    _minuit->Minimize();
    std::cout << "Done! \n";

    print_results();

    return 0;
};

//-----------------------------------------------------------------------
// Methods to print out relevant info to command line

void jpacPhoto::amplitude_fitter::data_info()
{
    std::cout << std::left << "Fitting amplitude (\"" << _amplitude->get_id() << "\") to " << _N << " data points: \n";
    new_line();
    std::cout << std::left << std::setw(30) << "DATA SET"            << std::setw(20) << "OBSERVABLE"      << std::setw(10) << "POINTS"  << std::endl;
    std::cout << std::left << std::setw(30) << "----------------"    << std::setw(20) << "--------------"  << std::setw(10) << "-------" << std::endl;
    for (int k = 0; k < _integrated_data.size(); k++)
    {
        std::cout << std::left << std::setw(30) << _integrated_data[k]._id  << std::setw(20)  << "integrated xs"   << std::setw(10) <<  _integrated_data[k]._s.size()  << std::endl;  
    };
    for (int k = 0; k < _differential_data.size(); k++)
    {   
        // std::stringstream ss;
        // ss << std::setprecision(3) << "s = "    << _differential_data[k]._fixed_s << " GeV2";
        std::cout << std::left << std::setw(30) << _differential_data[k]._id << std::setw(20) << "differential xs" << std::setw(10) <<  _differential_data[k]._t.size()/*  << std::setw(10) << ss.str() */ << std::endl;  
    };
};

// Print out a little table of the current status of parameters
void jpacPhoto::amplitude_fitter::variable_info(std::vector<double> starting_guess, bool opt)
{  
    std::cout << std::setprecision(10);

    std::string column_3;
    (opt) ? (column_3 = "FIT VALUE") : (column_3 = "START VALUE");

    if (!opt)  std::cout << std::left << "Fitting " + std::to_string(_minuit->NFree()) << " (out of " << std::to_string(_minuit->NDim()) << ") parameters:\n" << std::endl; 
    std::cout << std::left << std::setw(10) << "N" << std::setw(20) << "PARAMETER" << std::setw(10) << column_3 << std::endl;
    std::cout << std::left << std::setw(10) << "-----" << std::setw(20) << "----------" << std::setw(10) << "------------"<< std::endl;

    for (int i = 0; i < _pars.size(); i++)
    {
        std::string extra = "";
        if (_pars[i]._custom_limits && !opt)
        {   
            std::stringstream ss;
            ss << std::setprecision(5) << "[" << _pars[i]._lower_limit << ", " << _pars[i]._upper_limit << "]";
            extra = ss.str();
        };
        if (_pars[i]._fixed && !opt) extra = "FIXED";
        std::cout << std::left << std::setw(10) << i << std::setw(20) << _pars[i]._label << std::setw(20) << starting_guess[i] << std::setw(10) << extra << std::endl;
    };
};

void jpacPhoto::amplitude_fitter::set_up(std::vector<double> starting_guess)
{
    _minuit->Clear();
    _minuit->SetTolerance(1.E-6);
    _minuit->SetPrintLevel(_printLevel);
    _minuit->SetMaxFunctionCalls(_maxCalls);
    
    for (int a = 0; a < _pars.size(); a++)
    {   
        _minuit->SetVariable(a, _pars[a]._label, starting_guess[a], _pars[a]._step_size);

        if (_pars[a]._custom_limits)
        {
            _minuit->SetVariableLimits(a, _pars[a]._lower_limit, _pars[a]._upper_limit);
        }
        if (_pars[a]._fixed)
        {
            _minuit->FixVariable(a);
        }
    };

    fcn = ROOT::Math::Functor(this, &amplitude_fitter::chi2, _pars.size());
    _minuit->SetFunction(fcn);
};

void jpacPhoto::amplitude_fitter::print_results()
{
    // Fit results
    int dof = _N - _minuit->NFree();

    new_line();
    double chi2    = _minuit->MinValue();
    double chi2dof = chi2 / double(dof) ;
    std::vector<double> best_params = convert(_minuit->X());

    divider();
    std::cout << std::left << std::setw(10) << "chi2 = "      << std::setw(15) << chi2;
    std::cout << std::left << std::setw(10) << "chi2/dof = "  << std::setw(15) << chi2dof << std::endl;
    new_line();
    variable_info(best_params, 1);
    divider();
    new_line();
    
    // At the end update the amplitude parameters to include the fit results
    _amplitude->set_params(best_params);
}