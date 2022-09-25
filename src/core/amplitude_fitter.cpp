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

    // Grab the ith data set 
    data_set ith_data_set = _integrated_data[i];

    // Individual quantities
    std::vector<double> s     = ith_data_set._s;
    std::vector<double> sigma = ith_data_set._sigma;
    std::vector<double> error = ith_data_set._error;

    // Calculate chi2
    double chi2 = 0;
    for (int n = 0; n < s.size(); n++)
    {   
        double sigma_th = _amplitude->integrated_xsection(s[n]);
        double sigma_ex = sigma[n];
        double err      = error[n];

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
    std::vector<double> params;
    for (int n = 0; n < _pars.size(); n++)
    {
        params.push_back(par[n]);
    };

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
    double s                  = ith_data_set._fixed_s;
    std::vector<double> t     = ith_data_set._t;
    std::vector<double> sigma = ith_data_set._dsigma;
    std::vector<double> error = ith_data_set._derror;

    // Calculate chi2
    double chi2 = 0;
    for (int n = 0; n < t.size(); n++)
    {   
        double sigma_th = _amplitude->differential_xsection(s, t[n]);
        double sigma_ex = sigma[n];
        double err      = error[n];

        chi2 += pow( (sigma_th - sigma_ex) / err, 2.);
    }; 

    return chi2;
};

//-----------------------------------------------------------------------
// Fitting routines!
double jpacPhoto::amplitude_fitter::fit_integrated(std::vector<double> starting_guess)
{
    if (starting_guess.size() != _pars.size()) 
    {
        std::cout << "amplitude_fitter::fit_integrated() : Error! Size of starting vector not the expected size of parameters (" << std::to_string(_pars.size()) <<"). \n";
        std::cout << "Returning without fit..." << std::endl;
        return -1.;
    };

    // MINUIT object from root
    ROOT::Math::Minimizer* minuit = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
    minuit->SetMaxFunctionCalls(1E7);
    minuit->SetTolerance(1.E-6);
    minuit->SetPrintLevel(_nError);

    ROOT::Math::Functor fcn(this, &amplitude_fitter::chi2_integrated, _pars.size());
    minuit->SetFunction(fcn);

    for (int a = 0; a < _pars.size(); a++)
    {   
        minuit->SetVariable(a, _pars[a]._label, starting_guess[a], 0.1);

        if (_pars[a]._custom_limits)
        {
            minuit->SetVariableLimits(a, _pars[a]._lower_limit, _pars[a]._upper_limit);
        }
        if (_pars[a]._fixed)
        {
            minuit->FixVariable(a);
        }
    };

    new_line(); integrated_data_info();
    new_line(); divider(); 
    std::cout << std::left << "Fitting " + std::to_string(minuit->NFree()) << " parameters with starting values:\n";   
    new_line(); variable_info(starting_guess, 1);
    new_line(); divider(); new_line();

    std::cout << "Beginning fit...";
    minuit->Minimize();
    std::cout << "done! \n";
    new_line();

    // Fit results
    double chi2    = minuit->MinValue();
    double chi2dof = chi2 / (_N_int - minuit->NFree());
    std::vector<double> best_params = convert(minuit->X());

    divider();
    std::cout << std::left << std::setw(10) << "chi2 = " << std::setw(15) << chi2;
    std::cout << std::left << std::setw(10) << "chi2/dof = "   << std::setw(15) << chi2dof << std::endl;
    new_line();
    variable_info(best_params, 0);
    divider();
    
    // At the end update the amplitude parameters to include the fit results
    _amplitude->set_params(best_params);

    return 0;
};