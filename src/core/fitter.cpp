// Class which allows an amplitude to be fit to data based on chi2 minimization
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#include "fitter.hpp"
#include "data_set.hpp"

namespace jpacPhoto
{
    // -----------------------------------------------------------------------
    // Methods for managing data

    void fitter::add_data(data_set data)
    {
        switch (data._type)
        {
            case integrated_data:   _int_data.push_back(data); break;
            case differential_data: _dif_data.push_back(data); break;
            default:
            {
                warning("fitter::add_data", "data_set " + data._id + " of unsupported type!");
                return;
            }
        }

        _N += data._N;
    };

    void fitter::clear_data()
    {
        _N = 0;
        _int_data.clear();
        _dif_data.clear();
    };

    // -----------------------------------------------------------------------
    // Set limits, labels, and fix parameters

    void fitter::reset_parameters()
    {
        _pars.clear();
        // populate parameters vector of appropriate size
        for (int i = 0; i < _amplitude->N_pars(); i++)
        {
            _pars.push_back(i);
        };
        set_parameter_labels(_amplitude->parameter_labels());
    };

    // Give each parameter a label beyond their default par[i] name
    void fitter::set_parameter_labels(std::vector<std::string> labels)
    {
        if (labels.size() != _pars.size())
        {
            warning("fitter::set_parameter_labels", "Labels vector does not match number of parameters!");
            return;
        }

        for (int i = 0; i < _pars.size(); i++)
        {
            _pars[i]._label = labels[i];
        };
    };

    // Given a parameter label, find the corresponding index
    int fitter::find_parameter(std::string label)
    {
        for (auto par : _pars)
        {
            if (par._label == label) return par._i;
        };

        return error("fitter::find_parameter", "Cannot find parameter labeled " + label + "!", NaN<int>());
    };

    // Set limits and/or a custom stepsize
    void fitter::set_parameter_limits(int i, std::array<double,2> bounds, double step)
    {
        _pars[i]._custom_limits = true;
        _pars[i]._lower         = bounds[0];
        _pars[i]._upper         = bounds[1];
        _pars[i]._step          = step;
    };

    void fitter::set_parameter_limits(std::string label, std::array<double,2> bounds, double step)
    {
        return set_parameter_limits(find_parameter(label), bounds, step);
    }

    void fitter::fix_parameter(int i)
    {
        _pars[i]._fixed = true;
    };

    void fitter::fix_parameter(std::string label)
    {
        return fix_parameter(find_parameter(label));
    };

    void fitter::free_parameter(int i)
    {
        _pars[i]._fixed = false;
    };

    void fitter::free_parameter(std::string label)
    {
        return free_parameter(find_parameter(label));
    };

    std::vector<double> fitter::convert(const double * pars)
    {
        std::vector<double> result;
        for (int i = 0; i < _pars.size(); i++)
        {
            result.push_back(pars[i]);
        };
        return result;
    };

    // ---------------------------------------------------------------------------
    // Chi-squared functions

    // Total chi2, this combines all data sets and is the function that gets called 
    // by the minimizer
    double fitter::fit_chi2(const double * cpars)
    {
        // First convert the C string to a C++ vector
        std::vector<double> pars = convert(cpars);

        // Pass parameters to the amplitude
        _amplitude->set_parameters(pars);


        // Then sum over all data sets and observables
        double chi2 = 0;
        for (auto data : _int_data)
        {
            chi2 += chi2_int(data);
        };
        for (auto data : _dif_data)
        {
            chi2 += chi2_dif(data);
        };

        return chi2;
    };

    // Calculate the chi2 from integrated xsection data
    double fitter::chi2_int(data_set data)
    {
        // Sum over data points
        double chi2 = 0;
        for (int i = 0; i < data._N; i++)
        {
            double W = (data._lab) ? W_cm(data._w[i]) : data._w[i];
            double s = W*W;

            double sigma_th = _amplitude->integrated_xsection(s);
            double sigma_ex = data._obs[i];
            double error    = data._obserr[i];

            chi2 += pow((sigma_th - sigma_ex)/error, 2);
        };

        return chi2;
    };

    // Calculate chi2 from differential data set
    double fitter::chi2_dif(data_set data)
    {   
        double chi2 = 0;
        for (int i = 0; i < data._N; i++)
        {
            double W = (data._lab) ? W_cm(data._w[i]) : data._w[i];
            double s = W*W;

            double t = (data._negt) ? -data._t[i] : data._t[i];
            if (data._tprime) t += _amplitude->_kinematics->t_min(s);

            double sigma_th = _amplitude->differential_xsection(s, t);
            double sigma_ex = data._obs[i];
            double error    = data._obserr[i];

            chi2 += pow((sigma_th - sigma_ex)/error, 2);
        };

        return chi2;
    };
    
    // ---------------------------------------------------------------------------
    // Load parameter information into minuit

    void fitter::set_up(std::vector<double> starting_guess)
    {
        _minuit->Clear();
        _minuit->SetTolerance(_tolerance);
        _minuit->SetPrintLevel(_print_level);
        _minuit->SetMaxFunctionCalls(_max_calls);

        for (int a = 0; a < _pars.size(); a++)
        {   
            _minuit->SetVariable(a, _pars[a]._label, starting_guess[a], _pars[a]._step);

            if (_pars[a]._custom_limits)
            {
                _minuit->SetVariableLimits(a, _pars[a]._lower, _pars[a]._upper);
            }
            if (_pars[a]._fixed)
            {
                _minuit->FixVariable(a);
            }
        };

        fcn = ROOT::Math::Functor(this, &fitter::fit_chi2, _pars.size());
        _minuit->SetFunction(fcn);
    };

    // Actually do the fit given a vector of size amp->N_pars() as starting values
    // Prints results to command line and sets best fit parameters into _amplitude
    void fitter::do_fit(std::vector<double> starting_guess, bool show_data)
    {
        if (starting_guess.size() != _pars.size()) 
        {
            warning("fitter::do_fit", "Starting guess not the correct size!");
            return;
        };

        set_up(starting_guess);

        if (show_data) { line(); data_info(); };
        line(); divider(); 
        line(); parameter_info(starting_guess, true);
        line(); divider(); line();

        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Beginning fit..." << std::flush; 

        if (_print_level != 0) line();   
        _minuit->Minimize();
        if (_print_level != 0) line();   

        std::cout << "Done! \n";
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast< std::chrono::seconds>(stop - start);
        std::cout << std::left << "Finished in " << duration.count() << " s" << std::endl;
        line();
        print_results();
    };

    // Repeat do_fit N times and find the best fit
    // Parameters are randomly initialized each time on the interval [-5, 5] unless custom limits are set
    void fitter::do_fit(int N)
    {
        // Initial guess
        std::vector<double> guess;

        for (int i = 1; i <= N; i++)
        {
            guess.clear();

            // Whether this is the first iteration
            bool first_fit =  !(i-1);
            
            // Initialize the guess for each parameter
            for (auto par : _pars)
            {
                if (par._custom_limits)
                { guess.push_back(_guesser->Uniform(par._lower, par._upper)); }
                else                    
                { guess.push_back(_guesser->Uniform(-5., 5));}
            };

            // Do out fit with this random guess
            std::cout << std::left << "Fit #" + std::to_string(i) << std::endl;
            do_fit(guess, first_fit);

            // Compare with previous best and update
            if ( first_fit || chi2dof() < _best_chi2dof)
            {
                _best_chi2dof = _minuit->MinValue() / (_N - _minuit->NFree());
                _best_chi2    = _minuit->MinValue();
                _best_pars    = convert(_minuit->X());
            };
        };

        // After looping, set the best_pars to the amplitude
        _amplitude->set_parameters(_best_pars);

        // And set the global saved pars to the best_fit
        std::cout << std::left << "Best fit found after N = " + std::to_string(N) + " iterations" << std::endl;
        print_results(false);
    };

    // ---------------------------------------------------------------------------
    // Status messages printed to command line

    // Summary of data sets that have been recieved
    void fitter::data_info()
    {
        using std::cout; 
        using std::left;
        using std::setw;
        using std::endl;

        cout << left;
        cout << "Fitting amplitude (\"" << _amplitude->id() << "\") to " << _N << " data points:" << endl;
        line();
        cout << setw(30) << "DATA SET"         << setw(20) << "OBSERVABLE"     << setw(10) << "POINTS" << endl;
        cout << setw(30) << "----------------" << setw(20) << "--------------" << setw(10) << "-------" << endl;
        for (auto data : _int_data)
        {
            cout << setw(30) << data._id  << setw(20)  << "integrated xs"   << setw(10) << data._w.size()  << endl;  
        };
        for (auto data : _dif_data)
        {   
            cout << setw(30) << data._id  << setw(20) << "differential xs" << setw(10) << data._t.size() << endl;  
        };
    };

    // Print out a little table of the current status of parameters
    void fitter::parameter_info(std::vector<double> starting_guess, bool start)
    {  
        using std::cout; 
        using std::left;
        using std::setw;
        using std::endl;

        cout << std::setprecision(10);
        cout << left;

        // Print message at the beginning of the fit
        if (start)  cout << "Fitting " + std::to_string(_minuit->NFree()) << " (out of " << std::to_string(_minuit->NDim()) << ") parameters:" << endl; 
        line();
        std::string column_3;
        (start) ? (column_3 = "START VALUE") : (column_3 = "FIT VALUE");
        cout << setw(10) << "N"     << setw(20) << "PARAMETER"  << setw(10) << column_3       << endl;
        cout << setw(10) << "-----" << setw(20) << "----------" << setw(10) << "------------" << endl;

        for (int i = 0; i < _pars.size(); i++)
        {
            // Parse whether a parameter has extra options 
            // such as custom limits
            std::string extra = "";
            if (_pars[i]._custom_limits && start)
            {   
                std::stringstream ss;
                ss << std::setprecision(5) << "[" << _pars[i]._lower << ", " << _pars[i]._upper << "]";
                extra = ss.str();
            };
            // Or is fixed
            if (_pars[i]._fixed && start) extra = "FIXED";


            cout << left << setw(10) << i << setw(20) << _pars[i]._label << setw(20) << starting_guess[i] << setw(10) << extra << endl;
        };
    };

    // At the end of a fit, print out a table sumarizing the fit results
    void fitter::print_results(bool last_fit)
    {
        int dof        = _N - _minuit->NFree();
        double chi2    = (last_fit) ?  _minuit->MinValue() : _best_chi2;
        double chi2dof = (last_fit) ? _minuit->MinValue() / double(dof) : _best_chi2dof;
        std::vector<double> pars = (last_fit) ? convert(_minuit->X()): _best_pars;

        divider();
        std::cout << std::left << std::setw(10) << "chi2 = "      << std::setw(15) << chi2;
        std::cout << std::left << std::setw(10) << "chi2/dof = "  << std::setw(15) << chi2dof << "\n";
        line();
        parameter_info(pars, false);
        divider();
        line();
        
        // At the end update the amplitude parameters to include the fit results
        _amplitude->set_parameters(pars);

        _chi2     = chi2;
        _chi2dof  = chi2dof;
        _fit_pars = pars;

        // Let rest of the fitter that a fit result has been saved
        _fit = true;
    }
};