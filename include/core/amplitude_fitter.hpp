// Class which takes in an amplitude and some data and produces a fit based on chi2 minimization with MINUIT
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef FITTER
#define FITTER

#include "reaction_kinematics.hpp"
#include "amplitude.hpp"

#include <sstream> 
#include <chrono> 

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

namespace jpacPhoto
{
    class amplitude_fitter
    {
        // --------------------------------------------------------------------
        public:

        // Constructor requires a reaction_kinematics object
        amplitude_fitter(amplitude * amp)
        : _amplitude(amp)
        {
            // populate parameters vector of appropriate size
            for (int i = 0; i < amp->get_nParams(); i++)
            {
                _pars.push_back(i);
            };

            _minuit = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
        };

        // Constructor requires a reaction_kinematics object
        // Optional explicit choice of minimization strategy passes to minuit object
        amplitude_fitter(amplitude * amp, std::string strategy, double tolerance = 1.E-6)
        : _amplitude(amp), _tolerance(tolerance)
        {
            // populate parameters vector of appropriate size
            for (int i = 0; i < amp->get_nParams(); i++)
            {
                _pars.push_back(i);
            };

            _minuit = ROOT::Math::Factory::CreateMinimizer("Minuit2", strategy);
        };


        // Add vectors corresponding to integrated xsection data sets
        // String id optional parameter for feeding fit results to a plotter object
        inline void add_integrated_data(std::vector<double> s, std::vector<double> sigma, std::vector<double> errors, std::string id = "")
        {
            int n = s.size();
            if (sigma.size() != n || errors.size() != n) 
            {
                std::cout << "fitter::add_integrated_data() : Input vectors are incorrect sizes!" << std::endl;
                std::cout << "Continuing without adding data set..." << std::endl;
                return;
            };

            if ( id == "" ) id = "int_data[" + std::to_string(_differential_data.size()) + "]";
            data_set new_data(s, sigma, errors, id);
            new_data._useEgamma = _globaluseEgam;

            _integrated_data.push_back(new_data);

             // Add number of points to the running totals
            _N += n; _N_int += n;
        };
    
        // Add vectors corresponding to differential xsection data sets
        // in addition to vectors for t, dsigma/dt, and errors we require the fixed energy s 
        // String id optional parameter for feeding fit results to a plotter object
        inline void add_differential_data(double s, std::vector<double> t, std::vector<double> dsigma, std::vector<double> errors, bool use_tp = false){ add_differential_data(s, t, dsigma, errors, "", use_tp); };
        inline void add_differential_data(double s, std::vector<double> t, std::vector<double> dsigma, std::vector<double> errors, std::string id = "", bool use_tp = false)
        {
            int n = t.size();
            if (dsigma.size() != n || errors.size() != n) 
            {
                std::cout << "fitter::add_differential_data() : Input vectors are incorrect sizes!" << std::endl;
                std::cout << "Continuing without adding data set..." << std::endl;
                return;
            };

            if ( id == "" ) id = "dxs_data[" + std::to_string(_differential_data.size()) + "]";
            data_set new_data(s, t, dsigma, errors, id);
            new_data._useTPrime = use_tp;
            new_data._useEgamma = _globaluseEgam;

            _differential_data.push_back(new_data);

            // Add number of points to the running totals
            _N += n; _N_dif += n;
        };
        // If each bin in t has an associated avg energy pass a vector
        inline void add_differential_data(std::vector<double> s, std::vector<double> t, std::vector<double> dsigma, std::vector<double> errors, std::string id = "", bool use_tp = false)
        {
            int n = t.size();
            if (dsigma.size() != n || errors.size() != n || s.size() != n) 
            {
                std::cout << "fitter::add_differential_data() : Input vectors are incorrect sizes!" << std::endl;
                std::cout << "Continuing without adding data set..." << std::endl;
                return;
            };

            if ( id == "" ) id = "dxs_data[" + std::to_string(_differential_data.size()) + "]";
            data_set new_data(s, t, dsigma, errors, id);
            new_data._useTPrime = use_tp;
            new_data._useEgamma = _globaluseEgam;

            _differential_data.push_back(new_data);

            // Add number of points to the running totals
            _N += n; _N_dif += n;
        };
        inline void set_parameter_labels(std::vector<std::string> labels)
        {
            if (labels.size() != _pars.size())
            {
                std::cout << "fitter::set_parameter_labels() : Input vector is of incorrect size! Expected " << _pars.size() << " but recieved " << labels.size() << "!" << std::endl;
                std::cout << "Continuing without adding labels..." << std::endl;
                return;
            };

            for (int i = 0; i < _pars.size(); i++)
            {
                _pars[i]._label = labels[i];
            }
        };

        /// Set limits of a parameter specificed by index
        inline void set_parameter_limits(int i, std::array<double,2> ranges, double step = 0.1)
        {
            _pars[i]._custom_limits = true;
            _pars[i]._lower_limit = ranges[0];
            _pars[i]._upper_limit = ranges[1];
            _pars[i]._step_size   = step;
        };
        // or specified by label
        inline void set_parameter_limits(std::string name, std::array<double,2> ranges, double step = 0.1)
        {
            for (int i = 0; i < _pars.size(); i++)
            {
                if (_pars[i]._label != name)
                {
                    if (i == _pars.size() - 1) std::cout << "amplitude_fitter::set_parameter_limit() : Parameter named " << name << " not found!" << std::endl;
                    continue;
                }

                _pars[i]._custom_limits = true;
                _pars[i]._lower_limit = ranges[0];
                _pars[i]._upper_limit = ranges[1];
                _pars[i]._step_size   = step;
                break;
            }
        };

        // Specify a parameter is fixed
        inline void fix_parameter(int i)
        {
            _pars[i]._fixed = true;
         };
        // or specified by label
        inline void fix_parameter(std::string name)
        {
            for (int i = 0; i < _pars.size(); i++)
            {
                if (_pars[i]._label != name)
                {
                    if (i ==_pars.size() - 1) std::cout << "amplitude_fitter::fix_parameter() : Parameter named " << name << " not found!" << std::endl;
                    continue;
                }

                _pars[i]._fixed = true;
                break;
            }
        };
        
        // Unfix a parameter
        inline void free_parameter(int i)
        {
            _pars[i]._fixed = false;
         };
        // or specified by label
        inline void free_parameter(std::string name)
        {
            for (int i = 0; i < _pars.size(); i++)
            {
                if (_pars[i]._label != name)
                {
                    if (i ==_pars.size() - 1) std::cout << "amplitude_fitter::free_parameter() : Parameter named " << name << " not found!" << std::endl;
                    continue;
                }

                _pars[i]._fixed = false;
            }
        };

        // Run the main loop to fit to all the stored data
        // Need an initial guess for all the parameters (even the fixed ones)
        double do_fit(std::vector<double> starting_guess);

        //Utility to change print level in TMinuit, default is to surpress all messages
        inline void set_print_level(int n){ _printLevel = n; };

        // Supplied energies of a given data set are in Egamma (lab beam energy) not in s
        // Here we set it globally
        inline void use_beam_energy(bool x = true)
        {
            _globaluseEgam = x;
            for (int i = 0; i < _integrated_data.size();   i++) _integrated_data[i]._useEgamma   = x;
            for (int i = 0; i < _differential_data.size(); i++) _differential_data[i]._useEgamma = x;
        }
        // Or for individual data sets
        inline void use_beam_energy(std::string id, bool x = true)
        {   
            bool found = false;
            // Search for the data set with given id 
            for (int i = 0; i < _integrated_data.size(); i++)
            {
                if (_integrated_data[i]._id != id) continue;
                found = true;
                _integrated_data[i]._useEgamma = x;
            };

            // Search also in the differential data sets
            for (int i = 0; i < _differential_data.size(); i++)
            {
                if (_differential_data[i]._id != id) continue;
                found = true;
                _differential_data[i]._useEgamma = x;
            };

            if (!found)
            {
                std::cout << "amplitude_fitter::use_beam_energy() : Data set with id " << id << " not found! Continuing..." << std::endl;
            }
        };

        // Whether the supplied t values in the differential xsection is tprime = t - tmin
        inline void use_tprime(std::string id, bool x = true)
        {   
            bool found = false;

            // This function only needs to search through differential data sets
            for (int i = 0; i < _differential_data.size(); i++)
            {
                if (_differential_data[i]._id != id) continue;
                found = true;
                _differential_data[i]._useTPrime = x;
            };
            if (!found)
            {
                std::cout << "amplitude_fitter::use_beam_energy() : Data set with id " << id << " not found! Continuing..." << std::endl;
            }
        };

        inline void set_max_calls(int n){ _maxCalls = n; };

        // --------------------------------------------------------------------
        private:

        // Simple container struct to hold all the relevent info for each user-added data set to fit against
        struct data_set
        {
            data_set(std::vector<double> x, std::vector<double> fx, std::vector<double> err, std::string id)
            : _id(id), _s(x), _sigma(fx), _error(err)
            {};

            data_set(double s, std::vector<double> x, std::vector<double> fx, std::vector<double> err, std::string id)
            : _id(id), _fixed_s(s), _t(x), _dsigma(fx), _derror(err), _if_fixed_s(true)
            {};

            data_set(std::vector<double> s, std::vector<double> x, std::vector<double> fx, std::vector<double> err, std::string id)
            : _id(id), _s(s), _t(x), _dsigma(fx), _derror(err), _if_fixed_s(false)
            {};

            bool _useEgamma = false;
            bool _useTPrime = false;

            std::string _id;
            std::vector<double> _s, _sigma, _error;

            bool _if_fixed_s;
            double _fixed_s;
            std::vector<double> _t, _dsigma, _derror;
        };
        
        // Similar struct but for each parameter involved in the fit
        struct parameter
        {
            parameter(int i)
            : _i(i), _label("par["+std::to_string(i)+"]")
            {};

            int _i;
            std::string _label;

            bool _fixed = false;
            double _fixed_value;

            bool _custom_limits = false;
            double _upper_limit;
            double _lower_limit;
            double _step_size = 0.1;
        };

        // minimization function
        inline double chi2(const double *pars){ return chi2_integrated(pars) + chi2_differential(pars); };

        // Chi squared just from the integrated cross-section data
        double chi2_integrated_i(int i, std::vector<double> pars);
        double chi2_integrated(const double *par);
        std::vector<double> _integrated_chi2s;

        // samething but for differential
        double chi2_differential_i(int i, std::vector<double> pars);
        double chi2_differential(const double *par);
        std::vector<double> _differential_chi2s;

        // Amplitude being fit
        amplitude * _amplitude;
        bool _globaluseEgam = false;

        // MINUIT error code
        int _printLevel = 0;
        int _maxCalls   = 1E6;
        double _tolerance = 1.E-6;
        ROOT::Math::Minimizer * _minuit;
        ROOT::Math::Functor fcn;

        // Parameter handling
        std::vector<parameter> _pars;

        // convert double[] to vector<double>
        std::vector<double> convert(const double * par)
        {
            std::vector<double> result;
            for (int n = 0; n < _pars.size(); n++)
            {
                result.push_back(par[n]);
            };
            return result;
        };

        // Saved integrated and differential cross-section data
        int _N = 0;
        int _N_int = 0, _N_dif = 0;
        std::vector<data_set> _integrated_data;
        std::vector<data_set> _differential_data;

        // Methods for printing out messages
        inline void divider()
        {
            std::cout << "--------------------------------------------------------------" << std::endl;
        };
        inline void new_line()
        {
            std::cout << std::endl;
        };

        // Print out a little table of the current status
        void data_info();
        void variable_info(std::vector<double> pars, bool opt = 1);
        void set_up(std::vector<double> starting_guess);
        void print_results();
    };
};

#endif