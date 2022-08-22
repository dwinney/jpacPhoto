// This discontinuity relied on PW projections of the disconinuity saved to file 
// as a interpolation_2D grid
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef INTERP_DISC
#define INTERP_DISC 

#include "box_discontinuity.hpp"
#include "interpolation_2D.hpp"
#include <functional>

namespace jpacPhoto
{
    class interpolated_discontinuity : public box_discontinuity
    {
        public:

        // Constructor needs just a kinematics object to get masses and number of helicity combinations
        // and jmax for the number of terms in the PW expansion to expect
        interpolated_discontinuity(reaction_kinematics * xkinem)
        : box_discontinuity(xkinem), _jmax(0),  _nAmps(xkinem->num_amps())
        {};

        // Constructor needs just a kinematics object to get masses and number of helicity combinations
        // and jmax for the number of terms in the PW expansion to expect
        interpolated_discontinuity(reaction_kinematics * xkinem, int jmax, bool verbose = false)
        : box_discontinuity(xkinem), _jmax(jmax), _nAmps(xkinem->num_amps())
        {
            set_Jmax(jmax, verbose);
        };

        // Destructor needs to clean up all the pointers we created
        ~interpolated_discontinuity()
        {
            delete_pointers();
        };

        // Evaluate the discontinutiy by summing helicity amplitudes with corresponding d-functions
        std::complex<double> helicity_amplitude(std::array<int,4> helicities, double s, double t);

        // Evaluate the dispersion relation to spit out the partial-wave amplitude
        std::complex<double> helicity_pwa(int i, int j, double s);
        std::complex<double> helicity_pwa(std::array<int,4> helicities, int j, double s)
        {
            int index = _kinematics->helicity_index(helicities);
            return helicity_pwa(index, j, s);
        };

        // For a selected jmax we set up j * nAmps 2D interpolations 
        void set_Jmax(int J, bool verbose = false)
        {
            _jmax = J;

            if (_hpw_projections.size() != 0) delete_pointers();

            for (int j = 0; (2*j+1) <= J; j++)
            {
                std::vector<interpolation_2D*> jth_row;

                // Only need half of the amplitudes
                for (int i = 0; i < _kinematics->num_amps()/2 ; i++)
                {
                    jth_row.push_back( new interpolation_2D(verbose) );
                };

                _hpw_projections.push_back(jth_row);
            };
        };

        // Parameter setting and getting
        int  get_nParams(){ return 2; };
        void set_params(std::vector<double> params) { _scut = params[0]; _eta = params[1]; };

        // Set the file path and prefix for where to search for grids
        // Files are assumed to be in the format:
        // XXXJ_%_H_%.dat where XXX is the prefix and includes the directory path
        void set_import_prefix(std::string x){ _prefix = x; };

        // Clear all the saved data in the grids but keep the pointer instances
        void clear_data()
        {
            for (int j = 0; j < _hpw_projections.size() ; j++)
            {
                for (int i = 0; i < _hpw_projections[j].size(); i++)
                {
                    _hpw_projections[j][i]->clear_grid();
                };
            };

            _dataImported = false;
        }

        // Grab data from file in the preset _prefix
        void import_data()
        {
            for (int j = 0; j < _jmax; j++)
            {
                for (int i = 0; i < _hpw_projections[j].size(); i++)
                {
                    std::string filename = _prefix;
                    filename += "J_" + std::to_string(2*j+1);
                    filename += "_H_" + std::to_string(i);
                    filename += ".dat";

                    _hpw_projections[j][i]->import_grid(filename);
                };
            };

            _dataImported = true;
        };

        private:

        int    _jmax;   // Maximal s-channel spin in PWA expansion
        int    _nAmps;  // Number of helicity amplitudes to import
        double _scut;   // s-channel cut-off in dispersion relation
        double _eta;    // t-channel cut-off in form-factor

        std::string _prefix = "./";
        bool _dataImported = false;

        // Does exactly what it says
        void delete_pointers()
        {
            for (int j = 0; j < _hpw_projections.size() ; j++)
            {
                for (int i = 0; i < _hpw_projections[j].size(); i++)
                {
                    delete _hpw_projections[j][i];
                };
            };  
        };

        // Take in a function f and calculate the disperion integral of it
        std::complex<double> dispersion(std::function<double(double)> f, double s);

        // Vector of vector of interpolations
        // Rows correspond to J value
        // Columns are helicity amplitude index
        std::vector<std::vector<interpolation_2D*>> _hpw_projections;
    };
};

#endif  