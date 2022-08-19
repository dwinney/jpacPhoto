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

namespace jpacPhoto
{
    class interpolated_discontinuity : public box_discontinuity
    {
        public:

        // Constructor needs just a kinematics object to get masses and number of helicity combinations
        // and jmax for the number of terms in the PW expansion to expect
        interpolated_discontinuity(reaction_kinematics * xkinem, int jmax, bool verbose = false)
        : box_discontinuity(xkinem), _jmax(jmax), _nAmps(xkinem->num_amps())
        {
            for (int j = 0; (2*j+1) <= jmax; j++)
            {
                std::vector<interpolation_2D*> jth_row;

                // Only need half of the amplitudes
                for (int i = 0; i < xkinem->num_amps() /2 ; i++)
                {
                    jth_row.push_back( new interpolation_2D(verbose) );
                };

                _hpw_projections.push_back(jth_row);
            };
        };

        // Destructor needs to clean up all the pointers we created
        ~interpolated_discontinuity()
        {
            for (int j = 0; j < _hpw_projections.size() ; j++)
            {
                for (int i = 0; i < _hpw_projections[j].size(); i++)
                {
                    delete _hpw_projections[j][i];
                };
            };
        };

        // Evaluate the discontinutiy by summing helicity amplitudes with corresponding d-functions
        double eval(double s);

        // Parameter setting and getting
        int get_nParams(){ return 1; };
        void set_params(std::vector<double> params) { _eta = params[0]; };

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
        }

        // Grab data from file in the preset _prefix
        void import_data()
        {
            for (int j = 0; j < _hpw_projections.size() ; j++)
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
        };

        private:

        int _jmax;      // Maximal s-channel spin in PWA expansion
        int _nAmps;     // Number of helicity amplitudes to import
        double _eta;    // t-channel cut-off in form-factor

        std::string _prefix = "./";

        // Vector of vector of interpolations
        // Rows correspond to J value
        // Columns are helicity amplitude index
        std::vector<std::vector<interpolation_2D*>> _hpw_projections;
    };
};

#endif  