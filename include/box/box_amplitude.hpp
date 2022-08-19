// Vector production via a one-loop box diagram
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef BOX_AMP
#define BOX_AMP

#include "constants.hpp"
#include "amplitude.hpp"
#include "reaction_kinematics.hpp"
#include "box_discontinuity.hpp"

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "Math/GSLIntegrator.h"
#include "Math/GaussLegendreIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

namespace jpacPhoto
{
    class box_amplitude : public amplitude
    {
        // ---------------------------------------------------------------------------
        public: 
        // Constructor.
        // Need the parent reaction kinematics and pre-set sub-amplitudes
        box_amplitude(reaction_kinematics * xkinem, box_discontinuity * disc, std::string id = "Box Amplitude")
        : amplitude(xkinem, "box_amplitude", id), _disc(disc)
        {
            set_nParams(1 + _disc->get_nParams());
            check_JP(xkinem);
        };

        box_amplitude(reaction_kinematics * xkinem, amplitude * left, amplitude * right, std::string id = "Box Amplitude")
        : amplitude(xkinem, "box_amplitude", id)
        {
            _disc = new brute_force_discontinuity(xkinem, left, right);
            _needDelete = true;

            set_nParams(1);
            check_JP(xkinem);
        };

        // Destructor
        ~box_amplitude()
        {
            if (_needDelete) delete _disc;
        };

        // Setting utility for free parameters
        void set_params(std::vector<double> params)
        {
            check_nParams(params);
            _s_cut = params[0];
            if (_nParams > 1)
            {
                for (int i = 1; i < params.size(); i++) _disc_params.push_back(params[i]);
                if (_disc_params.size() != _disc->get_nParams())
                {
                    std::cout << "Error! # of params passed to box_amplitude, doesnt match the number needed by discontinuity. \n";
                    std::cout << "Results may vary..." << std::endl;
                } 
                _disc->set_params(_disc_params);
            }
        };

        // only vector and proton quantum numbers available currently
        inline std::vector<std::array<int,2>> allowed_meson_JP()
        {
            return { {1, -1} };
        };
        inline std::vector<std::array<int,2>> allowed_baryon_JP()
        {
            return { {1,  1} };
        };

        // Box always works with s-channel helicity projections
        inline helicity_channel helicity_CM_frame(){ return S; };

        // Evaluate the helicity amplitude by dispersing
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t);

        // Override the jpacPhoto::amplitude::integrated_xsection
        double integrated_xsection(double s);
        
        // ---------------------------------------------------------------------------
        private:        

        // Discontinutity given in terms of the two tree amplitudes
        box_discontinuity * _disc;
        bool _needDelete = false;

        // Integration momentum cutoff. Defaults to 2 GeV above threshold (an arbitrary value)
        double _s_cut;
        std::vector<double> _disc_params;
    };
};

#endif