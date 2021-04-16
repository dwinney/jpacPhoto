// Vector production via a one-loop box diagram
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _BOX_AMP_
#define _BOX_AMP_

#include "constants.hpp"
#include "amplitudes/amplitude.hpp"
#include "amplitudes/reaction_kinematics.hpp"
#include "box/box_discontinuity.hpp"

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "Math/GSLIntegrator.h"
#include "Math/GaussLegendreIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

namespace jpacPhoto
{
    class box_amplitude : public amplitude
    {
        public: 
        // Constructor.
        // Need the parent reaction kinematics and pre-set sub-amplitudes
        box_amplitude(reaction_kinematics * xkinem, box_discontinuity * disc, std::string id = "Box Amplitude")
        : amplitude(xkinem, id), _disc(disc)
        {};

        box_amplitude(reaction_kinematics * xkinem, amplitude * left, amplitude * right, std::string id = "Box Amplitude")
        : amplitude(xkinem, id)
        {
            _disc = new box_discontinuity(left, right);
            _needDelete = true;
        };

        // Destructor
        ~box_amplitude()
        {
            if (_needDelete) delete _disc;
        };

        // Setter for max cutoff in dispersion relation
        inline void set_cutoff(double s_cut)
        {
            _s_cut = s_cut;
        };

        // only vector available
        inline std::vector<std::array<int,2>> allowedJP()
        {
            return { {1, -1} };
        };

        // Evaluate the helicity amplitude by dispersing
        std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t);

        // Override the jpacPhoto::amplitude::integrated_xsection
        double integrated_xsection(double s);

        private:        

        // Discontinutity given in terms of the two tree amplitudes
        box_discontinuity * _disc;
        bool _needDelete = false;

        // Integration momentum cutoff. Defaults to 2 GeV (an arbitrary but sensible value)
        double _s_cut = 2.;
    };
};

#endif