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
        {};

        box_amplitude(reaction_kinematics * xkinem, amplitude * left, amplitude * right, std::string id = "Box Amplitude")
        : amplitude(xkinem, "box_amplitude", id)
        {
            _disc = new box_discontinuity(xkinem, left, right);
            _needDelete = true;
        };

        // Destructor
        ~box_amplitude()
        {
            if (_needDelete) delete _disc;
        };

        // Setter for max cutoff in dispersion relation
        inline void set_cutoff(double s_cut){ _s_cut = s_cut; };

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
        double _s_cut = _kinematics->sth() + 2.;
    };
};

#endif