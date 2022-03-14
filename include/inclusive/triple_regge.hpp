
// Basis class for the invariant cross-section from a triple regge interaction.
// Contains inclusve_kinematics objects as well as dynamical objects with either
// 'JPAC' or 'Field & Fox' models.
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef TRIP_REGGE
#define TRIP_REGGE

#include "inclusive_kinematics.hpp"

namespace jpacPhoto
{
    class triple_regge
    {
        public:
        // Constructor only needs a kinematics object
        triple_regge(double mass, std::string id = "")
        : _identifier(id)
        {
            _kinematics = new inclusive_kinematics(mass);
        };

    };
};

#endif 