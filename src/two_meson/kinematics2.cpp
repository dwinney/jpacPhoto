// Class to contain all relevant kinematic quantities. 
// (e.g. angles, momenta, helicity info)
//
// Unlike single_meson::kinematics this class represents a 2->3 process
// with two pseudo-scalars and a proton in the final state
// Generalizations to more complicated final states can be implemented lated
//
// Since we dont need to specify spins, everything is determined by each external mass:
// _mB = beam particle mass (default 0)
// _mT = target recoil mass (default M_PROTON)
// _m1 = "meson 1"
// _m2 = "meson 2"
// _mR = recoil baryon mass (default M_PROTON)
//
// --------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC)
//               University of Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// --------------------------------------------------------------------------------

#include "kinematics2.hpp"

namespace jpacPhoto
{
    namespace two_meson
    {
        // ---------------------------------------------------------------------------
        // Given a set of helicities, return the index its saved under
        int raw_kinematics::helicity_index(std::array<int,3> helicities)
        {
            auto iterator = std::find(_helicities.begin(), _helicities.end(), helicities);
            if (iterator != _helicities.end()) return iterator - _helicities.begin();
            return error("find_helicity", "Cannot find helicities: " + jpacPhoto::print_helicities(helicities) + "!", -1);
        };        
    };
};