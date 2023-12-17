// Methods for indexing helicity amplitudes 
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef HELICITIES_HPP
#define HELICITIES_HPP

#include <algorithm>
#include <iostream>
#include <vector>
#include <array>
#include "debug.hpp"

namespace jpacPhoto
{
    // Each amplitude needs to be able to tell which frame its helicities are defined in
    enum helicity_frame{ HELICITY_ERROR, HELICITY_INDEPENDENT, S_CHANNEL, T_CHANNEL, U_CHANNEL };

    // Output a string of a given helicity set in format e.g. {+,+,+,+}
    inline std::string print_helicities(std::array<int,4> lam)
    {
        std::array<std::string,4> lams;
        for (int i = 0; i < 4; i++)
        {
            std::string sign = (sgn(lam[i]) == +1) ? "+" : "-";
            lams[i] = sign + std::to_string(std::abs(lam[i]));
        };
        return "[ " + lams[0] + ", " + lams[1] + ", " + lams[2] + ", " + lams[3] + "]";
    };

    // Generate a vector containing all the helicity combinations
    inline std::vector<std::array<int, 4>> get_helicities(int mJ, int bJ, bool is_massless = true)
    {
        std::vector<std::array<int,4>> output;

        // Inital state is always the same
        std::vector<int> hg = (is_massless) ? std::vector<int>({1, -1}) : std::vector<int>({1, 0, -1});
        std::vector<int> ht = {1, -1};

        // Final state depends on the spins 
        std::vector<int> hx, hr; 
        for (int i = 0;  mJ +   i >= - mJ; i--) hx.push_back(mJ +   i);
        for (int i = 0;  bJ + 2*i >= - bJ; i--) hr.push_back(bJ + 2*i);

        for (auto g : hg){
            for (auto t : ht){
                for(auto x : hx){
                    for(auto r : hr){
                        output.push_back({g, t, x, r});
                    };
                };
            };
        };
        return output;
    };

    // Given a set of helicities, find its helicity index
    inline int find_helicity(std::array<int, 4> helicities, int mj, int bj, bool is_massless = true)
    {
        std::vector<std::array<int,4>> hels = get_helicities(mj, bj, is_massless);

        auto iterator = std::find(hels.begin(), hels.end(), helicities);
        
        if (iterator != hels.end())
        {
            return iterator - hels.begin();
        }
        
        return error("find_helicity", "Cannot find helicities: " + print_helicities(helicities) + "!", -1);
    };
};

#endif