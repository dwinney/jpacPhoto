// Arrays containing helicity combinations for indexing.
// Moved to a seperate file to not clutter up reaction_kinematics.hpp
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _HELIC_COMBO_
#define _HELIC_COMBO_

#include <iostream>
#include <vector>
#include <array>

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // Massless helicity combinations
    const std::vector< std::array<int, 4> > SPIN_ZERO_HELICITIES =
    {
    //  {  γ,  p,  S,  p'}
        {  1, -1,  0, -1}, // 0
        {  1, -1,  0,  1}, // 1
        {  1,  1,  0, -1}, // 2
        {  1,  1,  0,  1}, // 3
        { -1, -1,  0, -1}, // 4
        { -1, -1,  0,  1}, // 5
        { -1,  1,  0, -1}, // 6
        { -1,  1,  0,  1}  // 7
    };

    // POS / NEG refer to ordering whether hel[0] is +-1 everything else equal
    const std::vector<int> SPIN_ZERO_POS_ITERS = {0, 1, 2, 3, 4, 5, 6, 7};
    const std::vector<int> SPIN_ZERO_NEG_ITERS = {4, 5, 6, 7, 0, 1, 2, 3};

    const std::vector< std::array<int, 4> > SPIN_ONE_HELICITIES =
    {
    //  {  γ,  p,  V,  p'}
        {  1, -1,  1, -1}, // 0
        {  1, -1,  1,  1}, // 1
        {  1, -1,  0, -1}, // 2
        {  1, -1,  0,  1}, // 3
        {  1, -1, -1, -1}, // 4
        {  1, -1, -1,  1}, // 5
        {  1,  1,  1, -1}, // 6
        {  1,  1,  1,  1}, // 7
        {  1,  1,  0, -1}, // 8
        {  1,  1,  0,  1}, // 9
        {  1,  1, -1, -1}, // 10
        {  1,  1, -1,  1}, // 11
        { -1, -1,  1, -1}, // 12
        { -1, -1,  1,  1}, // 13
        { -1, -1,  0, -1}, // 14
        { -1, -1,  0,  1}, // 15
        { -1, -1, -1, -1}, // 16
        { -1, -1, -1,  1}, // 17
        { -1,  1,  1, -1}, // 18
        { -1,  1,  1,  1}, // 19
        { -1,  1,  0, -1}, // 20
        { -1,  1,  0,  1}, // 21
        { -1,  1, -1, -1}, // 22
        { -1,  1, -1,  1}  // 23
    };

    // Correspond to amplitudes with hel[2] = +1
    // POS / NEG refer to ordering whether hel[0] is +-1 everything else equal
    const std::vector<int> SPIN_ONE_POS_ITERS = {0, 1, 6, 7, 12, 13, 18, 19};
    const std::vector<int> SPIN_ONE_NEG_ITERS = {12, 13, 18, 19, 0, 1, 6, 7};

    const std::vector< std::array<int, 4> > SPIN_TWO_HELICITIES =
    {
    //  {  γ,  p,  V,  p'}
        {  1, -1,  2, -1}, //  0
        {  1, -1,  2,  1}, //  1
        {  1, -1,  1, -1}, //  2
        {  1, -1,  1,  1}, //  3
        {  1, -1,  0, -1}, //  4
        {  1, -1,  0,  1}, //  5
        {  1, -1, -1, -1}, //  6
        {  1, -1, -1,  1}, //  7
        {  1, -1, -2, -1}, //  8
        {  1, -1, -2,  1}, //  9
        {  1,  1,  2, -1}, // 10
        {  1,  1,  2,  1}, // 11
        {  1,  1,  1, -1}, // 12
        {  1,  1,  1,  1}, // 13
        {  1,  1,  0, -1}, // 14
        {  1,  1,  0,  1}, // 15
        {  1,  1, -1, -1}, // 16
        {  1,  1, -1,  1}, // 17
        {  1,  1, -2, -1}, // 18
        {  1,  1, -2,  1}, // 19
        { -1, -1,  2, -1}, // 20
        { -1, -1,  2,  1}, // 21
        { -1, -1,  1, -1}, // 22
        { -1, -1,  1,  1}, // 23
        { -1, -1,  0, -1}, // 24
        { -1, -1,  0,  1}, // 25
        { -1, -1, -1, -1}, // 26
        { -1, -1, -1,  1}, // 27
        { -1, -1, -2, -1}, // 28
        { -1, -1, -2,  1}, // 29
        { -1,  1,  2, -1}, // 30
        { -1,  1,  2,  1}, // 31
        { -1,  1,  1, -1}, // 32
        { -1,  1,  1,  1}, // 33
        { -1,  1,  0, -1}, // 34
        { -1,  1,  0,  1}, // 35
        { -1,  1, -1, -1}, // 36
        { -1,  1, -1,  1}, // 37
        { -1,  1, -2, -1}, // 38
        { -1,  1, -2,  1}  // 39
    };

    // Correspond to amplitudes with hel[2] = +1
    // POS / NEG refer to ordering whether hel[0] is +-1 everything else equal
    const std::vector<int> SPIN_TWO_POS_ITERS = {2, 3, 12, 13, 22, 23, 32, 33};
    const std::vector<int> SPIN_TWO_NEG_ITERS = {22, 23, 32, 33, 2, 3, 12, 13};

    // ---------------------------------------------------------------------------
    // If we allow massive "beam"
    const std::vector< std::array<int, 4> > MASSIVE_SPIN_ZERO_HELICITIES =
    {
    //  {  γ,  p,  S,  p'}
        {  1, -1,  0, -1},
        {  1, -1,  0,  1},
        {  1,  1,  0, -1},
        {  1,  1,  0,  1},
        {  0, -1,  0, -1},
        {  0, -1,  0,  1},
        {  0,  1,  0, -1},
        {  0,  1,  0,  1},
        { -1, -1,  0, -1},
        { -1, -1,  0,  1},
        { -1,  1,  0, -1},
        { -1,  1,  0,  1}
    };
  
    const std::vector< std::array<int, 4> > MASSIVE_SPIN_ONE_HELICITIES =
    {
    //  {  γ,  p,  V,  p'}
        {  1, -1,  1, -1},
        {  1, -1,  1,  1},
        {  1, -1,  0, -1},
        {  1, -1,  0,  1},
        {  1, -1, -1, -1},
        {  1, -1, -1,  1},
        {  1,  1,  1, -1},
        {  1,  1,  1,  1},
        {  1,  1,  0, -1},
        {  1,  1,  0,  1},
        {  1,  1, -1, -1},
        {  1,  1, -1,  1},
        {  0, -1,  1, -1},
        {  0, -1,  1,  1},
        {  0, -1,  0, -1},
        {  0, -1,  0,  1},
        {  0, -1, -1, -1},
        {  0, -1, -1,  1},
        {  0,  1,  1, -1},
        {  0,  1,  1,  1},
        {  0,  1,  0, -1},
        {  0,  1,  0,  1},
        {  0,  1, -1, -1},
        {  0,  1, -1,  1},
        { -1, -1,  1, -1},
        { -1, -1,  1,  1},
        { -1, -1,  0, -1},
        { -1, -1,  0,  1},
        { -1, -1, -1, -1},
        { -1, -1, -1,  1},
        { -1,  1,  1, -1},
        { -1,  1,  1,  1},
        { -1,  1,  0, -1},
        { -1,  1,  0,  1}, 
        { -1,  1, -1, -1},
        { -1,  1, -1,  1}
    };

    inline std::vector<std::array<int, 4>> get_helicities(int J, double m = 0.)
    {
        if (m > 1.E-3) 
        {
            switch (J)
            {   
                case 0: return MASSIVE_SPIN_ZERO_HELICITIES;
                case 1: return MASSIVE_SPIN_ONE_HELICITIES;
                default:
                {
                    std::cout << "Error! Amplitudes for spin J = " << J << " not yet implemented. Quitting...\n";
                    exit(0);
                }
            };
        }
        else
        {
            switch (J)
            {   
                case 0: return SPIN_ZERO_HELICITIES;
                case 1: return SPIN_ONE_HELICITIES;
                case 2: return SPIN_TWO_HELICITIES;
                default:
                {
                    std::cout << "Error! Amplitudes for spin J = " << J << " not yet implemented. Quitting...\n";
                    exit(0);
                }
            };
        }
        
        return {};
    };

    inline std::array<std::vector<int>, 2> get_iters(int J)
    {
        switch (J)
        {
            case 0: return {SPIN_ZERO_POS_ITERS, SPIN_ZERO_NEG_ITERS};
            case 1: return {SPIN_ONE_POS_ITERS,  SPIN_ONE_NEG_ITERS};
            case 2: return {SPIN_TWO_POS_ITERS,  SPIN_TWO_NEG_ITERS};
            default:
            {
                std::cout << "Error! Amplitudes for spin J = " << J << " not yet implemented. Quitting...\n";
                exit(0);
            }
        };
    };

    inline int find_helicity(std::array<int, 4> helicities, int j, double m = 0.)
    {
        std::vector<std::array<int,4>> hels = get_helicities(j, m);

        for (int i = 0; i < hels.size(); i++)
        {
            if (helicities == hels[i])
            {
                return i;
            }
        };

        std::cout << "Error cannot find helicities! \n";  
        return -1;
    };
};

#endif