
#include "kmatrix/one_half.hpp"
#include "plotter.hpp"
#include "kinematics.hpp"

void test()
{
    using namespace jpacPhoto;
    
    kinematics kin = new_kinematics(M_JPSI, M_PROTON);
    kin->set_meson_JP( {1, -1} );

    auto hs = kin->helicities();
    for (int i = 0; i < hs.size(); i++) print(i, print_helicities(hs[i]));
    exit(1);
    divider();
    partial_wave one_half = new_partial_wave<kmatrix::one_half>(kin);
    one_half->partial_wave({1,1,1,1}, 25.);
};
