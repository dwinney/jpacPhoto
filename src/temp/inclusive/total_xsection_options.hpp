// Phenomenological expressions for the total cross-sections.
// We use a generic class callable by double total_xsection(double) to select different
// parameterizations or reactions
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef SIGMA_TOT_OPT
#define SIGMA_TOT_OPT

#include "sigma_JPAC.hpp"
#include "sigma_PDG.hpp"
#include "total_xsection.hpp"

namespace jpacPhoto
{
    // All the total cross-sections implemented so far
    enum sigma_option { PDG_pipp_onlyRegge,       PDG_pimp_onlyRegge, 
                        JPAC_pipp_onlyRegge,      JPAC_pimp_onlyRegge,
                        JPAC_pipp_withResonances, JPAC_pimp_withResonances,
                        JPAC_pipp_onlyDelta };

    // Return a pointer to a new instance of desired cross-section
    total_xsection * get_total_xsection(sigma_option opt);
};

#endif
  
  