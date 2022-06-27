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
                        JPAC_pipp_withResonances, JPAC_pimp_withResonances };

    // Return a pointer to a new instance of desired cross-section
    total_xsection * get_total_xsection(sigma_option opt)
    {
        total_xsection * sigma_tot;

        switch(opt)
        {
            case PDG_pipp_onlyRegge:
            {
                sigma_tot = new PDG_parameterization(M_PION, M_PROTON, {+1., 1., 9.56, 1.767, 18.75});
                break;
            } 
            case PDG_pimp_onlyRegge:
            {
                sigma_tot = new PDG_parameterization(M_PION, M_PROTON, {-1., 1., 9.56, 1.767, 18.75});
                break;
            }
            case JPAC_pipp_onlyRegge:
            {
                sigma_tot = new JPAC_parameterization(+1, false);
                break;
            }
            case JPAC_pimp_onlyRegge:
            {
                sigma_tot = new JPAC_parameterization(-1, false);
                break;
            }
            case JPAC_pipp_withResonances:
            {
                sigma_tot = new JPAC_parameterization(+1, true);
                break;
            }
            case JPAC_pimp_withResonances:
            {
                sigma_tot = new JPAC_parameterization(-1, true);
                break;
            }
            default:
            {
                sigma_tot = new zero_xsection();
                break;
            };
        };

        return sigma_tot;
    };
};

#endif
  
  