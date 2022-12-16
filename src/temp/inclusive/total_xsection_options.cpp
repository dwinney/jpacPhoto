// Phenomenological expressions for the total cross-sections.
// We use a generic class callable by double total_xsection(double) to select different
// parameterizations or reactions
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "total_xsection_options.hpp"

jpacPhoto::total_xsection * jpacPhoto::get_total_xsection(sigma_option opt)
{
        total_xsection * sigma_tot = NULL;

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
            case  JPAC_pipp_onlyDelta:
            {
                sigma_tot = new JPAC_parameterization(+1, 1);
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