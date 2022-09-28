// Abstract amplitude class specifically for the open-charm box diagrams relevant for jpsi photoproduction
// This should be used as a basis for implementations for specific D and D* intermediate states
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef DLAMBDAC_BOX
#define DLAMBDAC_BOX

#include "open_charm_box.hpp"
#include "pseudoscalar_exchange.hpp"
#include "vector_exchange.hpp"
#include "dirac_exchange.hpp"
#include "amplitude_sum.hpp"

namespace jpacPhoto
{
    class DLambdac_box : public open_charm_box
    {
        // ---------------------------------------------------------------------------

        public:

        DLambdac_box(reaction_kinematics * xkinem, std::string id = "DLambdac_box")
        : open_charm_box(xkinem, D, id)
        {
            // The pseudo-scalar D intermediate state only requires total J = 1/2 
            // This is because we're dominated by S-wave which is only this spin projection
            _Jmax = 1;

            _gam_dsEx   = new vector_exchange(_kGam, M_DSTAR, "D* exchange");
            _gam_lamcEx = new dirac_exchange(_kGam, M_LAMBDAC, "#Lambda_{c} exchange");

            _gam_dsEx->force_covariant(true);
            _gam_lamcEx->force_covariant(true);

            _gam_sum    = new amplitude_sum(_kGam, {_gam_dsEx, _gam_lamcEx}, "Sum");
            _gamp_amp   = _gam_sum;

            _psi_dEx    = new pseudoscalar_exchange(_kPsi, M_D, "D exchange");
            _psi_dsEx   = new vector_exchange(_kPsi, M_DSTAR, "D* exchange");
            _psi_lamcEx = new dirac_exchange(_kPsi, M_LAMBDAC, "#Lambda_{c} exchange");

            _psi_dEx->force_covariant(true);
            _psi_dsEx->force_covariant(true);
            _psi_lamcEx->force_covariant(true);

            _psi_sum    = new amplitude_sum(_kPsi, {_psi_dEx, _psi_dsEx, _psi_lamcEx}, "Sum");
            _psip_amp   = _psi_sum;
        };

        ~DLambdac_box()
        {
            delete _gam_dsEx;
            delete _gam_lamcEx;
            delete _gam_sum;
            delete _psi_dEx;
            delete _psi_dsEx;
            delete _psi_lamcEx;
            delete _psi_sum;
        };

        // ---------------------------------------------------------------------------

        protected:

        // Photon amplitudes
        vector_exchange       * _gam_dsEx;
        dirac_exchange        * _gam_lamcEx;
        amplitude_sum         * _gam_sum;

        void update_gamp(double eta)
        {
            _gam_dsEx->set_params(   {_gGamDDs, _gDsNL, 0.} );
            _gam_lamcEx->set_params( { _gGamLL,  _gDNL, 0.} );

            _gam_dsEx->set_formfactor(  2, M_DSTAR   + eta * _lambdaQCD);
            _gam_lamcEx->set_formfactor(2, M_LAMBDAC + eta * _lambdaQCD);
        };

        // Psi amplitudes
        pseudoscalar_exchange * _psi_dEx;
        vector_exchange       * _psi_dsEx;
        dirac_exchange        * _psi_lamcEx;
        amplitude_sum         * _psi_sum;

        void update_psip(double eta)
        {
            _psi_dEx->set_params(    { _gPsiDD, _gDNL     } );
            _psi_dsEx->set_params(   {_gPsiDDs, _gDsNL, 0.} );
            _psi_lamcEx->set_params( { _gPsiLL,  _gDNL, 0.} );
 
            _psi_dEx->set_formfactor(   2, M_D       + eta * _lambdaQCD);
            _psi_dsEx->set_formfactor(  2, M_DSTAR   + eta * _lambdaQCD);
            _psi_lamcEx->set_formfactor(2, M_LAMBDAC + eta * _lambdaQCD);
        };
    };
};

#endif