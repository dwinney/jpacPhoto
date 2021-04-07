// Form of the invariant cross-section from a triple regge interaction.
//
// Author:       Daniel Winney (2021)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef TRIP_REGGE
#define TRIP_REGGE

#include "inclusive_kinematics.hpp"
#include "ffTripleRegge.hpp"
#include "jpacTripleRegge.hpp"
#include "regge_trajectory.hpp"
#include "sigma_tot.hpp"

#include "Math/IntegratorMultiDim.h"
#include "Math/GSLIntegrator.h"
#include "Math/IntegrationTypes.h"
#include "Math/Functor.h"

#include <functional>
#include <vector>
#include <tuple>

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

        // Destructor to clean up pointers
        ~triple_regge()
        {
            delete _kinematics;

            for (int i = 0; i < _termsFF.size(); i++)
            {
                delete _termsFF[i];
            }
            for (int i = 0; i < _termsJPAC.size(); i++)
            {
                delete _termsJPAC[i];
            }
        };

        //--------------------------------------------------------------------
        // Methods to add field and fox like terms
        inline void add_term(std::array<regge_trajectory*, 3> trajectories, std::function<double(double)> coupling)
        {
            auto new_term = new ffTripleRegge(_kinematics, trajectories, coupling);
            _termsFF.push_back(new_term);
        };

        //--------------------------------------------------------------------
        // Methods to add terms following Vincent's normalization
        inline void add_term(regge_trajectory* trajectory, const std::function<double(double)>& coupling, const sigma_tot * sigmatot)
        {
            auto new_term = new jpacTripleRegge(_kinematics, trajectory, coupling, sigmatot);
            _termsJPAC.push_back(new_term);
        };

        inline void set_formfactor(bool ifuse, double b)
        {
            _useFF  = ifuse;
            _cutoff = b;
        };

        //--------------------------------------------------------------------
        // dsigma / dt dM2
        double invariant_xsection(double s, double t, double M2);

        // (t, M2)
        double dsigma_dt(double s, double t);     // integrated over M2
        double dsigma_dM2(double s, double M2);   // integrated over t

        // (x, pT2)
        double dsigma_dpT2(double s, double pT2); // integrated over x
        double dsigma_dx(double s, double x);     // integrated over pT2

        // Fully integrated
        double integrated_xsection(double s);   

        inclusive_kinematics * _kinematics;
        std::string _identifier;


        //--------------------------------------------------------------------
        protected:

        bool   _useFF = false;
        double _cutoff = 0.;
        inline double form_factor_squared(double s, double t, double M2)
        {
            return exp(2. * t  / (_cutoff * _cutoff));
        };

        // Cross-section build out of sum of triple regge interaction terms
        std::vector<ffTripleRegge*>   _termsFF;
        std::vector<jpacTripleRegge*> _termsJPAC;
    };
};

#endif 