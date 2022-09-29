// Abstract amplitude class specifically for the open-charm box diagrams relevant for jpsi photoproduction
// This should be used as a basis for implementations for specific D and D* intermediate states
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef OPENCHARM_BOX
#define OPENCHARM_BOX

#include "box_amplitude.hpp"
#include "box_discontinuity.hpp"
#include "interpolated_discontinuity.hpp"
#include "helicity_PWA.hpp"

namespace jpacPhoto
{
    enum charmed_meson {D, Dstar};

    class open_charm_box : public box_amplitude
    {
        // ---------------------------------------------------------------------------

        public:

        open_charm_box(reaction_kinematics * xkinem, charmed_meson meson, std::string id = "open_charm_box")
        : box_amplitude(xkinem, id), _charmedMeson(meson),
          _kGam(new reaction_kinematics(M_D, M_LAMBDAC)),
          _kPsi(new reaction_kinematics(M_JPSI, M_PROTON, M_D, M_LAMBDAC))
        {
            if (meson == D)
            {
                _kGam->set_meson_JP(0, -1);
                _kPsi->set_meson_JP(0, -1);
                _file_id = "dlam";
            }
            else 
            {
                _kGam->set_meson_mass(M_DSTAR);
                _kGam->set_meson_JP(1, -1);
                _kPsi->set_meson_mass(M_DSTAR);
                _kPsi->set_meson_JP(1, -1);
                _file_id = "dslam";
            }

            _nGam = _kGam->num_amps()/2; 
            _nPsi = _kPsi->num_amps()/2;
            _m    = _kGam->get_meson_mass();

            _open_charm_disc = new interpolated_discontinuity(xkinem);
            _open_charm_disc->set_intermediate_threshold(_m + M_LAMBDAC + _thOffset);
            set_discontinuity(this->_open_charm_disc);
        };

        ~open_charm_box()
        {
            delete _open_charm_disc;
            delete _kGam; delete _kPsi;
        };

        // Import data from an existing file
        // This looks for files of the type "path/FILEID_J_%_H_%.dat"
        inline void import_grid(std::string path)
        {
            _open_charm_disc->import_data( path + _file_id + "_disc_");
        };

        // Else generate a new one
        void generate_grid(double Wmax, std::array<double,2> eta_bounds, std::array<int,2> ns, std::string path = "");

        void set_params(std::vector<double> params)
        {
            check_nParams(params);
            _open_charm_disc->set_params({W_cut(params[0]), params[1]});
        };

        // Display all the messages with grid/file handling
        void set_verbose(bool x){ _verbose = x; };

        // Whether to use the alternative phase for the psi-opencharm lagrangians
        void use_alternative_phase(bool x){ _altPhase = x; };

        // ---------------------------------------------------------------------------

        protected:

        int _Jmax;

        // Related to the D or Ds meson
        charmed_meson _charmedMeson;
        double _m;

        // ID tag to add to filenames to differentiate the two amplitudes
        std::string _file_id;

        // Integral cutoff as a function of qmax in GeV
        inline double W_cut(double qmax)
        {
            return sqrt(qmax*qmax + /* _m*_m */ M_D*M_D) + sqrt(qmax*qmax + M_LAMBDAC*M_LAMBDAC);
        };

        bool _verbose = false;

        // The discontinuity comes from pre-generated tables 
        interpolated_discontinuity * _open_charm_disc;

        // The sub-process amplitudes
        reaction_kinematics * _kGam, * _kPsi;
        amplitude * _gamp_amp, * _psip_amp;

        // These functions are called at each set of eta in the grid generation
        // they should be overwritten with somethint to update the sub-process amplitudes at each step
        // this is amplitude dependent, by default we do nothing!
        virtual void update_gamp(double eta){};
        virtual void update_psip(double eta){};

        // Number of amplitudes to expect
        int _nGam, _nPsi, _nBox = 12; // Box asssumed to always has jpsi in final state so always 12

        // Two-body phase-space
        inline double phase_space(double s)
        {
            if   ( sqrt(s) < (M_LAMBDAC + _m) ) return 0.;
            return sqrt( Kallen(s, _m*_m, M_LAMBDAC*M_LAMBDAC ) ) / sqrt(4. * s) / (8. * PI * sqrt(s));
        };

        // Grid parameters
        std::array<int,   2> _ns;
        std::array<double,2> _wbounds, _etabounds;
        std::string          _path;

        bool _altPhase = false;
        inline double psi_phase(){ return -1. * (!_altPhase) + 1. * (_altPhase); };

        // Universal couplings
        double _lambdaQCD = 0.25;
        double _gGamDDs   = 0.134;
        double _gGamDsDs  = 0.641;
        double _gPsiDD    =  7.44;
        double _gPsiDDs   =  3.84;
        double _gPsiDsDs  =  7.99;
        double _gPsiLL    = -1.4;
        double _gGamLL    = E;
        double _gDsNL     = -4.3;
        double _gDNL      = -13.2;

        double _thOffset = 1.E-4;

        // In order to generate the grids for the discontinuity we need them for the B and C amplitudes
        void generate_gamp_grid();
        void generate_psip_grid();
    };
};
#endif 