// Methods for quickly parsing CLAS12 2026 data files into the data_set format
// Data corresponds to jpsi photoproduction near threshold off the proton
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               Instituto de Ciencias Nucleares (ICN),
//               Universidad Nacional Autonoma de México (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------
// References:
// [1] - 10.1103/lqv1-kf2n
// ------------------------------------------------------------------------------

#ifndef CLAS12_DATA_HPP
#define CLAS12_DATA_HPP

#include "jpsi/fitter.hpp"
#include "constants.hpp"
#include "kinematics.hpp"

namespace jpacPhoto
{
    namespace jpsi
    {
        namespace clas
        {
            // Integrated cross-section data  as a function of Egamma
            inline data_set integrated()
            {
                std::string id = "CLAS (2026)";
    
                auto raw = import_data<7>(data_dir()+"/clas/clas2026_int.tsv");
                int  N   = check<7>(raw, id);
    
                data_set out;
                out._id   = id;
                out._type = kIntegrated;
                out._N    = N;
                out._add_to_legend = true;
    
                out._x       = raw[2]; // <Egamma>
                out._xerr[0] = raw[3] /* raw[2] - raw[0] */; // Bin size lower bound
                out._xerr[1] = raw[3] /* raw[1] - raw[2] */; // bin size upper bound
    
                out._z       = raw[4]; // sigma
                out._zerr[0] = add_in_quadrature(raw[5], raw[6]);
                out._zerr[1] = add_in_quadrature(raw[5], raw[6]);
    
                return out;
            };
    
            // Parse the three slices available for the differential cross-section
            inline data_set slice(int sliceid)
            {
                std::string id, filename;
                double Eavg;
                switch (sliceid)
                {
                    case 0:
                    {
                        Eavg = 8.91;
                        id = "CLAS (E = 8.91 GeV)";
                        filename = data_dir()+"/clas/clas2026_diff_E0891.tsv";
                        break;
                    };
                    case 1:
                    {
                        Eavg = 9.67;
                        id = "CLAS (E = 9.67 GeV)";
                        filename = data_dir()+"/clas/clas2026_diff_E0967.tsv";
                        break;
                    };
                    case 2:
                    {
                        Eavg = 10.27;
                        id = "CLAS (E = 10.27 GeV)";
                        filename = data_dir()+"/clas/clas2026_diff_E1027.tsv";
                        break;
                    };
                    default: return data_set();
                };
    
                auto raw = import_data<9>(filename);
                int  N   = check<9>(raw, id);
    
                data_set out;
                out._id = id;
                out._N  = N;
                out._type = kDifferential;
                out._add_to_legend = true;
    
                out._x       = raw[2];
                out._xerr[0] = raw[2] - raw[0];
                out._xerr[1] = raw[1] - raw[2];
                out._y       = raw[4];
                out._z       = raw[6];
                out._zerr[0] = add_in_quadrature(raw[7], raw[8]);
                out._zerr[1] = add_in_quadrature(raw[7], raw[8]);
                out._extras["Eavg"] = Eavg;
    
                return out;
            };
    
            inline std::vector<data_set> differential()
            {
                std::vector<data_set> data; 
                for (int i = 0; i < 3; i++) data.push_back(clas::slice(i));
                return data;
            }
    
            // Indices correspond to:
            // [0] slice-0 differential
            // [1] slice-1 differential
            // [2] slice-2 differential
            // [3] integrated xsection
            inline std::vector<data_set> all()
            {
                std::vector<data_set> data = differential();
                data.push_back( clas::integrated() );
                return data;
            };
        }; // namespace clas
    }; // namespace jpsiPhoto
}; // namespace jpacphoto

#endif