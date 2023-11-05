// Simple data writer for a 2->2 production process with no subsequent decays
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------

#ifndef TWOBODYWRITER_HPP
#define TWOBODYWRITER_HPP

#include <string>
#include "print.hpp"
#include "DataWriter.hpp"
#include "TLorentzVector.h"

namespace jpacPhoto
{
    class TwoBodyWriter : public DataWriter<2>
    {
        public:

        TwoBodyWriter( std::array<std::string,2> labels, const std::string & outfile )
        : DataWriter(labels)
        {
            initialize(outfile);
        };

        protected: 

        // Only extra variables are invariants
        inline void setupExtras()
        {
            _tree->Branch( "s", &_s, "s/D");
            _tree->Branch( "t", &_t, "t/D");
        };

        inline void calculateExtras(std::vector<TLorentzVector> fvecs)
        {
            // Invariant total energy
            _s = (fvecs[0] + fvecs[1]).M2();

            // Invariant momentum transfer
            TLorentzVector q(0., 0., _Ebeam, _Ebeam);
            _t = (q - fvecs[0]).M2();
        };

        // Mandelstam Variables
        double _s, _t;
    };
};

#endif 