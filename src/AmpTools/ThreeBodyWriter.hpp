// Simple data writer for a 2->3 production process with no subsequent decays
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU),
//               Helmholtz-Institut für Strahlen-und Kernphysik (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef THREEBODYWRITER_HPP
#define THREEBODYWRITER_HPP

#include <string>
#include "print.hpp"
#include "DataWriter.hpp"
#include "TLorentzVector.h"

namespace jpacPhoto
{
    class ThreeBodyWriter : public DataWriter<3>
    {
        public:

        ThreeBodyWriter( std::array<std::string,3> labels, const std::string & outfile )
        : DataWriter(labels)
        {
            initialize(outfile);
        };

        protected:

        // Calculate five invariants
        inline void setupExtras()
        {
            _tree->Branch( "s",    &_s, "s/D" );
            _tree->Branch( "t",    &_t, "t/D" );
            _tree->Branch( subscript("s", 3).c_str(), &_s12, subscriptD("s", 3).c_str());
            _tree->Branch( subscript("s", 0).c_str(), &_s1,  subscriptD("s", 0).c_str());
            _tree->Branch( subscript("s", 1).c_str(), &_s2,  subscriptD("s", 1).c_str());
            _tree->Branch( subscript("t", 0).c_str(), &_t1,  subscriptD("t", 0).c_str());
            _tree->Branch( subscript("t", 1).c_str(), &_t1,  subscriptD("t", 0).c_str());

            _tree->Branch( "cosGJ",   &_cosGJ,   "cosGJ/D");
            _tree->Branch( "thetaGJ", &_thetaGJ, "thetaGJ/D");
            _tree->Branch( "phiGJ",   &_phiGJ,   "phiGJ/D");
        };

        inline void calculateExtras(std::vector<TLorentzVector> fvecs)
        {
            // Inital state reconstructed from the internally calcualted _Ebeam
            TLorentzVector beam(  0., 0., _Ebeam, _Ebeam);
            TLorentzVector target(0., 0., 0.,     jpacPhoto::M_PROTON);
            TLorentzVector m1(fvecs[0]);
            TLorentzVector m2(fvecs[1]);
            TLorentzVector recoil(fvecs[2]);

            _s   = (beam + target).M2();
            _t   = (recoil - target).M2();
            _s12 = (m1 + m2).M2();

            // Subsystem invariant masses
            _s1  = (m1 + recoil).M2();
            _s2  = (m2 + recoil).M2();

            // Subsystem momentum transfers
            _t1 = (beam - m1).M2();
            _t2 = (beam - m2).M2();

            // Boost only the initalstate and particle 1 to GJ frame to calculate angles
            beam.Boost(  -(m1 + m2).BoostVector());
            target.Boost(-(m1 + m2).BoostVector());
            m1.Boost(    -(m1 + m2).BoostVector());

            TVector3 z = beam.Vect().Unit();
            TVector3 y = (beam.Vect().Cross(recoil.Vect())).Unit();
            TVector3 x = y.Cross(z);

            TVector3 GJ( m1.Vect().Dot(x), m1.Vect().Dot(y), m1.Vect().Dot(z));
            _thetaGJ = GJ.Theta(); _phiGJ = GJ.Phi(); _cosGJ = cos(_thetaGJ);
        };

        inline std::string subscript(std::string x, int i)
        {
            return (i < 3) ? x + "_" + _labels[i] : x + "_" + _labels[0] + _labels[1];
        };
        inline std::string subscriptD(std::string x, int i)
        {
            return subscript(x,i) + "/D";
        };

        double _s, _t, _s12, _s1, _s2, _t1, _t2, _cosGJ, _thetaGJ, _phiGJ;
    };
};

#endif
