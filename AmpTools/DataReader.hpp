// Abstract data reader class for an N particle final state
// This is the analogue of the jpacPhoto::DataWriter class and similarly
// should be overriden in order to read extra quantities
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               University of Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------

#ifndef DATA_READER_HPP
#define DATA_READER_HPP

#include <array>
#include <cassert>
#include <iostream>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSystem.h"
#include "IUAmpTools/Kinematics.h"
#include "IUAmpTools/DataReader.h"
#include "IUAmpTools/UserDataReader.h"
#include "constants.hpp"
#include "print.hpp"

namespace jpacPhoto
{
    template<int N>
    class DataReader : public UserDataReader<DataReader<N>>
    {
        public:

        // Basic constructor
        DataReader<N>()
        : UserDataReader<DataReader<N>>()
        {};

        // Constructor with the input file name and tree name
        DataReader<N>(const std::vector<std::string>& args)
        : UserDataReader<DataReader<N>>(args), _nCounter(0)
        {
            initialize(args);
        };

        inline Kinematics* getEvent()
        {
            if (_nCounter > numEvents() || _ERROR) return nullptr;

            _tree->GetEntry(_nCounter++);
            std::vector<TLorentzVector> particleList;
            particleList.push_back(TLorentzVector(0., 0., _Ebeam, _Ebeam));
            for (int i = 0; i < N; i++)
            {
                particleList.push_back(TLorentzVector(_Px[i], _Py[i], _Pz[i], _E[i]));
            };
            
            return new Kinematics(particleList, _weight);
        };

        // While getEvent() allows you to get the 4vecs per event, any saved extras
        // would beed to be recalculated. Instead grab them from the file here.
        // Extra variables should be indexed by an integer id
        enum extras {};
        virtual double getExtra(extras id){ return std::nan(""); };

        virtual void plot(std::string outfile){};

        inline void resetSource(){ _nCounter = 0; };

        inline unsigned int numEvents() const
        {
            return (!_tree) ? 0 : static_cast<unsigned int>(_tree->GetEntries());
        };

        inline std::string name() const
        {
            return "DataReader<" + std::to_string(N) + ">";
        };

        protected:

        inline void initialize(const std::vector<std::string>& args)
        {
            if (args.size() != N + 1)
            {
                warning("DataReader::initialize()", "Unexpected number of arguments passes!");
                _ERROR = true;
                return;
            };  

            // First N arguments are the particle labels
            for (int i = 0; i < N; i++)
            {
                _labels[i] = args[i];
            };
            // Last one is the file name
            std::string infile_name = args.back();

            TH1::AddDirectory( kFALSE );
            gSystem->Load( "libTree" );

            ifstream file_exists (infile_name.c_str());
            if (!file_exists)
            {
                warning("DataReader::initialize()", "Couldn't find input file " + infile_name);
                _file = nullptr; _tree = nullptr; _ERROR = true;
                return;
            };

            _file = new TFile( infile_name.c_str(), "READ" );
            _tree = static_cast<TTree*>(_file->Get(default_treename.c_str()));

            if (!_tree)
            {
                warning("DataReader::initialize()", "Couldn't find tree " + default_treename + " in file " + infile_name);
                _tree = nullptr;
                _ERROR = true;
                return;
            };

            _tree->SetBranchAddress( "E_beam", &_Ebeam);
            _tree->SetBranchAddress( "weight", &_weight);

            for (int i = 0; i < N; i++)
            {
                _tree->SetBranchAddress( ("E_" +_labels[i]).c_str(), &_E[i]);
                _tree->SetBranchAddress( ("Px_"+_labels[i]).c_str(), &_Px[i]);
                _tree->SetBranchAddress( ("Py_"+_labels[i]).c_str(), &_Py[i]);
                _tree->SetBranchAddress( ("Pz_"+_labels[i]).c_str(), &_Pz[i]);
            };

            setupExtras();
        };

        // Each implentation can allow additional branches to read from the file
        // This function is called at the end of initialize and does nothing by default
        virtual void setupExtras(){ print("noooo");};

        // Whether sometihng went wrong inthe file opening
        bool _ERROR = false; 

        // Particle labels to appear in branches
        std::array<std::string, N> _labels;

        // (Final state) Four vector elements
        std::array<double, N> _E, _Px, _Py, _Pz;

        // Beam energy
        double _Ebeam;

        // Incase using weighted events
        double _weight;

        // Root structures
        std::string default_treename = "nt";
        TFile* _file;
        TTree* _tree;
        
        int _nCounter = 0;
    };

    class TwoBodyReader : public DataReader<2>
    {
        public: 

        TwoBodyReader(const std::vector<std::string>& args)
        : DataReader<2>()
        {
            initialize(args);
        };

        enum extras{ kS, kT };
        inline double getExtra(extras var)
        {
            switch (var)
            {
                case kS: return _s;
                case kT: return _t;
                default: return std::nan("");
            };
            return std::nan("");
        };

        protected:

        double _s, _t;

        // Two extra branches for the two mandelstam variables
        inline void setupExtras()
        {
            _tree->SetBranchAddress( "s", &_s);
            _tree->SetBranchAddress( "t", &_t);
        };
    };

    class ThreeBodyReader : public DataReader<3>
    {
        public: 

        ThreeBodyReader(const std::vector<std::string>& args)
        : DataReader<3>()
        {
            initialize(args);
        };

        enum extras{ kS, kT, kS12, kS1, kS2, kT1, kT2, kCosGJ, kThetaGJ, kPhiGJ };
        inline double getExtra(extras var)
        {
            if (_ERROR) return std::nan("");
            switch (var)
            {
                case kS:       return _s;
                case kT:       return _t;
                case kS1:      return _s1;
                case kS2:      return _s2;
                case kS12:     return _s12;
                case kT1:      return _t1;
                case kT2:      return _t2;
                case kCosGJ:   return _cosGJ;
                case kThetaGJ: return _thetaGJ;
                case kPhiGJ:   return _phiGJ;
                default: return std::nan("");
            };
            return std::nan("");
        };

        protected:

        inline std::string subscript(std::string x, int i, int j = -1)
        {
            return  (j == -1) ? x + "_" + _labels[i] : x + "_" + _labels[i] + _labels[j];
        };

        double _s, _t, _s12, _s1, _s2, _t1, _t2, _cosGJ, _thetaGJ, _phiGJ;

        // Two extra branches for the two mandelstam variables
        inline void setupExtras()
        {
            _tree->SetBranchAddress( "s",    &_s);
            _tree->SetBranchAddress( "t",    &_t);
            _tree->SetBranchAddress( subscript("s", 0, 1).c_str(), &_s12);
            _tree->SetBranchAddress( subscript("s", 0, 2).c_str(), &_s1);
            _tree->SetBranchAddress( subscript("s", 1, 2).c_str(), &_s2);
            _tree->SetBranchAddress( subscript("t", 0).c_str(), &_t1);
            _tree->SetBranchAddress( subscript("t", 1).c_str(), &_t2);
            _tree->SetBranchAddress( "cosGJ",   &_cosGJ);
            _tree->SetBranchAddress( "thetaGJ", &_thetaGJ);
            _tree->SetBranchAddress( "phiGJ",   &_phiGJ);
        };
    };
};

#endif