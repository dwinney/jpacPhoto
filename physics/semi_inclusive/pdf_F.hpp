// Implementation of unpolarized structure functions based on CTEQ's leading
// order PDFs in [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        daniel.winney@iu.alumni.edu
//               dwinney@scnu.edu.cn
// ------------------------------------------------------------------------------
// REFERENCES:
// [1] - https://arxiv.org/abs/220 5.00137
// ------------------------------------------------------------------------------

#ifndef PDFS_HPP
#define PDFS_HPP

#include "constants.hpp"
#include "data_set.hpp"
#include "inclusive_function.hpp"
#include <Math/Interpolator.h>

namespace jpacPhoto
{
    class pdf_interpolator
    {
        public: 

        pdf_interpolator(){};
        ~pdf_interpolator(){ clear_data(); };

        inline double evaluate(double x, double y)
        {
            std::vector<double> fys;
            for (int i = 0; i < _ys.size(); i++) fys.push_back( _slices[i]->Eval(x) );
            _evaler.SetData(_ys, fys);
            return _evaler.Eval(y);
        };

        inline void set_data(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &f)
        {
            if (x.size() * y.size() != f.size()) warning("interpolator_2D", "Data dimensions don't match!");

            clear_data(); _ys = y;
            int dx = x.size(), dy = y.size();

            for (int i = 0; i < dy; i++)
            {
                std::vector<double> fi;
                for (int j = 0; j < dx; j++) fi.push_back(f[i + dy*j]);
                _slices.push_back( new ROOT::Math::Interpolator(x, fi, ROOT::Math::Interpolation::kCSPLINE) );
            };
        };

        private: 

        std::vector<double> _ys; 
        std::vector<ROOT::Math::Interpolator*> _slices;
        ROOT::Math::Interpolator _evaler;

        inline void clear_data()
        { 
            for (auto x : _slices) delete x;
            _slices.clear(); _ys.clear();            
        };
    };

    enum class flavor: int {bbar = 0, cbar = 1, sbar = 2, ubar = 3, dbar = 4,
                            d    = 5, u    = 6, s    = 7, c    = 8, b    = 9, g = 10};

    // Import the LO PDF grid and interpolate
    class xPDF
    {
        public:

        // Constructor doesnt need any parameters
        xPDF(){ initialize(); };
        
        inline double operator()(flavor iflav, double x, double q){ return _pdfs[int(iflav)].evaluate(x, q); };
        inline double operator()(int    iflav, double x, double q)
        { 
            if (iflav < 0 || iflav > 11) return error("Invalid flavor index!", NaN<double>());
            return _pdfs[iflav].evaluate(x, q); 
        };

        private:

        inline void initialize()
        {
            std::string file_path = "/data/CTEQ/";

            std::vector<double> xs = import_transposed<1>(file_path + "xs.dat")[0];
            std::vector<double> qs = import_transposed<1>(file_path + "qs.dat")[0];
            auto fs = import_data<11>(file_path + "CT18LO.dat");
            
            for (int i = 0; i < 11; i++) _pdfs[i].set_data(xs, qs, fs[i]);
        };

        std::array<pdf_interpolator, 11> _pdfs; // 11 flavors 
    };

    class pdf_F : public raw_inclusive_function
    {
        public: 

        pdf_F(int mode)
        : _mode(mode)
        { 
            if (mode != 1 && mode != 2) warning("pdf_F", "Invalid mode recieved!");
        };

        inline double evaluate(double M2, double t)
        {
            double x = xB(M2, t), q = sqrt(-t);
            // Keep within bounds of grid
            if (x < 1E-9  || x > 1)   return 0.;
            if (q < 1.295 || q > 1E5) return 0.;

            double sum = 0;
            for (int i = 0; i < 11; i++) sum += pow(_charges[i], 2)* _xpdf(i, x, q);
            return (_mode == 2) ? sum : sum/2/x;
        };

            inline double evaluate_xq(double x, double q)
        {
            // Keep within bounds of grid
            if (x < 1E-9  || x > 1)   return 0.;
            if (q < 1.295 || q > 1E5) return 0.;

            double sum = 0;
            for (int i = 0; i < 11; i++) sum += pow(_charges[i], 2)* _xpdf(i, x, q);
            return (_mode == 2) ? sum : sum/2/x;
        };

        private:

        int _mode = 0;

        std::array<double,11> _charges = {+1/3., -2/3., +1/3., -2/3., +1/3., -1/3., +2/3., -1/3., +2/3., -1/3., 0.};

        inline double xB(double M2, double t){ return -t/(M2 - M_PROTON*M_PROTON - t);}

        xPDF _xpdf;
    };

}; // jpacPhoto

#endif