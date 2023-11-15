// Class to generate histograms with the JPAC style
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2023)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               University of Bonn (HISKP)
// Email:        winney@hiskp.uni-bonn.de
// ------------------------------------------------------------------------------

#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "colors.hpp"
#include "combinable.hpp"
#include <array>
#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TError.h>
#include <TLatex.h>
#include <string>

namespace jpacPhoto
{
    class plotter;

    class histogram_1D : public combinable
    {
        public:

        ~histogram_1D()
        {
            if (_canvas)
            {
                _canvas->Close();
            }
        };

        // Outputting function, generates plot and saves it to file
        void save(std::string filename, bool cumulative = false);

        // Fill the saved histogram
        inline void fill(double x, double w = 1){ _hist->Fill(x, w); };

        // Set the x labels
        inline void set_labels(std::string x, std::string y){ _xlabel = x; _ylabel = y; };

        // Set custom bounds for both axes
        inline void set_range( double xl, double xh)
        {   
            _xbounds = {xl, xh}; _customranges = true; 
        };

        protected:

        histogram_1D(TCanvas* canvas, std::string label = "", std::array<double,2> bounds = {1., 0.})
        : _canvas(canvas), _xlabel(label)
        {
            _hist = new TH1D("h", "h", 100, bounds[0], bounds[1]);
        };

        friend class plotter;

        // Draw the histogram onto whatver is the active canvas
        void draw();
        void draw_cumulative();
        void combine_draw(double scale);

        // Canvas that the plot actually gets drawn on
        TCanvas* _canvas;
        
        // Each histogram saves a single histogram instance
        TH1D * _hist;

        // Filename of where to produce the desired plot
        std::string _filename;

        // Axis labels
        std::string _xlabel = "", _ylabel = "";

        // Custom bounds for the different axes. 
        bool _customranges = false;
        std::array<double,2> _xbounds, _ybounds;

        // Things related to drawing the JPAC logo
        bool _add_logo = true;
        std::array<double,2> _logo_coords = {0.93, 0.885};
        double _logo_scale = 1;
        inline void add_logo()
        {
            int red  = +jpacColor::Red;
            int blue = +jpacColor::Blue;
            
            std::string JPAC = "#scale[1.3]{#font[32]{#color[" + std::to_string(blue) + "]{J}}"
                      + "^{#scale[0.8]{#font[32]{" + "#color[" + std::to_string(blue) + "]{P}"
                                                   + "#color[" + std::to_string(red) +  "]{A}"
                                                   + "#color[" + std::to_string(blue) + "]{C}}}}}";

            TLatex *logo = new TLatex(_logo_coords[0], _logo_coords[1], JPAC.c_str());

            logo->SetNDC();
            logo->SetTextSize(2/30. * _logo_scale);
            logo->SetTextAlign(32);
            logo->Draw();
        };
    };

    class histogram_2D : public histogram_1D
    {
        public:

        ~histogram_2D()
        {
            if (_canvas)
            {
                _canvas->Close();
            }
        };

        // Outputting function, generates plot and saves it to file
        void save(std::string filename);

        // Fill the saved histogram
        inline void fill(double x, double y, double w = 1){ _hist->Fill(x, y, w); };

        inline void set_range( std::array<double,2> x, std::array<double,2> y)
        { _xbounds = x; _ybounds = y; _customranges = true; };

        protected:

        histogram_2D(TCanvas* canvas, std::array<std::string,2> labels)
        : histogram_1D(canvas, ""), 
          _xlabel(labels[0]), _ylabel(labels[1])
        {
            _hist = new TH2D("h", "h", 60, 1., 0., 60, 1., 0.);
        };

        friend class plotter;

        // Draw the histogram onto whatver is the active canvas
        void draw();
        void combine_draw(double x) override;

        // Axis labels
        std::string _xlabel = "", _ylabel = "";
    
        // Each histogram saves a single histogram instance
        TH2D * _hist;
    };
};

#endif