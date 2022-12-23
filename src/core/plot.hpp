// Reworking of the jpacStyle as internal library special for jpacPhoto
// This defines a single plot object, which is what actually prints out the files
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef PLOT_HPP
#define PLOT_HPP

#include <array>
#include <vector>
#include <iostream>   
#include <sstream> 
#include <functional>

#include <TROOT.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>

#include "data_set.hpp"
#include "amplitude.hpp"
#include "colors.hpp"

namespace jpacPhoto
{
    class plotter;

    struct entry_style
    {
        jpacColor _color    = jpacColor::DarkGrey;  // Color code 
        int  _style         = 0;                    // Either linestyle or markerstyle code
        bool _add_to_legend = true;                 // Whether to add this curve to the legend
        std::string _label  = "";                   // Label to add to Legend
    };

    // Each entry represents a curve to draw as a TGraph
    struct plot_entry 
    {
        // Constructor with just the graph
        plot_entry(TGraph* graph, bool is_data)
        : _graph(graph), _isdata(is_data)
        {};

        // Constructor with an already set up style settings
        plot_entry(TGraph* graph, entry_style xstyle, bool is_data)
        : _graph(graph), _style(xstyle), _isdata(is_data)
        {
            apply_style();
        };

        inline void set_style(entry_style xstyle)
        {
            _style = xstyle; 
            apply_style();
        }

        // Drawing options
        entry_style _style;              

        // Whether this is a data points entry or a curve
        bool _isdata = false; 

        // Graphics object
        TGraph * _graph = NULL;      

        // Apply relevant settings from _style to _graph
        inline void apply_style()
        {
            _graph->SetLineWidth(_linewidth);
            _graph->SetLineColor(+_style._color);
            _graph->SetLineStyle(_style._style);
            _graph->SetMarkerColor(+_style._color);
            _graph->SetMarkerStyle(_style._style);
        };

        // Default line-width
        static const int _linewidth = 3;
    };  

    // This class contains the entries, data, and options of producing a single plot/file
    // These can be generated from the plotter->make_plot() method which applies
    // all global settings
    class plot 
    {
        public:

        // Outputting function, generates plot and saves it to file
        void save(); 

        // Add a plot entry
        inline void add_entry(plot_entry curve)
        {
            _entries.push_back(curve);
        };
        
        // -----------------------------------------------------------------------
        // Methods to add data points to your plot

        // Convert a data_set object to a plot_entry
        void add_data(data_set data);

        // This second function can be used if you want the data_set to have
        // a different string id in the legend than the one saved in the data_set
        inline void add_data(data_set data, std::string different_id)
        {
            data_set copy(data);
            copy._id = different_id;
            add_data(copy);
        };

        // -----------------------------------------------------------------------
        // Methods to add curves from amplitudes to your plot

        // Basic function which uses the raw vectors and a string id
        void add_curve(std::vector<double> x, std::vector<double> fx, entry_style style);
        void add_curve(std::vector<double> x, std::vector<double> fx, std::string id = "");
        
        // Take in a lambda an evaluation range to get the vectors
        void add_curve(int N, std::array<double,2> bounds, std::function<double(double)> F, entry_style style);
        void add_curve(int N, std::array<double,2> bounds, std::function<double(double)> F, std::string id = "");

        // Curves added by these functions appear as dashed, not on the legend, and synced with the
        // colors of the "full" curves
        void add_dashed(std::vector<double> x, std::vector<double> fx);
        void add_dashed(int N, std::array<double,2> bounds, std::function<double(double)> F);

        // -----------------------------------------------------------------------
        // OPTION SETTERS

        // Add string labels to the axes, follows TLatex 
        inline void set_labels(std::string x, std::string y){ _xlabel = x; _ylabel = y;};
        
        // Set if the x and/or y axes is in logscale
        inline void set_logscale(bool x, bool y){ _xlog = x; _ylog = y; };

        // Set custom bounds for both axes
        inline void set_ranges( std::array<double,2> x,  std::array<double,2> y)
        { _xbounds = x; _ybounds = y; _customranges = true; };

        inline void set_legend(double x, double y){ _legendxcord = x; _legendycord = y; };
        inline void add_header(std::string x){ _header = x; _addheader = true; };
        inline void add_header(std::string variable, double value, std::string units = "")
        {
            std::ostringstream ss;
            ss << std::setprecision(3) << variable + " " << value << " " + units;
            _header = ss.str(); 
            _addheader = true; 
        };

        // -----------------------------------------------------------------------
        
        private: 
        
        // Constructor is private, only creatable through plotter
        plot(TCanvas* canvas, std::string file)
        : _canvas(canvas), _filename(file)
        {};

        friend class plotter;

        // Canvas that the plot actually gets drawn on
        TCanvas* _canvas;

        // Filename of where to produce the desired plot
        std::string _filename;

        inline void add_logo()
        {
            int red  = +jpacColor::Red;
            int blue = +jpacColor::Blue;
            std::string JPAC = "#scale[1.3]{#font[32]{#color[" + std::to_string(blue) + "]{J}}"
                      + "^{#scale[0.8]{#font[32]{" + "#color[" + std::to_string(blue) + "]{P}"
                                                   + "#color[" + std::to_string(red) +  "]{A}"
                                                   + "#color[" + std::to_string(blue) + "]{C}}}}}";
            TLatex *logo = new TLatex(.93, .885,  JPAC.c_str());

            logo->SetNDC();
            logo->SetTextSize(2/30.);
            logo->SetTextAlign(32);
            logo->Draw();
        };

        // -----------------------------------------------------------------------
        // AXIS SETUP 

        // Whether to use logscale of either axis
        bool _xlog = false, _ylog = false;
            
        // Axis labels
        std::string _xlabel = "", _ylabel = "";
        
        // Custom bounds for the different axes. 
        bool _customranges = false;
        std::array<double,2> _xbounds, _ybounds;

        // -----------------------------------------------------------------------
        // LEGEND SET UP

        bool   _addlegend     = true;
        double _legendxcord   = 0.3, _legendycord   = 0.7;
        double _legendxoffset = 0.3, _legendyoffset = 0.15;

        // Number of entries to expect on the legend
        // used to calculate the offset to be visually appealing
        int _Nlegend = 0; 
        
        bool _addheader = false;
        std::string _header = "";

        // ------------------------------------------------------------------------
        // ENTRY MANAGEMENT

        // Number of entries which are data and which are curves
        // Count used for choosing colors, etc
        int _Ndata = 0, _Ncurve = -1; 

        // List of all entries (theoretical curves) to be plotted
        std::vector<plot_entry> _entries;
    };

};

#endif