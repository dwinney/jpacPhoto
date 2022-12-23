// Reworking of the jpacStyle as internal library special for jpacPhoto
// Define color set used for JPAC style plots
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#include "plotter.hpp"

namespace jpacPhoto
{
    // Set up the colors we can use
    void plotter::initialize_colors()
    {
        // make sure each color gets a free index
        jpacBlue     = new TColor(+jpacColor::Blue,     0.12156862745098039, 0.4666666666666667,  0.7058823529411765);
        jpacRed      = new TColor(+jpacColor::Red,      0.8392156862745098,  0.15294117647058825, 0.1568627450980392);
        jpacGreen    = new TColor(+jpacColor::Green,    0.0,                 0.6196078431372549,  0.45098039215686275);
        jpacOrange   = new TColor(+jpacColor::Orange,   0.8823529411764706,  0.4980392156862745,  0.054901960784313725);
        jpacPurple   = new TColor(+jpacColor::Purple,   0.5803921568627451,  0.403921568627451,   0.7411764705882353);
        jpacBrown    = new TColor(+jpacColor::Brown,    0.5490196078431373,  0.33725490196078434, 0.29411764705882354);
        jpacPink     = new TColor(+jpacColor::Pink,     0.8901960784313725,  0.4666666666666667,  0.7607843137254902);
        jpacGold     = new TColor(+jpacColor::Gold,     0.7372549019607844,  0.7411764705882353,  0.13333333333333333);
        jpacAqua     = new TColor(+jpacColor::Aqua,     0.09019607843137255, 0.7450980392156863,  0.8117647058823529);
        jpacGrey     = new TColor(+jpacColor::Grey,     0.4980392156862745,  0.4980392156862745,  0.4980392156862745);   
        jpacDarkGrey = new TColor(+jpacColor::DarkGrey, 0.31372549019,       0.31372549019,       0.31372549019);
    };

    // Set up all the default style options
    void plotter::initialize_style()
    {
        // Remove output messages when printing to file
        gErrorIgnoreLevel = kWarning;

        // remove info box
        _style->SetOptStat(0);

        // Centre title
        _style->SetTitleAlign(22);
        _style->SetTitleX(.5);
        _style->SetTitleY(.95);
        _style->SetTitleBorderSize(0);

        // set background colors to white
        _style->SetFillColor(10);
        _style->SetFrameFillColor(10);
        _style->SetCanvasColor(10);
        _style->SetPadColor(10);
        _style->SetTitleFillColor(0);
        _style->SetStatColor(10);

        // Don't put a colored frame around the plots
        _style->SetFrameBorderMode(0);
        _style->SetCanvasBorderMode(0);
        _style->SetPadBorderMode(0);

        // No border on legends
        _style->SetLegendBorderSize(0);
        _style->SetLegendTextSize(0.03);

        // Axis titles
        _style->SetNdivisions(  506, "xy");
        _style->SetTitleSize(  .045, "xyz");
        _style->SetTitleOffset(1.25, "xz");

        // More space for y-axis to avoid clashing with big numbers
        _style->SetTitleOffset(1.5, "y");

        // This applies the same settings to the overall plot title
        _style->SetTitleSize(.055, "");
        _style->SetTitleOffset(.8, "");

        // Axis labels (numbering)
        _style->SetLabelSize(  .04, "xyz");
        _style->SetLabelOffset(.01, "xyz");

        // Thicker lines
        _style->SetFrameLineWidth(2);

        // Set the tick mark style
        _style->SetPadTickX(1);
        _style->SetPadTickY(1);

        _style->SetStatFont(kjpacFont);
        _style->SetLabelFont(kjpacFont, "xyz");
        _style->SetTitleFont(kjpacFont, "xyz");
        _style->SetTitleFont(kjpacFont, ""); // Apply same setting to plot titles
        _style->SetTextFont( kjpacFont);
        _style->SetLegendFont(kjpacFont);
    };
};