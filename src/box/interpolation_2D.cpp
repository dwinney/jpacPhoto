// Utility class to create interpolations of 2D data arrays
//
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "interpolation_2D.hpp"

// ---------------------------------------------------------------------------
// Take in a function of signature double(double, double) and generate grid of points
// based on saved grid parameters
void jpacPhoto::interpolation_2D::generate_grid(std::function<double(double, double)> f, bool skip_interp)
{
    // Make sure we're starting with clean slate
    clear_grid();
    if (!_limits_set)
    {
        std::cout << "interpolation_2D: Error! Grid limits not set, generate_grid() returning without change!" << std::endl;
        return;
    };
    
    // Sum over the x variable
    for (int i=0; i < _xN; i++)
    {
        // We take generate curves of y at fixed x
        double x_i = _xmin + (_xmax - _xmin) * double(i) / double(_xN - 1);
        
        // Store vectors of f(x_i, y) 
        std::vector<double> y_i, f_i;
        for (int j=0; j < _yN; j++)
        {
            double y_j = _ymin + (_ymax - _ymin) * double(j) / double(_yN - 1);
            double f_ij = f(x_i, y_j);

            y_i.push_back(y_j);
            f_i.push_back(f_ij);
        }

        _x_values.push_back(x_i);
        _y_values.push_back(y_i);
        _f_values.push_back(f_i);
    }

    if (!skip_interp) set_up_slices();
};

void jpacPhoto::interpolation_2D::set_up_slices()
{
    if ( _y_slices.size() != 0 )
    {
        for (int i = 0; i < _y_slices.size(); i++) delete _y_slices[i];
        _y_slices.clear();
    }

    for (int i=0; i < _xN; i++)
    {
        _y_slices.push_back(new ROOT::Math::Interpolator(_y_values[i], _f_values[i], ROOT::Math::Interpolation::kLINEAR));
    };
    _slices_made = true;
}

// ---------------------------------------------------------------------------
// Evaluate the interpolation at given x and y
double jpacPhoto::interpolation_2D::eval(double x, double y)
{
    // Error checks
    if (x > _xmax || x < _xmin || y > _ymax || y < _ymin)
    {
        std::cout << "interpolation_2D: Evaluting outside grid (" << x << ", " << y << "). Returning 0!" << std::endl;
        return 0.;
    }
    if (_x_values.size() == 0 || _y_slices.size() == 0)
    {
        std::cout << "interpolation_2D: No grid saved. Returning 0!" << std::endl;
        return 0.;
    }

    // Set up vector of values from the stored slices at fixed y 
    std::vector<double> fy_values;
    for (int i = 0; i < _xN; i++)
    {
        double fy = _y_slices[i]->Eval(y);
        fy_values.push_back(fy);
    }
    
    // Generate new interpolation across x
    ROOT::Math::Interpolator fxy(_x_values, fy_values, ROOT::Math::Interpolation::kLINEAR);
    
    // Return value
    return fxy.Eval(x);
};

// ---------------------------------------------------------------------------
// Can be useful to export grid to a file 
void jpacPhoto::interpolation_2D::export_grid(std::string filename)
{
    if (_x_values.size() == 0 || _y_values.size() == 0 || _f_values.size() == 0)
    {
        std::cout << "interpolation_2D: No grid is stored! export_grid() returning  without change!" << std::endl;
        return;
    }

    std::ofstream output;
    output.open(filename.c_str());

    // First line contains grid parameters
    output << std::left;
    output << std::setw(10) << _xmin << std::setw(10) << _xmax << std::setw(10) << _xN;
    output << std::setw(10) << _ymin << std::setw(10) << _ymax << std::setw(10) << _yN << std::endl;

    // Sum over the x variable
    for (int i=0; i < _xN; i++)
    {
        double x_i  = _x_values[i];

        std::vector<double> yi, fi;
        for (int j=0; j < _yN; j++)
        {
            double y_ij = _y_values[i][j];
            double f_ij = _f_values[i][j];

            output << std::left;
            output << std::setw(15) << x_i;
            output << std::setw(15) << y_ij;
            output << std::setw(15) << f_ij << std::endl;
        }
    }
    output.close();

    // Command-line messages
    if (_verbose) std::cout << "interpolation_2D: grid exported to " << filename << std::endl;
};

// ---------------------------------------------------------------------------
// Can be useful to export grid to a file 
void jpacPhoto::interpolation_2D::import_grid(std::string filename)
{
    // Open file
    std::ifstream data(filename.c_str());

    if (data.fail())
    {
        std::cout << "interpolation_2D:: Could not open file " << filename << "\n"; 
        return;
    }

    if (_verbose) std::cout << "interpolation_2D: importing grid from " << filename << std::endl;
    
    // On new import delete all saved data
    clear_grid();

    // Read first line which will give us grid parameters
    data >> _xmin >> _xmax >> _xN >> _ymin >> _ymax >> _yN; 

    // Then read the rest of the file
    double x_i, y_j, f_ij, trash;
    for (int i=0; i<_xN; i++)
    {
        // Import the zero-th row to save x_i
        data >> x_i >> y_j >> f_ij;
        _x_values.push_back(x_i);
        
        std::vector<double> ys, fs;
        ys.push_back(y_j);
        fs.push_back(f_ij);

        // Then the rest throwing away the x value
        for (int j=1; j<_yN; j++)
        {
            data >> trash >> y_j >> f_ij;
            ys.push_back(y_j);
            fs.push_back(f_ij);
        }

        if ( ys.size( )!= _yN || fs.size() != _yN ) 
        {
            std::cout << "interpolation_2D: Error! Imported data doesnt match expected size (Expected " << _yN << " but got " << ys.size() << ")...\n";
            std::cout << "Returning without change." << std::endl;
            return;
        }

        _y_values.push_back(ys);
        _f_values.push_back(fs);
    }

    if ( _y_values.size()!= _xN || _f_values.size() != _xN ) 
    {
        std::cout << "interpolation_2D: Error! Imported data doesnt match expected size (Expected " << _xN << " but got " << _y_values.size() << ")...\n";
        std::cout << "Returning without change." << std::endl;
        return;
    }

    // In addition to saving the actual grid points, at each fixed x_i we interpolate across y values
    for (int i=0; i < _xN; i++)
    {
        _y_slices.push_back(new ROOT::Math::Interpolator(_y_values[i], _f_values[i], ROOT::Math::Interpolation::kLINEAR));
    };
};