// This provides implementation of the 4-momenta of a two-body system in their CM frame
//
// Here masses the two masses of the particles are all thats needed to specify the system
//
// Author:       Daniel Winney (2019)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "four_momenta.hpp"

namespace jpacPhoto
{
    // ---------------------------------------------------------------------------
    // The four momentum of the meson

    complex four_momenta::q(lorentz_index mu, double s, double theta)
    {
        if (updated(s)) recalculate(s);

        switch (mu)
        {
            case t: return _meson_energy;
            case x: return _momentum * sin(theta);
            case y: return 0.;
            case z: return _momentum * cos(theta);
        };
    };

    // ---------------------------------------------------------------------------
    // The four momenta of the baryon, goes backwards compared to q
    
    complex four_momenta::p(lorentz_index mu, double s, double theta)
    {
        if (updated(s)) recalculate(s);

        switch (mu)
        {
            case t: return _baryon_energy;
            case x: return - _momentum * sin(theta);
            case y: return 0.;
            case z: return - _momentum * cos(theta);
        };
    };

    // -----------------------------------------------------------------------
    // Save a cache of the energy dependence to avoid recalculating so much

    void recalculate(double s)
    {
        _momentum      = sqrt(Kallen(s - IEPS, _mM2, _mB2)) / (2.*sqrt(s));
        _meson_energy  = (s + _mM2 - _mB2) / (2. * sqrt(s));
        _baryon_energy = (s - _mM2 + _mB2) / (2. * sqrt(s));
        _cached_s      = _s;
    };
};
