// Abstract class for a regge trajectory. Any user-defined class for a specific
// for of the RT can be used in amplitudes (e.g. pomeron_exchange) as long as it
// is derived from this class
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _REGGE_TRAJ_
#define _REGGE_TRAJ_

#include <complex>
#include <vector>
#include <string>

class regge_trajectory
{
    public:

    // constructor
    regge_trajectory(std::string name = "")
    : _parent(name)
    {};

    regge_trajectory(int sig, std::string name = "")
    : _signature(sig), _parent(name)
    {};

    // copy constructor
    regge_trajectory(const regge_trajectory & old)
    : _parent(old._parent), _signature(old._signature)
    {};

    virtual ~regge_trajectory() = default;

    // Only need a function to evaluate the trajectory at some s
    virtual std::complex<double> eval(double s) = 0;

    virtual std::complex<double> slope(double s = 0.){return 0.;};

    virtual std::complex<double> intercept(){return 0.;};

    virtual void set_params(std::vector<double> pars){};

    // These parameters define the trajectory
    // name, spin, and mass of the lowest lying resonance on the parent trajectory
    std::string _parent;
    int _signature;

    void set_minJ(int j) { _minJ = j; };
    int _minJ = 0;
};


// Basic linear regge_trajectory
class linear_trajectory : public regge_trajectory
{
    private:

    // Intercept and slope
    double _a0, _aprime;

    public:

    // Empty constructor
    linear_trajectory(){};

    // Parameterized constructor
    linear_trajectory(int sig, double inter, double slope, std::string name = "")
    : regge_trajectory(sig, name),
      _a0(inter), _aprime(slope)
    {};

    // Parameterized constructor with minJ also
    linear_trajectory(int sig, int minJ, double inter, double slope, std::string name = "")
    : regge_trajectory(sig, name),
      _a0(inter), _aprime(slope)
    {
        set_minJ(minJ);
    };

    // Setting utility
    inline void set_params(std::vector<double> pars){ if (pars.size() == 2) set_params(pars[0], pars[1]); };
    inline void set_params(double inter, double slope)
    {
        _a0 = inter; _aprime = slope;
    };

    inline std::complex<double> eval(double s)
    {
        return _a0 + _aprime * s;
    };

    inline std::complex<double> slope(double s = 0.)
    {
        return _aprime;
    };

    inline std::complex<double> intercept(){return _a0;};
};

#endif
