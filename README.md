#   jpacPhoto
Framework for amplitude analysis involving single meson production via quasi-elastic scattering on a nucleon target. Focus on expandability and easy interfacing with other libraries / analysis code. 

<p align="center">
  <img width="300" src="./doc/FeynmanDiagram.png">
</p>


Compilation of the base library requires only [ROOT](https://root.cern.ch/) (tested with versions 6.17 and 6.24) with [*MathMore*](https://root.cern.ch/mathmore-library) and [Boost C++](https://www.boost.org/) (version $\geq$ 1.68) libraries.

##  INSTALLATION
To install, clone normally and use:
```bash
cd jpacPhoto
mkdir build && cd build
cmake ..
cmake --build . --target install
```
This will create the core library `/lib/libJPACPHOTO.so` with the linkable library as well as ROOT dictionary (.pcm) files. 

Additionally a [scripting executable](./src/cling/jpacPhoto.cpp) will be installed into `/bin/jpacPhoto` which short-cuts loading the libraries into an interactive ROOT session and running a .cpp file as an interpreted script.   This executable requires the environment variable `JPACPHOTO` to be set to the top-level directory in order to find auxilary files. This can be done as such:
```bash
export JPACPHOTO=/path/to/jpacPhoto # for bash
setenv JPACPHOTO /path/to/jpacPhoto # for csh
```

##  USAGE

The primary use case is to reproduce results from recent [[1-5]](#references) and older analyses [[6-8]](#references) from the JPAC collaboration in a general and unified framework for amplitude analysis in photoproduction. The compiled library contains the framework of amplitudes and observables as abstract classes with specific amplitude models needed to be imported at run-time as header-only libraries. 

The compiled jpacPhoto executable pipes an analysis script, relevent amplitude files, and the compiled library into ROOT's cling interpeter to run. 
This set up mimics a Python-like environment without requiring recompilation of the whole library when changes are made to amplitude files. Amplitudes and scripts relevant for JPAC papers are located in [`/physics`](./physics) and [`/scripts`](./scripts) respectively. To run a script located in the bin directory simply run 
```bash
jpacPhoto my_script.cpp
```
or add the bin directory to $PATH to call `jpacPhoto` from any directory. 

The linkable library can be added to an existing CMake project using `find_library()` and linked as normal:
```cmake
find_library(JPACPHOTO NAMES JPACPHOTO libJPACPHOTO
                       HINTS "$ENV{JPACPHOTO}/lib")
target_link_libraries( myTarget JPACPHOTO)
```

###  AMPLITUDES
The main object of interest in the core library is the abstract [`amplitude`](./src/amplitude.hpp) and derived implementations defined by the user for specific physics models. All amplitudes available so far may be found in [/physics](./physics) as well as a [template file](./physics/template.hpp) to guide adding new classes. Amplitudes are calculated on a per-helicity basis which allows one to compute an array of observables (units of GeV and nb assumed where appropriate):

| Observable                           |                                                                | Callable `amplitude` function                                                                             
|--------------------------------------|----------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------|
| Unpolarized probability distribution | $\sum_{\{\lambda\}} \|A_{\{\lambda\}}\|^2$                     | `probability_distribution(double s, double t)`                                                                          | 
| Differential cross section           | $d\sigma/dt$                                                   | `differential_xsection(double s, double t)`                                                                             |
| Polarized differential cross section | $d\sigma_{\perp,\parallel}/dt$                                 | `polarized_dxsection(double pm, double s, double t)`                                                                    |
| Integrated total cross section       | $\sigma$                                                       | `integrated_xsection(double s)`                                                                                         |
| Double polarization asymmetries      | $A_{LL},K_{LL}$                                                | `A_LL(double s, double t)` <br /> `K_LL(double s, double t)`                                                            |
| Meson spin density matrix elements   | $\rho^{\alpha}_{\lambda,\lambda^\prime}$                       | `mSDME_H(int a, int lam, int lamp, double s, double t)` <br /> `mSDME_GJ(int a, int lam, int lamp, double s, double t)` |
| Baryon spin density matrix elements  | $\rho^{\alpha}_{\frac{\lambda}{2},\frac{\lambda^\prime}{2}}$   | `bSDME_H(int a, int lam, int lamp, double s, double t)` <br /> `bSDME_GJ(int a, int lam, int lamp, double s, double t)` |
| Beam asymmetry                       | $\Sigma_{4\pi}$                                                | `beam_asymmetry_4pi(double s, double t)`                                                                                |

All kinematics are passed around by the [`kinematics`](./src/kinematics.hpp) class which allows arbitrary masses for all particles and arbitrary quantum numbers for the produced final state meson & baryon. Multiple amplitudes may describe the same process and share the same kinematics instance.

The basic usage is:
```c++
// Set up kinematics
kinematics myKin = new_kinematics(M_MESON, M_BARYON);
myKin->set_meson_JP(1, -1);  // J = 1  , P = -1
myKin->set_baryon_JP(1, 1);  // J = 1/2, P =  1

// Set up amplitude
amplitude myAmp1 = new_amplitude<my_implementation>(myKin, /* additional parameters */);
myAmp1->set_parameters{ {/* couplings etc */} };

amplitude myAmp2 = new_amplitude<my_other_implementation>(myKin, /* additional parameters */);
myAmp2->set_parameters{ {/* couplings etc */} };

// Evaluate observables
myAmp1->integrated_xsection(s);
myAmp2->mSDME_GJ(0, 1, -1, s, t);
```
Interfering sums of amplitues and/or partial wave projections can be constucted out of any combination of compatible amplitude instances. These are all treated as amplitudes themselves and have access to all observables with the same syntax:
```c++
// Two different amplitudes but share same kinematics
amplitude amp1, amp2;

// Sum together
amplitude sum = amp1 + amp2;
sum->set_parameters{ /* couplings for both amp1 and amp2 */};

// Take the J = 1 partial wave
amplitude pwave = project(1, sum);

// Access observables
amp1->integrated_xsection(s);  // Individual term
sum->integrated_xsection(s);   // Interfering sum
pwave->integrated_xsection(s); // Only P-wave contribution of sum
```

### SEMI-INCLUSIVE DISTRIBUTIONS
Methods for semi-inclusive processes can be added via the [`semi_inclusive`](./src/semi_inclusive.hpp) class. These are used for example in [[3-4]](#references) to investigate inclusive XYZ production. 
Because semi-inclusive models are implemented at the cross section level and lose explicit helicity dependence, only unpolarized observables are available (Units of GeV and nb assumed where appropriate):
| Observable                                       |   | Callable `semi_inclusive` function |
|--------------------------------------------------|---|---------------------------------------|
| Lorentz-invariant cross section | $E_{\mathcal{Q}} \frac{d^3\sigma}{d^3\mathbf{q}}$  | `invariant_xsection(double s, double t, double M2)`  |
| Double-differential cross section | $\frac{d^2\sigma}{dt  dM^2}$, <br /> $\frac{d^2\sigma}{dt dx}$, <br /> $\frac{d^2\sigma}{dx dy^2}$,  | `dsigma_dtdM2(double s, double t, double M2)` <br />  `dsigma_dtdx(double s, double t, double x)`<br /> `dsigma_dxdy2(double s, double x, double y2)` |
| Single-differential cross section | $\frac{d\sigma}{dt}$, <br /> $\frac{d\sigma}{dM^2}$, <br /> $\frac{d\sigma}{dx}$, <br /> $\dots$ | `dsigma_dt(double s, double t)` <br />  `dsigma_dM2(double s, double M2)`<br /> `dsigma_dx(double s, double x)` <br /> $\dots$ |
| Fully integrated cross section | $\sigma$ | `integrated_xsection(double s)`|

Much of the syntax is the same as with exclusive amplitudes although models are implemented at the level of amplitude-squared. 
``` c++
// Use same interface with amplitudes
kinematics kZ = new_kinematics(M_Z);

// Semi-inclusive Z production
semi_inclusive Z_inc = new_semi_inclusive<inclusive_model>(kZ, /* other parameter */);
Z_inc->set_parameters({ /* couplings, etc */ });

// Some processes we want to add the exclusive amplitude
amplitude Z_exc = new_amplitude<exclusive_model>(kZ, /* etc */);
Z_inc += Z_exc; // Add to the purely inclusive distribution

// Get observables
Z_inc->integrated_xsection(s);
```
### ANALYSIS TOOLS
Tools to fit amplitudes to experimental data are available through the [`fitter`](./src/fitter.hpp) and [`plotter`](./src/plotter.hpp) classes. Data may be imported using the [`data_set`](./src/data_set.hpp) class as interface. Arbitrarily many data sets may be imported into a fitter where one must specify the minimazition function per data type. An end-to-end example used in [[4]](#references) may be found in the appropriate [scripts directory](./scripts/jpsi_p/fit.cpp).

A schematic analysis may look like:
```c++
struct my_analysis
{
    // Strings to identify data types in fitter's terminal messages
    static std::string data_type(int i)
    {
        if (i == 0) "Allowed type";
        else        "Not allowed type";
    };

    // Minimization funcion
    static double fcn(const std::vector<data_set> &data, amplitude amp)
    {
        double chi2 = 0;
        for (auto datum : data) if (data._type == 0) chi2 += /* calculate chi2*/;
        return chi2;
    };
};

int main()
{
  // Set up amplitude
  amplitude my_amplitude;

  // Do fit
  fitter fitter<my_analysis>(my_amplitude);
  fitter.add_data(my_data);
  fitter.set_parameter_limits("parameter 1", {0., 1});
  fitter.fix_parameter("fixed parameter", 0.5);
  fitter.sync_parameter("parameter 2", "parameter 3"); // enforce 2 == 3 with 3 free
  fitter.do_fit();

  // Plot results
  plotter plotter;
  plot p = plotter.new_plot();
  p.add_data(my_data);
  p.add_curve({min, max}, [&](double x){ return my_amplitude->observable(x); }, "label");
  p.save("outfile.pdf");
};

```

### AmpTools
In addition to the built-in analysis tools, the aim of this project is to be usable with existing tools. A relatively simple interface with [AmpTools](https://github.com/mashephe/AmpTools) is provided in [/AmpTools](./AmpTools) as well as scripts in the analogous [/scripts](./scripts/AmpTools/) directory. More complicated processes involving unstable particles can be incorporated by wrapping the jpacPhoto helicity amplitudes into AmpTools classes which decay lineshapes. 

To run these scripts with the `jpacPhoto` executable requires having AmpTools to be built and the `AMPTOOLS` environment variable set to the top level install directory. Since the cling interpreter can only load _dynamic_ libraries (which is not the default installation mode of AmpTools), the optional cmake flag `-DDYNAMIC_AMPTOOLS=TRUE` is provided to repackage an existic static library to one shared version. 

###  REFERENCES
+ [1] [Double Polarization Observables in Pentaquark Photoproduction](https://arxiv.org/abs/1907.09393)
+ [2] [XYZ spectroscopy at electron-hadron facilities: Exclusive processes](https://arxiv.org/abs/2008.01001)
+ [3] [XYZ spectroscopy at electron-hadron facilities II: Semi-inclusive processes with pion exchange](https://arxiv.org/abs/2209.05882)
+ [4] [XYZ spectroscopy at electron-hadron facilities III: Semi-inclusive processes with vector exchanges](https://arxiv.org/abs/2404.05326)
+ [5] [Dynamics in near-threshold J/ψ photoproduction](https://arxiv.org/abs/2305.01449)
+ [6] [Features of πΔ Photoproduction at High Energies](https://arxiv.org/abs/1710.09394)
+ [7] [Vector Meson Photoproduction with a Linearly Polarized Beam](https://arxiv.org/abs/1802.09403)
+ [8] [Exclusive Tensor Meson Photoproduction](https://arxiv.org/abs/2005.01617)
+ [9] [JPAC Website](https://www.jpac-physics.org/)

<p align="center">
  <img width="275" src="./doc/JPAClogo.png">
</p>
