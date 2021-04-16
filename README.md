#   jpacPhoto
Framework for amplitude analysis involving single meson production via quasi-elastic scattering of a real photon on a nucleon target. Focus on expandability and easy interfacing with Monte-Carlo tools and event generators.

<p align="center">
  <img width="400" src="./doc/FeynmanDiagram.png">
</p>

Such processes are of interest at many experiments at JLab and the future EIC.

Base library requires only [ROOT](https://root.cern.ch/) (tested with version 6.17) with [*MathMore*](https://root.cern.ch/mathmore-library) libraries installed. 

Optional dependencies include [Boost C++](https://www.boost.org/) libraries for the extension jpacBox library to compute one-loop box diagrams and the [jpacStyle](https://github.com/dwinney/jpacStyle) library to build executables which reproduce plots in e.g. to reproduce plots in [[1]](https://arxiv.org/abs/1907.09393) and [[2]](https://arxiv.org/abs/2008.01001).

##  INSTALLATION
To install clone normally and use:
```bash
mkdir build && cd build
cmake ..
cmake --build . --target install
```
This will create a `jpacPhoto/lib/jpacPhotolib.so` with the linkable library. If Boost is found in PATH, `jpacBoxlib.so` will also be built.


To build the suite of executables the [jpacStyle](https://github.com/dwinney/jpacStyle) library must be installed with environment variables set as such:
```bash
# for bash
export JPACSTYLE=/path/to/jpacStyle

# for csh
setenv JPACSTYLE /path/to/jpacStyle
```

##  AMPLITUDES
The main object of interest is the abstract [`amplitude`](./include/amplitudes/amplitude.hpp) class. Reactions are calculated on a per helicity amplitude basis  which allows one to compute an array of [observables](./src/amplitudes/observables.cpp):

* Probability distribution ( Σ_λ | A |^2 )
* Differential cross section ( dσ / dt )
* Integrated total cross section ( σ )
* Polarization asymmetries ( A_LL and K_LL )
* Spin density matrix elements ( ρ^α_λ,λ' )
* Integrated beam asymmetry ( Σ_4pi )
* Beam asymmetry in the y-direction ( Σ_y )
* Parity asymmetry ( P_σ )

All kinematics are passed around by the `reaction_kinematics` class which allows all masses to float, allowing amplitudes to handle non-elastic processes such as ψ p -> D Λc with minimal change.

Available amplitudes, so far, include:

### s-channel:
* [Baryon resonance](./include/amplitudes/baryon_resonance.hpp) - following [[1]](https://arxiv.org/abs/1907.09393), Breit-Wigner resonance for a baryon of spin up to 5/2 and arbitrary parity. 
 
### t-channel:
* [Pomeron exchange](./include/amplitudes/pomeron_exchange.hpp) - three parameterizations of Pomeron exchange are available: [helicity-conserving](https://arxiv.org/abs/1606.08912), [exponential-pomeron](https://arxiv.org/abs/1907.09393), and [dipole-pomeron](https://arxiv.org/abs/1508.00339).
* [(fixed-spin and reggeized) Charged pseudo-scalar meson exchange](./include/amplitudes/pseudoscalar_exchange.hpp) - allows production of mesons with pseudo-scalar, vector, and axial-vector quantum numbers
* [(fixed-spin and reggeized) Vector meson exchange](./include/amplitudes/vector_exchange.hpp) - allows production of mesons with scalar, pseudo-scalar, vector, and axial-vector quantum numbers
* [Primakoff effect off nuclear target](./include/amplitudes/primakoff_effect.hpp) - special amplitude for investigation of X(3872) production via Primakoff effect

### u-channel:
* [(fixed-spin) Dirac fermion exchange](./include/amplitudes/dirac_exchange.hpp) - allows production of mesons with pseudo-scalar and vector quantum numbers
* [(fixed-spin) Rarita-Schwinger fermion exchange](./include/amplitudes/rarita_exchange.hpp) - allows production of mesons with pseudo-scalar and vector quantum numbers

Incoherent (interfering) sums of amplitudes may be constructed through the [`amplitude_sum`](./include/amplitudes/amplitude_sum.hpp) class.

##  BOX AMPLITUDE
The optional library `jpacBox` allows exchange amplitudes (in the t and u-channels) to be combined into a [box diagram](./include/box/box_amplitude.hpp) of the form:

<p align="center">
  <img width="400" src="./doc/BoxDiagram.png">
</p>

The calculation is done via a dispersion relation and integrating over the entire intermediate phase-space. The `box_amplitude` class requires the `gauss_kronrod` integration method from Boost C++ which can natively handle complex integrands and thus makes it particularly efficient in computing dispersion relations. 

##  REFERENCES
+ [1] [Double Polarization Observables in Pentaquark Photoproduction](https://arxiv.org/abs/1907.09393)
+ [2] [XYZ spectroscopy at electron-hadron facilities: Exclusive processes](https://arxiv.org/abs/2008.01001)
+ [3] [JPAC Website](http://cgl.soic.indiana.edu/jpac/index.php)

<p align="center">
  <img width="275" src="./doc/JPAClogo.png">
</p>