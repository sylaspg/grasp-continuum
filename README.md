# GRASPC - GRASP package adapted for the generation of continuum orbitals wave functions

This repository, GRASPC, is a fork of
**The General-purpose Relativistic Atomic Structure Package (GRASP)** -
a set of Fortran 90 programs for performing fully relativistic electron structure calculations of atoms.

[Link to GRASP repository](https://github.com/compas/grasp)


## About GRASPC

![GitHub Release](https://img.shields.io/github/v/release/sylaspg/grasp-continuum)
![GitHub Release Date](https://img.shields.io/github/release-date/sylaspg/grasp-continuum)
![GitHub last commit](https://img.shields.io/github/last-commit/sylaspg/grasp-continuum)

GRASPC is a fork of GRASP, intended to incorporate the scattered (continuum) electron calculations
in elastic (single-channel) and inelastic (multi-channel) processes.
> **Please note:** Currently, only elastic processes are implemented.

[Link to GRASPC paper in
Computer Physics Communications](https://www.sciencedirect.com/science/article/pii/S0010465525001936)

The main idea is to use as many computational apparatus as they are implemented in GRASP,
by adapting them to calculate continuum orbitals wave functions.

GRASPC is entirely transparent for usual calculations in GRASP (bound states and their properties).
Only when calculations involving continuum orbital are requested (during the execution of the `rmcdhf` program), the default flow is changed, and different outputs are produced.
The only differences are listed below:
- the default number of points in the radial (NNNP) is increased to 5000 in the `src/lib/libmod/parameter_def_M.f90` file, since continuum orbitals have to be calculated on the far-reaching computational grid;
for larger systems, it is recommended to magnify the grid size even several times
- `rmcdhf`: in addition to the other properties, $\langle r^3\rangle$, $\langle r^5\rangle$ and $\langle r^7\rangle$ are calculated for each orbital, since they are part of the model polarization potential
- `rwfnplot`: generation of the input file for the `gnuplot` plotting program has been added
- `CMakeList.txt`:
The file has been reorganized so that the additional user settings in the `CMakeList.user` file are correctly taken into account
- `CMakeList.user`: file with `-fallow-argument-mismatch` flag, required for recent versions of the `gfortran` compiler; this file is not part of the original GRASP repository

All modifications in the original source files are marked in the following way:
```
! PS      (or # PS in CMakeList.* )
  ... code modified/added
! PS END  (or # PS END in CMakeList.* )
```


## Theoretical background

### Generation of the continuum orbital wave function

The basics of the relativistic multiconfiguration Dirac-Hartree-Fock method (**RMCDHF**)
applied to scattering are described in [1].
In short, the scattering system is constructed as $N+1$ electron system,
where $N$ electrons are bound (with discrete, negative energy levels),
and one electron is from the continuum spectrum, with given (positive) energy $\epsilon$ and quantum number $\kappa$.

To obtain the large $P_{\kappa\epsilon}(r)$ and small $Q_{\kappa\epsilon}(r)$ components of the continuum orbitals,
the Dirac-Hartree-Fock equations

$\big(\frac{\textrm{d}}{\textrm{d}r}+\frac{\kappa}{r}\big)P_{\kappa\epsilon}=\big(\frac{2}{\alpha}+\alpha(\epsilon-V-V_{pol})\big)Q_{\kappa\epsilon}-X^{(Q)}$

$\big(\frac{\textrm{d}}{\textrm{d}r}-\frac{\kappa}{r}\big)Q_{\kappa\epsilon}=-\alpha\big(\epsilon-V-V_{pol}\big)P_{\kappa\epsilon}+X^{(P)}$

are solved using the _outward_ integration method implemented in GRASP. In the above equations,
$\alpha$ is the fine structure constant,
$\epsilon$ is the energy of the scattered particle,
$V = V(r)$ is the Coulomb potential,
$V_{pol} = V_{pol}(r)$ is a polarization potential,
and $X^{(P)}$ and $X^{(Q)}$ are the exchange terms [2].

> **Please note:** From the perspective of the GRASP source code, bound orbital energies are _positive_.
Thus, the continuum electron energy is considered _negative_.
This does not affect the user's actions, for him the continuum electron energy is always positive (or zero in a special case).

If the calculated wave function is to be coupled
with some other function (e.g., the bound state),
it should be normalized first.
The currently implemented _per energy_ continuum wave function normalization
is described in [3].

> **Please note:** For calculations of phase shifts and scattering lengths,
> normalization is usually not required.

### Polarization potential

At greater distances, the polarization potential coincides with dipole, quadrupole, and octupole polarization but is limited near the nucleus.
It is currently modeled and implemented in the following form:

$V_{pol}\left(r\right)=-\frac{1}{2}\frac{\alpha_d r^2}{\left(r^3+\langle r_0^3\rangle\right)^2}
-\frac{1}{2}\frac{\alpha_q r^4}{\left(r^5+\langle r_0^5\rangle\right)^2}
-\frac{1}{2}\frac{\alpha_o r^6}{\left(r^7+\langle r_0^7\rangle\right)^2}$,

where $\alpha_d$, $\alpha_q$ and $\alpha_o$ represent the static dipole, quadrupole and octupole polarizabilities, respectively;
$\langle r_0^3\rangle$, $\langle r_0^5\rangle$ and $\langle r_0^7\rangle$ are the cut-off parameters.

For atoms, a great source for static dipole polarizabilities is [4]
(except for the Livermorium atom, atomic number 116).
There is no such a similar source for static quadrupole and octupole polarizabilities, unfortunately.

For atoms, cut-offs $\langle r_0^3\rangle$, $\langle r_0^5\rangle$ and $\langle r_0^7\rangle$ can be taken
from bound states calculations, assuming that $\langle r_0\rangle$ is the radius
of the outermost orbital of the target atom.

Polarization potential can also be provided
in numerical form.
It has to be provided in a text file named `vpol` located in the working directory, containing space-separated pairs $r\ \ V_{pol}(r)$ in rows
(in Bohr radii $a_0$, and hartree units, respectively), ordered by increasing $r$ values, e.g.
```
1.00000E-05 -4.10660E-07
1.05127E-05 -4.31715E-07
1.10517E-05 -4.76698E-07
(...)
```

### Phase shift

Relativistic phase shift $\delta_l$ is obtained by comparing
the numerical solution with the asymptotic one at
large $r$ [1],
determined to the nearest $2\pi$:

$P_\kappa(r)/r = j_l(kr)\cos(\delta_l) - n_l(kr)\sin(\delta_l)$,

where $j_l(kr)$ and $n_l(kr)$ are the spherical Bessel functions
of the first and the second kind, respectively,
and $k=\sqrt{2\epsilon+\alpha^2\epsilon^2}$ is the wavenumber.
Orbital angular momentum quantum number $l$ is related to $\kappa$ such that $\kappa = l$ or $\kappa = -l-1$, whereby $\kappa \neq 0$.
Therefore, for the $l$-th partial wave (except the $l=0$ case) we obtain two phase shifts, sometimes referred to $\delta_l^+$ and $\delta_l^-$.

>**Please note:**
> To correctly determine the phase shift, the continuum orbital
> has to be calculated far enough from the origin,
> to ensure that atom-electron interaction can be neglected.

>**Please note:**
> The phase shift is determined to the
nearest $2\pi$, and in many scientific publications is given in the $[−\pi, \pi]$ range.
> Thus, the calculated value is shifted by $2\pi$ if necessary, giving the second,
equivalent value in the output

### Scattering length

Scattering length is one of the most useful parameters for describing low-energy electron-atom collisions.
It is defined as the radius of a rigid sphere in the zero-energy total cross-section.
The sign of a scattering length represents the type of interaction:
positive for repulsion and negative for attraction.

The typical method for determining the scattering length $a$ involves
analyzing the asymptotic behavior of the wave function:

$a = -\lim_{k \to 0}\frac{\tan(\delta_0)}{k}$,

where $\delta_0$ is the phase shift for $l=0$ (_s_-wave, $\kappa=-1$), calculated for (energy-dependent) wavenumber tending towards zero in the limit.

In addition to the above method
(which is subject to the inaccuracy associated with the non-zero energy),
the calculation of the scattering length is implemented as the intersection
of the asymptote of the zero-energy ($\epsilon=0$) wave function with the $r$-axis.
The details of that approach are described in [5].

> **Please note:**
> _Zero energy_ wave function, used for accurate scattering length calculations
> is not a typical oscillating wave function.

More detailed explanations of the theoretical background
and computational procedure
can be found in [6].

### References

1. P. Syty and J.E. Sienkiewicz, Relativistic Multiconfiguration Dirac-Hartree-Fock in the scattering of electrons from argon atoms,
_J. Phys. B: At. Mol. Opt. Phys._ 38 (2005) 2859. https://doi.org/10.1088/0953-4075/38/16/001
2. I.P. Grant, B.J. McKenzie, P.H. Norrington, D.F. Mayers, and N.C. Pyper. An atomic multiconfigurational Dirac-Fock package. Comput. Phys. Commun. 21 (1980) 207–231. https://doi.org/10.1016/0010-4655(80)90041-7
3. R.D. Cowan, The Theory of Atomic Structure and Spectra,
_University of California Press, Berkeley_ (1981) pp. 516–521. https://www.ucpress.edu/book/9780520038219/the-theory-of-atomic-structure-and-spectra
4. P. Schwerdtfeger and J.K. Nagle,
2018 Table of static dipole polarizabilities of the neutral elements in the periodic table,
_Molecular Physics_ 117 (2019) 9-12. https://doi.org/10.1080/00268976.2018.1535143
5. P. Syty, M.P. Piłat, J.E. Sienkiewicz,
Calculation of electron scattering lengths on Ar, Kr, Xe, Rn and Og atoms,
_J. Phys. B: At. Mol. Opt. Phys._ 57  (2024) 175202.
https://doi.org/10.1088/1361-6455/ad4fd1
6. P. Syty, M. Piłat, J.E. Sienkiewicz, _GRASPC – GRASP package adapted for the generation of continuum orbitals wave functions_, Comput. Phys. Commun. 315 (2025) 109691.
https://doi.org/10.1016/j.cpc.2025.109691


## Current status

### Implemented and tested
- Calculations of the continuum wave function of electrons _elastically_ scattered from atoms and ions, with the model or numerical polarization term
- Phase shift calculations (with grid control to ensure correct results)
- Calculations of electronic scattering lengths using the _zero energy_ wave function
- Normalization of the calculated continuum orbital

### Changelog

- **2025-03-13**
  - `co_add_polarization_potential` Added: Octupole term in the model polarization potential (this entailed a slight change in the input reception)
  - `co_scattering_length` Added: Calculation of relative error of straight line determination at the tail, in scattering length calculation (in the zero-energy case)
  - Several minor code and test-case improvements

- **2024-07-18**
  - `rmcdhf` Added: phase shift calculation
  - `rwfnestimate` Added: Dedicated method for initial estimate of the radial wave function for continuum electron
  - `grasptest/continuum` Added: new example for electron-argon scattering (_d_-wave calculation)
  - `README` Updated: most sections, the most important updates for _User guide_ and _Theoretical background_
  - `rmcdhf` Now, the energy of the continuum electron should be entered as a positive value (or zero)
  - Several minor code and test-case improvements

- **2024-06-25** (Initial release)
  - Continuum orbitals generator with scattering length calculations, together with examples for _e-Ar_ and _e-Sr_ scattering


### TO-DO list

- Positron scattering (work in progress)
- Validation of user-input parameters (work in progress)
- **Inelastic scattering** (work at the final concept phase)


## User guide

1. **Bound states**

   Optimize bound states of the selected target (atom/ion)
   in a usual way (`rnucleus` => `rcsfgenerate` => `rwfnestimate` => `rangular` => `rmcdhf` => `rsave`
   in the simplest case), or take nuclear data (`isodata`)
   and radial wave functions (`rwfn.out` / `.w`) files from previous calculations.

   GRASPC may also be used for that calculation,
   since it works _exactly_ as the original GRASP for bound states.

    > **Please note:**
    > Since only radial orbitals are used (and not the expansion coefficients),
    > Calling the `rci` code does not improve anything.


2. **CSF list**

   By invoking `rcsfgenerate`, create a special `rcsf.inp` file with only one CSF,
   where 'core' is the configuration of the target atom or ion,
   and 'peel' consists of _one_ additional _inactive_ electron.
   This electron will be treated as a scattered one;  its principal quantum number will be ignored,
   and its quantum number $\kappa$ will be determined from the subshell designation
   and final $J$ value, resulting from coupling with the core.
    > **Examples for argon-electron scattering:**
    > - _Core_ subshells: `1s 2s 2p- 2p 3s 3p- 3p`,
    > - _Additional_ subshell:
    >   - `4s(1,i)`, coupled to $J=1/2$ means _s_-wave electron of $\kappa=-1$
    >   - `4p(1,i)`, coupled to $J=1/2$ means _p_-wave electron of $\kappa=1$
    >   - `4p(1,i)`, coupled to $J=3/2$ means _p_-wave electron of $\kappa=-2$
    >   - `4d(1,i)`, coupled to $J=3/2$ means _d_-wave electron of $\kappa=2$
    >   - `4d(1,i)`, coupled to $J=5/2$ means _d_-wave electron of $\kappa=-3$, etc.

   Remember to provide the $J$ value as $2J, 2J$ range (e.g. _1,1_ for $J=1/2$), and to provide _0_ as the number of excitations.

    Example `rcsf.inp` file for electron-argon scattering (_s_-wave):
  ```
Core subshells:
  1s   2s   2p-  2p   3s   3p-  3p
Peel subshells:
  4s
CSF(s):
  4s ( 1)
      1/2
       1/2+
  ```

3. **Angular coefficients**

   Run `rangular` as usual (with _Full interaction_ option enabled).

4. **Initial estimations of the radial wave functions**

   Run `rwfnestimate` using the previously calculated radial wave functions
   as initial estimation for the 'core' orbitals (option _1 -- GRASP92 File_).
   To estimate the additional electron designed to be a _continuum_ one, use the new option _5 -- Continuum orbital_.
   Any other method (options _2 - 4_) should also work
   (try them if you encounter convergence problems during the actual calculations in the next step).

5. **SCF procedure**

    Invoke `rmcdhf`, then
    - answer _n_ when asked about _Default settings_
    - answer _y_ when asked about _Perform continuum wave function calculations?_
    - provide continuum electron energy in hartree
    (should be positive or zero for accurate scattering length calculation)
    - decide if polarization potential should be included.
      - _0_ - do not include polarization potential
      - _1_ - include the dipole term of the model potential with default parameters:
            $\alpha_d$ taken from [2], and cut-off $\langle r_0^3\rangle$ taken from bound state calculations
            as the size of the outermost orbital; here, the quadrupole and octupole terms are omitted
      - _2_ - include model potential with all parameters provided manually by the user; in turn:
            $\alpha_d$, $\langle r_0^3\rangle$, $\alpha_q$, $\langle r_0^5\rangle$, $\alpha_o$ and $\langle r_0^7\rangle$
      - _3_ - include numerical potential from the file named `vpol`
    - decide (_y_/_n_) if the calculated continuum wave function should be normalized
    - answer _y_ when asked _Change default speed of light or radial grid parameters?_
    - answer _y_ when asked _Revise default radial grid parameters?_
    - enter new _RNT_ and _H_ values (firstly, the defaults might be kept)
    - enter new _HP_; this value is the step size ath the tail; use non-zero value to force the linearly-logarithmic grid,
    which ensures adequate grid density far from the scattering center;
    _1.0_ or less is a good choice for a first-try
    - enter new _N_; in general, use as a big number as possible to ensure as long grid as possible;
    _5000_ is the default, but even tens of thousands of points might be required in some cases, which would require code modification (_NNNP_ value in `src/lib/libmod/parameter_def_M.f90`) and recompilation

    > **Please note:**
    > If calculations do not converge (_Maximal iterations exceeded_),
    experiment with the other grid parameters,
    > and/or try the other method to estimate the radial wave function for the continuum electron.

    > **Please note:**
    > Only the `rmcdhf` code has been adapted for continuum orbital calculations, not
    `rmcdhf_mem`, `rmcdhf_mem_mpi` or `rmcdhf_mpi`.

6. **Output**

    The calculated continuum orbital wave function will be stored in the `rwfn.out` file
    (together with the bound orbitals), and also in a text-formatted file `continuum.csp`.

    If the grid is long enough and electron energy is not zero, calculated phase shift and scattering length
    are written to the screen and the `rmcdhf.sum` summary file. The accuracy of the scattering length strongly depends
    on the electron energy.
    (smaller energy means better accuracy since it should be calculated in the $k\to0$ limit).

    If the electron energy is set to zero, only the scattering length is calculated in a more accurate approach,
    as the intersection of the asymptote of the zero-energy wave function with the r-axis.
    In that case, additional information is given, which may be helpful for error estimation.
    The first one is the relative error of straight line determination at the tail.
    The second one is the percentage difference between the scattering length
    calculated from the last two points on the grid
    and the penultimate and the one before the penultimate point.
    In both cases, a lower value means better accuracy.

## Examples

### Elastic scattering of electrons from argon atoms

Bound states estimation, continuum orbital wave functions generation
(for two different partial waves), phase shift calculations, scattering length calculations using the _zero energy_ approach;
calculations with default, numerical and manually entered parameters of polarization potential.
The energies of the incident electron were chosen to compare the results with those in
Y. Cheng, S. Liu, S. B. Zhang, and Y.-B. Tang, Phys. Rev. A 102, 012824 (2020),
https://doi.org/10.1103/PhysRevA.102.012824 easily.

Files in `/grasptest/continuum/argon-electron_scattering` directory:
- `1_Ar_continuum_electron_s-wave_function` - script calculating continuum orbital of $\kappa = -1$ (_s_-wave, $J=1/2$) and electron energy $\epsilon=0.00124997$ hartree, using the dipole term of the polarization potential with default parameters
- `2_Ar_continuum_electron_d-wave_function`- script calculating _normalized_ continuum orbital of $\kappa = 2$ (_d_-wave, $J=3/2$) and electron energy  $\epsilon=0.04499880$ hartree, using polarization potential read from a file
- `3_Ar_electronic_scattering_length` - script calculating electronic scattering length using the zero-energy wave function
- `Ar_bound_DF` - auxiliary script creating the nuclear data and calculating the bound states of argon atom (very simple case without correlations); there's no need to call it by hand; it is invoked automatically by all of the above scripts
- `clean` - auxiliary script, deleting all of the output files created during program execution
- `vpol` - auxiliary file with the polarization potential given in numerical form; it is used by the second script above

Scripts realizing the main calculations (beginning with a number) can be run independently and in any order.

The results stay in very good agreement with the results from Cheng when it comes to phase shifts
in the first two examples ($\delta_0 = 0.03432$ vs. $0.03999$,  and $\delta_2^+ = 0.03241$ vs. $0.03253$),
and electronic scattering length in the third example ($a_0 = -1.320$ vs. $-1.394$).

### Electronic scattering length of strontium

Calculation of scattering length with default parameters of polarization potential,
for a given `.w` file with optimized bound states.

Files in the `/grasptest/continuum/strontium_electronic_scattering_length` directory:
- `Sr_scattering_length` - script calculating electronic scattering length using the zero-energy wave function
- `isodata, Sr.w` - nuclear data and previously calculated bound states of strontium
- `clean` - script cleaning the output files

## Contributors
#### Code development and testing; preparing of the documentation; preparing and scripting the test cases; maintaining the repository; co-authoring the publication
- Paweł Syty, Gdańsk University of Technology, pawel.syty@pg.edu.pl
#### Giving the ideas and theoretical background; proposing and improving the test cases; improving the documentation;
co-authoring the publication
- Józef E. Sienkiewicz, Gdańsk University of Technology
- Michał Piłat, Gdańsk University of Technology

## How to cite

P. Syty, M. Piłat, J.E. Sienkiewicz, _GRASPC – GRASP package adapted for the generation of continuum orbitals wave functions_,
Comput. Phys. Commun. 315 (2025) 109691.
https://doi.org/10.1016/j.cpc.2025.109691


## Installation

> **Please note:**
> All the installation instructions for the original GRASP
> are also valid for GRASPC, with the repository address changed to:
`https://github.com/sylaspg/grasp-continuum.git`.
Further in this section are original instructions adapted for the GRASPC repository. **All the credits goes to the GRASP contributors, who prepared the original installation procedure.**

To compile and install GRASPC, first clone this Git repository:

```sh
git clone https://github.com/sylaspg/grasp-continuum.git
```

There are two ways to build GRASPC: either via [CMake](https://cmake.org/) or via the
`Makefile`s in the source tree. Either works, and you will end up with the GRASPC binaries in the
`bin/` directory.

CMake is the recommended way to build GRASPC. The `Makefile`-based workflow is still there to
make it smoother to transition from `Makefile`s to a modern build system.

### CMake-based build

The first step with CMake is to create a separate out-of-source build directory. The
`configure.sh` script can do that for you:

```sh
cd grasp-continuum/ && ./configure.sh
```

This will create a `build/` directory with the default _Release_ build
configuration. However, `configure.sh` is just a simple wrapper around a `cmake`
call, and if you need more control over the build, you can always invoke `cmake`
yourself (see [CMake documentation](https://cmake.org/documentation/) for more
information).

To then compile GRASPC, you need to go into the out-of-source build directory and
call `make`:

```sh
cd build/ && make install
```

Remarks:

* Running `make install` instructs CMake to actually _install_ the resulting binaries into
  the conventional `bin/` directory at the root of the repository.

  When you run `make`, the resulting binaries will end up under the `build/` directory
  (specifically in `build/bin/`). This is useful when developing and debugging, as it allows
  you to compile many versions of the binaries from the same source tree with different
  compilation options (e.g. build with debug symbols enabled) by using several out of source
  build directories.

* With CMake, GRASPC also supports parallel builds, which can be enabled by passing the `-j`
  option to `make` (e.g., `make -j4 install` to build with four processes).

* The CMake-based build allows running the (non-comprehensive) test suite by calling `ctest`
  in the `build/` directory. The configuration and source files for the tests are under
  `test/`/

### `Makefile`-based build

The legacy `Makefile`-based build can be performed by simply calling the `make` in the top
level directory:

```sh
make
```

In this case, each library and program is compiled in their
respective directory under `src/`, and the built artifacts are stored in the source tree.
The resulting binaries and libraries will directly get installed under the `bin/` and `lib/`
directories.

To build a specific library or binary, you can pass the path to the source directory as the
Make target:

```sh
# build libmod
make src/lib/libmod
# build the rci_mpi binary
make src/appl/rci90_mpi
```

Note that any necessary library dependencies will also get built automatically.

**WARNING:** The `Makefiles' do not know about the dependencies between the source files, so
parallel builds (i.e., calling `make` with the `-j` option) do not work.

#### Customizing the build

By default, the `Makefile` is designed to use `gfortran`. The variables affecting GRASPC
builds are defined and documented at the beginning of the `Makefile`.

It should never be necessary for the user to modify the `Makefile` itself. Rather, a
`Make.user` file can be created next to the main `Makefile` where the build variables can be
overridden. E.g., to use the Intel Fortran compiler instead, you may want to create the
following `Make.user` file:

```make
export FC = ifort
export FC_FLAGS = -O3 -save -mkl=sequential
export FC_LD =
export FC_MPI = mpiifort
export OMPI_FC=${FC}
```
where `-mkl=sequential` should be set depending on what version of ifort you can access.

Alternatively, to customize the GNU gfortran build to e.g., use a specific version of the compiler, you can create a `Make.user` file such as
```make
export FC = gfortran-9
export FC_FLAGS = -O3 -fno-automatic
export FC_LD =
export FC_MPI= mpifort
export OMPI_FC=${FC}
```

To set up a linker search path for the BLAS or LAPACK libraries, you can
set `FC_LD` as follows:

```make
export FC_LD = -L /path/to/blas
```

The repository also contains the `Make.user.gfortran` and `Make.user.ifort` files, which can be used as templates for your own `Make.user` file.

