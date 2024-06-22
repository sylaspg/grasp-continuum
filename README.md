# GRASP - The General-purpose Relativistic Atomic Structure Package

![Tests][tests-badge]
[![][doxygen-badge]][doxygen-url]
[![][manual-badge]][manual-pdf]

The General-purpose Relativistic Atomic Structure Package (GRASP) is a set of Fortran 90
programs for performing fully-relativistic electron structure calculations of
atoms.

## About this fork

This fork is intended for incorporating the continuum electron calculations
in both elastic (single channel) and inelastic (multi-channel) processes into GRASP.
> **Please note:** Currently, only elastic processes are implemented.

The main idea is to use as many computational processes as they are implemented in GRASP,
and simply adapting them to calculate electron from continuum spectra.

This fork is completely neutral for 'normal' type of calculations (bound states and their properties).
Only when calculations involving continuum orbital are requested (during execution of the `rmcdhf` program),
the default flow changes. The only differences are listed below:
- `rmcdhf`: number of points in the radial (NNNP) has been increased to 5000
- `rmcdhf`: $\langle r^3\rangle$ and $\langle r^5\rangle$ are also calculated for each orbital
- `rwfnplot`: generation of input file for the `gnuplot` plotting program has been added
- `CMakeList.txt`: file has been reorganized, because
in the original file flags introduced in the `CMakeList.user` are not applied correctly
- `CMakeList.user`: flag `-fallow-argument-mismatch` has been specified,
to allow code compilation using recent versions of the GNU Fortran (`gfortran`) compiler

All modifications in the original source code are clearly marked in the following way:
```
! PS      (or # PS in CMakeList.* )
  ... code modified/added
! PS END  (or # PS END in CMakeList.* )
```


### Theory

_This section is a draft only. It will be updated soon._

#### Continuum orbital wave function generation

Basics od the Relativistic Multiconfiguration Dirac-Hartree-Fock method (RMCDHF)
applied to scattering are described in [1].
In short, scattering system is constructed
as $N+1$-electron system.
Then the Dirac-Fock equations are being solved
to obtain the large and small components
of the continuum electron wave function.

If calculated wave function is to be coupled
with the other function (e.g. the bound state),
it should be normalized.
For calculations of phases shifts and scattering lengths,
normalization is not required.
_Per energy_ continuum wave function normalization is described in [2].

#### Polarization potential

The polarization potential coincides with dipole polarization at greater distances
but limited near the nucleus.
For atoms, it is currently implemented in the following form:

$V_{pol}\left(r\right)=-\frac{1}{2}\frac{\alpha_d r^2}{\left(r^3+\langle r_0^3\rangle\right)^2}
-\frac{1}{2}\frac{\alpha_q r^4}{\left(r^5+\langle r_0^5\rangle\right)^2}
$,

where $\alpha_d$ and $\alpha_q$ represents the static dipole
and quadrupole polarizability, respectively;
$\langle r_0^3\rangle$ and $\langle r_0^5\rangle$ are the cut-off parameters.

One possible source for atomic dipole polarizabilities is [3],
except for Livermorium atom (Lv, atomic number 116).
There is no single aggregate source for atomic quadrupole polarizabilities.

For atoms, cut-offs $\langle r_0^3\rangle$ and $\langle r_0^5\rangle$ can be calculated
from bound states, assuming that $\langle r_0\rangle$ is the radius
of the outermost orbital of the target atom.

>**Please note:**
> Polarization potential can be also provided
> in numerical form, which may be useful e.g. for ions.
> It has to be provided in a text file named `vpol`,
>containing pairs of numbers 'r  Vpol(r)' written in lines, e.g.
```
1.00000E-05 -4.10660E-07
1.05127E-05 -4.31715E-07
1.10517E-05 -4.76698E-07
(...)
```

#### Scattering length

Scattering length is one of the most useful parameters for describing low-energy electron-atom collisions.
It is defined as the radius of a rigid sphere in the zero-energy total cross-section.
The sign of a scattering length represents the type of interaction:
positive for repulsion and negative for attraction.

Calculation of the scattering length is implemented as the intersection of the asymptote
of the zero-energy wave function with the r-axis [4].

> **Please note:**
> When scattering lengths calculations are performed,
> the resulting _zero energy_ wave function is not a typical,
> oscillating wave function.

##### References

- [1] P. Syty and J.E. Sienkiewicz, Relativistic Multiconfiguration Dirac-Hartree-Fock in scattering,
_J. Phys. B: At. Mol. Opt. Phys._ 38 2859 (2005), https://doi.org/10.1088/0953-4075/38/16/001
- [2] R.D. Cowan, The Theory of Atomic Structure and Spectra,
_University of California Press, Oakland_ pp. 522–524 (1981)
- [3] P. Schwerdtfeger and J.K. Nagle,
2018 Table of static dipole polarizabilities of the neutral elements in the periodic table,
_Molecular Physics_ 117 9-12 (2019), https://doi.org/10.1080/00268976.2018.1535143
- [4] P. Syty, M.P. Piłat, J.E. Sienkiewicz,
Calculation of electron scattering lengths on Ar, Kr, Xe, Rn and Og atoms,
_J. Phys. B: At. Mol. Opt. Phys._ (accepted manuscript),
https://doi.org/10.1088/1361-6455/ad4fd1


### Current status

Implemented and tested:
- Calculations of continuum wave function (s-wave, $\kappa = -1$)
  of low and very low energy electron elastically scattered from atoms and ions
  (with model or numerical polarization)
- Calculations of electronic scattering lengths using the _zero energy_ wave function
- Normalization of the calculated wave function

### TO DO

- Test of higher partial waves scattering
- Phase shifts determination
- Positron scattering
- Validation of user-input parameters
- **Inelastic scattering**


### User guide

1. Optimize bound states of the selected target (atom / ion)
   in a normal way (`rnucleus` => `rcsfgenerate` => `rwfnestimate` => `rangular` => `rmcdhf` => `rsave`
   in the simplest case), or just take nucleus properties (`isodata`)
   and radial wave functions (`rwfn.out` / `.w`) files from any previous calculations.

   This fork may be also used for that calculations,
   since it works _exactly_ as the original GRASP for bound states.

2. By invoking `rcsfgenerate` create a special `rcsf.inp` file with only one CSF,
   where 'core' is the configuration of the target atom or ion,
   and 'peel' consists of _one_ additional _inactive_ electron.
   This electron will be treated as continuum one;  its principal quantum number will be ignored,
   and its quantum number kappa will be determined from the subshell designation.
   As an example, for electron-Argon scattering 'core' subshells are
   `1s 2s 2p- 2p 3s 3p- 3p`, and that additional electron may be provided
   as `4s(1,i)`, which stands for one _inactive_ electron of $\kappa = -1$ (s-wave), and $J = 1/2$.
   Next, provide proper $2*J$ range (e.g. _1,1_ if scattering from neutral atom)
   and enter _0_ as number of excitations.

Example `rcsf.inp` file for electron-Argon scattering:
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

3. Run `rangular` as usual.
4. Run `rwfnestimate`, use the previously calculated radial wave functions
   as initial estimation for the 'core' orbitals (option _1 -- GRASP92 File_),
   and any onther method (options _2 - 4_) for initial estimation of the 'peel' continuum electron.
5. Invoke `rmcdhf`, then
    - answer _n_ when asked _Default settings?_
    - answer _y_ when asked _Perform continuum wave function calculations?_
    - provide continuum electron energy in hartree
    (should be negative according to convention used in GRASP,
    or zero for scattering length calculation)
    - decide if polarization potential should be used.
      - _0_ - do not use polarization potential
      - _1_ - use model potential with default parameters:
            $\alpha_d$ taken from [2], and cut-off $\langle r_0^3\rangle$ taken from bound state calculations
            as the size of the outermost orbital; with that option the quadrupole term is omitted
      - _2_ - use model potential with all parameters
            provided manually by the user; in turn:
            $\alpha_d$, $\langle r_0^3\rangle$, $\alpha_q$ and $\langle r_0^5\rangle$
      - _3_ - use numerical potential from file named `vpol`
    - decide (_y_/_n_), if the calculated continuum wave function should be normalized
    - answer _y_ when asked _Change default speed of light or radial grid parameters?_
    - answer _y_ when asked _Revise default radial grid parameters?_
    - enter new _RNT_ and _H_ values (firstly, they might be the same as defaults)
    - enter new _HP_; use non-zero value to force the linearly-logarithmic grid,
    which ensures adequate grid density far from the scattering centre;
    _1.0_ or less is the good choice for a first try
    - enter new _N_; in general, use as big number as possible to ensure as long grid as possible;
    _5000_ is the default

      Answers to the other questions should be obvious to any GRASP user.

    > **Please note:**
    > If calculations do not converge (_Maximal iterations exceeded_),
    experiment with the other grid parameters.

6. The calculated continuum orbital wave function will saved in the `rwfn.out` file,
and also in a text-formatted file `continuum.csp`.

    If scattering length is calculated, the result will be
    written to screen and to `rmcdhf.sum` file. Moreover, additional parameter 'diff'
    is calculated and shown, specifying the percentage difference between the scattering length
    calculated from the last two points on the grid
    and that calculated from the penultimate
    and the one before the penultimate point (lower value means better accuracy).

    Summary file, `rmcdhf.sum` will be supplemented by
    some additional info about perfomed calculations, and (in zero energy case) scattering length.

### Examples

#### Electron-Argon elastic scattering

See `/grasptest/continuum/argon-electron_scattering` directory.

#### Electronic scattering length of Strontium

See `/grasptest/continuum/strontium_electronic_scattering_length` directory.

### Authors

- Paweł Syty, Gdańsk University of Technology, pawel.syty@pg.edu.pl
- Michał Piłat, Gdańsk University of Technology
- Józef E. Sienkiewicz, Gdańsk University of Technology


## Installation

> **Please note:**
> The installation instructions here are for the _development version_ on the
> `master` branch.
>
> To install the _latest published release_ (2018-12-03), go to the
> ["Releases" page](https://github.com/compas/grasp/releases/tag/2018-12-03),
> download the tarball from there and refer to the instructions in the README in
> the tarball.

To compile and install GRASP, first clone this Git repository:

```sh
git clone https://github.com/compas/grasp.git
```

> **Please note:**
> All the installation instructions in that section are for the original GRASP.
They are also valid for that fork, just use the following for cloning the repository:
`
git clone https://github.com/sylaspg/grasp-continuum.git
`

There are two ways to build GRASP: either via [CMake](https://cmake.org/) or via the
`Makefile`s in the source tree. Either works and you end up with the GRASP binaries in the
`bin/` directory.

CMake is the recommended way to build GRASP. The `Makefile`-based workflow is still there to
make smoother to transition from `Makefile`s to a modern build system.

### CMake-based build

The first step with CMake is to create a separate out-of-source build directory. The
`configure.sh` script can do that for you:

```sh
cd grasp/ && ./configure.sh
```

This will create a `build/` directory with the default _Release_ build
configuration. However, `configure.sh` is just a simple wrapper around a `cmake`
call and if you need more control over the build, you can always invoke `cmake`
yourself (see [CMake documentation](https://cmake.org/documentation/) for more
information).

To then compile GRASP, you need to go into the out-of-source build directory and
simply call `make`:

```sh
cd build/ && make install
```

Remarks:

* Running `make install` instructs CMake to actually _install_ the resulting binaries into
  the conventional `bin/` directory at the root of the repository.

  When you run just `make`, the resulting binaries will end up under the `build/` directory
  (specifically in `build/bin/`). This is useful when developing and debugging, as it allows
  you to compile many versions of the binaries from the same source tree with different
  compilation options (e.g. build with debug symbols enabled) by using several out of source
  build directories.

* With CMake, GRASP also supports parallel builds, which can be enabled by passing the `-j`
  option to `make` (e.g. `make -j4 install` to build with four processes).

* The CMake-based build allows running the (non-comprehensive) test suite by calling `ctest`
  in the `build/` directory. The configuration and source files for the tests are under
  `test/`/

### `Makefile`-based build

The legacy `Makefile`-based build can be performed by simply calling the `make` in the top
level directory:

```sh
make
```

In this case, the compilation of each of the libraries and programs happens in their
respective directory under `src/` and the build artifacts are stored in the source tree.
The resulting binaries and libraries will directly get installed under the `bin/` and `lib/`
directories.

To build a specific library or binary you can pass the path to the source directory as the
Make target:

```sh
# build libmod
make src/lib/libmod
# build the rci_mpi binary
make src/appl/rci90_mpi
```

Note that any necessary library dependencies will also get built automatically.

**WARNING:** the `Makefile`s do not know about the dependencies between the source files, so
parallel builds (i.e. calling `make` with the `-j` option) does not work.

#### Customizing the build

By default the `Makefile` is designed to use `gfortran`. The variables affecting GRASP
builds are defined and documented at the beginning of the `Makefile`.

For the user it should never be necessary to modify the `Makefile` itself. Rather, a
`Make.user` file can be create next to the main `Makefile` where the build variables can be
overridden. E.g. to use the Intel Fortran compiler instead, you may want to create the
following `Make.user` file:

```make
export FC = ifort
export FC_FLAGS = -O3 -save -mkl=sequential
export FC_LD =
export FC_MPI = mpiifort
export OMPI_FC=${FC}
```
where `-mkl=sequential` should be set depending on what version of ifort you have access to.

Alternatively, to customize the GNU gfortran build to e.g. use a specific version of the compiler, you can create a `Make.user` file such as
```make
export FC = gfortran-9
export FC_FLAGS = -O3 -fno-automatic
export FC_LD =
export FC_MPI= mpifort
export OMPI_FC=${FC}
```

To set up a linker search path for the BLAS or LAPACK libraries you can
set `FC_LD` as follows:

```make
export FC_LD = -L /path/to/blas
```

The repository also contains the `Make.user.gfortran` and `Make.user.ifort` files, which can be used as templates for your own `Make.user` file.

## About GRASP

This version of GRASP is a major revision of the previous GRASP2K package by [P.
Jonsson, G. Gaigalas, J. Bieron, C. Froese Fischer, and I.P. Grant Computer
Physics Communication, 184, 2197 - 2203 (2013)][grasp2k-2013] written in FORTRAN
77 style with COMMON and using Cray pointers for memory management.  The present
version is a FORTRAN95 translation using standard FORTRAN for memory management.
In addition, COMMONS have been replaced with MODULES, with some COMMONS merged.
Some algorithms have been changed to improve performance for large cases and
efficiently.

The previous package, was an extension and modification of GRASP92 by [Farid
Parpia, Charlotte Froese Fischer, and Ian Grant. Computer Physics Communication,
94, 249-271 (1996)][grasp92-1996].

This version of GRASP has been published in:

> C. Froese Fischer, G. Gaigalas, P. Jönsson, J. Bieroń,
> "GRASP2018 — a Fortran 95 version of the General Relativistic Atomic Structure Package",
> Computer Physics Communications, 237, 184-187 (2018),
> https://doi.org/10.1016/j.cpc.2018.10.032

Development of this package was performed largely by:
|                           | email                         |
| ------------------------- | ------------------------------|
| Charlotte Froese Fischer  | cff@cs.ubc.ca                 |
| Gediminas Gaigalas        | Gediminas.Gaigalas@tfai.vu.lt |
| Per Jönsson               | per.jonsson@mau.se            |
| Jacek Bieron              | jacek.bieron@uj.edu.pl        |

Supporters include:
|                           | email                         |
| ------------------------- | ------------------------------|
| Jörgen Ekman              | jorgen.ekman@mah.se           |
| Ian Grant                 | ian.grant@maths.ox.ac.uk      |

The GitHub repository is maintained by:
|                           | email                         |
| ------------------------- | ------------------------------|
| Jon Grumer                | jon.grumer@physics.uu.se

Please contact the repository manager should you have any questions with regards
to bugs or the general development procedure. Contact the leading developer for
specific questions related to a certain code.

## Structure of the Package

The package has the structure shown below where executables, after successful
compilation, reside in the `bin` directory. Compiled libraries are in the `lib`
directory. Scripts for example runs and case studies are in folders under
`grasptest`. Source code is in the `src` directory and divided into applications
in the `appl` directory, libraries in the `lib` directory and tools in the
`tool` directory.

```
   |-bin
   |-grasptest
   |---case1
   |-----script
   |---case1_mpi
   |-----script
   |-----tmp_mpi
   |---case2
   |-----script
   |---case2_mpi
   |-----script
   |-----tmp_mpi
   |---case3
   |-----script
   |---example1
   |-----script
   |---example2
   |-----script
   |---example3
   |-----script
   |---example4
   |-----script
   |-------tmp_mpi
   |---example5
   |-----script
   |-lib
   |-src
   |---appl
   |-----HF
   |-----jj2lsj90
   |-----jjgen90
   |-----rangular90
   |-----rangular90_mpi
   |-----rbiotransform90
   |-----rbiotransform90_mpi
   |-----ris4
   |-----rci90
   |-----rci90_mpi
   |-----rcsfgenerate90
   |-----rcsfinteract90
   |-----rcsfzerofirst90
   |-----rdensity
   |-----rhfs90
   |-----rhfszeeman95
   |-----rmcdhf90
   |-----rmcdhf90_mpi
   |-----rmcdhf90_mem
   |-----rmcdhf90_mem_mpi
   |-----rnucleus90
   |-----rtransition90
   |-----rtransition90_phase
   |-----rtransition90_mpi
   |-----rwfnestimate90
   |-----sms90
   |---lib
   |-----lib9290
   |-----libdvd90
   |-----libmcp90
   |-----libmod
   |-----librang90
   |-----mpi90
   |---tool
```


## Program Guide and Compilation

The software is distributed with a practical guide to [GRASP2018 in PDF-format
(click here to download)][manual-pdf]. The guide, which is under Creative
Commons Attribution 4.0 International (CC BY 4.0) license, contains full
information on how to compile and install the package.


## Acknowledgements

This work was supported by the Chemical Sciences, Geosciences and Biosciences
Division, Office of Basic Energy Sciences, Office of Science, U.S. Department of
Energy who made the Pacific Sierra translator available and the National
Institute of Standards and Technology. Computer resources were made available by
Compute Canada.  CFF had research support from the Canadian NSERC Discovery
Grant 2017-03851.  JB acknowledges financial support of the European Regional
Development Fund in the framework of the Polish Innovation Economy Operational
Program (Contract No. POIG.02.01.00-12-023/08).


## Copyright & license

The code in this repository is distributed under the [MIT license](LICENSE).
The accompanying guide  "A practical guide to GRASP2018" is licensed separately
under [the CC-BY-4.0 (Creative Commons Attribution 4.0 International) license][cc-by].

[manual-pdf]: https://github.com/compas/grasp2018/releases/download/2018-12-03/GRASP2018-manual.pdf
[manual-badge]: https://img.shields.io/badge/manual-pdf-blue.svg
[doxygen-url]: https://compas.github.io/grasp/
[doxygen-badge]: https://img.shields.io/badge/documentation-doxygen-blue.svg
[tests-badge]: https://github.com/compas/grasp/workflows/Tests/badge.svg
[grasp92-1996]: https://doi.org/10.1016/0010-4655(95)00136-0
[grasp2k-2013]: https://doi.org/10.1016/j.cpc.2013.02.016
[cc-by]: https://creativecommons.org/licenses/by/4.0/legalcode
