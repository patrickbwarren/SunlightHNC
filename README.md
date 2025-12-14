## General purpose HNC code

The current release is version 1.13.

Version 1.13 allows dissimilar and non-additive sphere diameters in the RPM and hard sphere potentials.  
Version 1.12 adds direct calculation of HNC free energy (see documentation).  
Version 1.11 fixes another bug: missing mean-field contribution to thermodynamics for off-SYM RPM potential.  
Version 1.10 fixes a bug: missing contact virial pressure contribution in the HNC free energy.  
Version 1.9 brings the source code up to modern FORTRAN standards and corrects a problem with EXP in v1.8.  
Version 1.8 includes the MSA closure and various potentials for Slater smeared charges.  
Version 1.7 handles an arbitrary number of components (rather than a maximum of three, as for previous versions).

A small update was made to the documentation recently (Dec 2025)
because it became unclear whether it was really true that the HNC free
energy satisfies *a* = ∑<sub>i</sub> ρ<sub>i</sub> μ<sub>i</sub> −
*p*, where *p* is the virial pressure.  In an exact theory, this would
be an identity, and it appears to hold accurately in HNC for soft
potentials.  I thank Josh Robinson for drawing my attention to this.

The experimental module `oz_aux.py` makes a start at implementing the
road map below, by wrapping the initialisation stages and the HNC
solver, and providing support for RPM potentials.

#### Features

* FORTRAN 90 based, with example python driver scripts;
* HNC, MSA, RPA, and EXP closures;
* fast Ng solver, and Ng decomposition for long range potentials;
* multicomponent (arbitrary number of components);
* hard core (RPM-like) and soft core (DPD-like) potentials;
* full structural thermodynamics;
* fully open source.

See [PDF documentation](oz_doc.pdf "oz_doc.pdf") for details.

A related pure python code with a more limited feature set is
available in a [separate GitHub
project](https://github.com/patrickbwarren/python3-HNC-solver).  It
might also be helpful to browse the documentation for this project
too.

#### Citation

You are of course expected to use and modify this code for your own
projects!  If you end up publishing results based on this code, or a
derived version, please cite :

*Screening properties of Gaussian electrolyte models, with application
to dissipative particle dynamics,* P. B. Warren, A.  Vlasov, L. Anton
and A. J. Masters, [J. Chem. Phys. **138**, 204907
(2013)](http://jcp.aip.org/resource/1/jcpsa6/v138/i20/p204907_s1 "AIP
link").

If you have difficulty getting hold of this article, contact me and I
will send a reprint. Alternatively, a preprint version can be found at
[arXiv:1303.0891](http://front.math.ucdavis.edu/1303.0891 "arXiv link").

Here is a BibTeX entry for the same :

```
@article{WVA+13,
  author         = {Warren, P. B. and Vlasov, A. and Anton, L. and Masters, A. J.},
  title          = {Screening properties of {G}aussian electrolyte models,
                    with application to dissipative particle dynamics},
  journal        = {J. Chem. Phys.},
  volume         = {138},
  pages          = {204907},
  year           = 2013
}
```

#### Installation Notes

The main code is written as a FORTRAN 90 module, and can be compiled
using `gfortran` from the [GNU compiler
collection](https://gcc.gnu.org/ "GNU website").  A basic Makefile is
provided.  The code is designed to be run from python using the `f2py`
interface from [SciPy](http://www.scipy.org/ "SciPy website"), though
some examples of 'vanilla' FORTRAN 90 driver codes are included.  The
default target in the `Makefile` builds the required files for the
`f2py` interface.

The [LAPACK](http://www.netlib.org/lapack/ "LAPACK webpage")
linear algebra library and [FFTW](http://www.fftw.org/ "FFTW website")
fast Fourier transform library are both required.

The compiler tools including `f2py`, and pre-built versions of the
`lapack` and `fftw` libraries, are available for most modern GNU/Linux
distributions.

#### Roadmap

Whilst the code is very functional there still remain legacy issues
from its origins in monolithic FORTRAN 90:

* aspects of the initialisation (potential model definitions) and
  post-processing (thermodynamic calculations) do not necessarily have
  to be done in FORTRAN;

* many parameters are hard-wired into the code in the
  potential functions, with in some cases overloaded definitions;

* adding new potentials or changing / verifying the current ones is
  clunky and requires access to the underlying FORTRAN code;

* the solver grid is allocated at startup, and cannot currently be
  re-sized, nor can multiple grids be used in the same application;

* certain features of the current code, such as reporting the solver
  status, involve string manipulations which are certainly better
  handled directly in python;

* overall, the current interface doesn't take full advantage of what
  might be considered modern python functionality and design patterns.

A longer term aim therefore is to split out the central
Ornstein-Zernike solver with the closure approximations into a
streamlined standalone suite of FORTRAN 90 kernel functions.  This
would particularly preserve the technically complex and bespoke Ng
accelerated convergence scheme.  The initialisation of specific
potential models (DPD, RPM, etc) would then be implemented using
[NumPy](https://numpy.org/) in a separate python module.  This
would make the potential model definitions and the associated
parameter specifications more accessible, user-friendly, and better in
accord with pythonic idioms.

Initially, the downstream thermodynamic calculations could be left in
FORTRAN but in time these too could be moved into python, using for
example the `trapz` function in NumPy to handle the integrations, along
the lines of the related pure python Ornstein-Zernike-HNC solver
[here](https://github.com/patrickbwarren/python3-HNC-solver).

Thus, functionally, the code would split into two:

* a `sunlighthnc-core` module, wrapping the core FORTRAN 90 functions;
* a `sunlighthnc` module, providing the application user interface.

If and when this comes to pass, this will be version 2.0 of the code !!

#### Copying

SunlightHNC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SunlightHNC is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see
<http://www.gnu.org/licenses/>.

#### Copyright

SunlightDPD is based on an original code copyright &copy; 2007 Lucian
Anton, with modifications copyright &copy; 2008, 2009 Andrey Vlasov, and
additional modifications copyright &copy; 2009-2018 Unilever UK Central
Resources Ltd (Registered in London number 29140;
Registered Office: Unilever House, 100 Victoria Embankment, London EC4Y 0DY, UK).
Later modifications copyright &copy; 2020-2025 Patrick B Warren (STFC).

#### Contact

Send email to patrick.warren{at}stfc.ac.uk
