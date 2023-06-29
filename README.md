## General purpose HNC code

The current release is version 1.13.

Version 1.13 allows dissimilar and non-additive sphere diameters in the RPM and hard sphere potentials.  
Version 1.12 adds direct calculation of HNC free energy (see documentation).  
Version 1.11 fixes another bug: missing mean-field contribution to thermodynamics for off-SYM RPM potential.  
Version 1.10 fixes a bug: missing contact virial pressure contribution in the HNC free energy.  
Version 1.9 brings the source code up to modern FORTRAN standards and corrects a problem with EXP in v1.8.  
Version 1.8 includes the MSA closure and various potentials for Slater smeared charges.  
Version 1.7 handles an arbitrary number of components (rather than a maximum of three, as for previous versions).

#### Features

* FORTRAN 90 based, with example python driver scripts;
* HNC, MSA, RPA, and EXP closures;
* fast Ng solver, and Ng decomposition for long range potentials;
* multicomponent (arbitrary number of components);
* hard core (RPM-like) and soft core (DPD-like) potentials;
* full structural thermodynamics;
* fully open source.

See [PDF documentation](oz_doc.pdf "oz_doc.pdf") for details.

#### Citation

You are of course expected to use and modify this code for your own
projects!  
If you end up publishing results based on this code, or a
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

#### Contact

Send email to patrick.warren{at}stfc.ac.uk
