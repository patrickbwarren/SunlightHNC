## General purpose HNC code

Currently at version 1.6.  Previous versions can be found at
http://sunlightdpd.sourceforge.net/.

#### Features

* FORTRAN 90 based, with example python driver scripts ;
* fast Ng solver, and Ng decomposition for electrostatics ;
* multicomponent (up to three components) ;
* hard and soft core (DPD) potentials ;
* full structural thermodynamics ;
* fully open source.

See [PDF documentation](oz_doc.pdf "oz_doc.pdf") for details.

#### Citation

You are of course welcome to use and modify this code for your own
projects. If you end up publishing results based on this code, or a
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

The main code is written as a FORTRAN 90 module, and can be compiled using
`gfortran` from the [GNU compiler collection](https://gcc.gnu.org/
"GNU website").  It is designed to work with `f2py` from
[SciPy](http://www.scipy.org/ "SciPy website").  The default target in
the `Makefile` builds the required files for this.

The linear algebra library
[LAPACK](http://www.netlib.org/lapack/ "LAPACK webpage") and fast
Fourier transform library [FFTW](http://www.fftw.org/ "FFTW website")
should be installed.

The compiler tools including `f2py`, and pre-built versions of the
libraries, are available for most modern GNU/Linux distributions.

#### Copying

SunlightHNC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
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
additional modifications copyright &copy; 2009-2015 Unilever UK Central
Resources Ltd (Registered in England & Wales, Company No 29140;
Registered Office: Unilever House, Blackfriars, London, EC4P 4BQ, UK).

#### Contact

Send email to patrick{dot}warren{at}unilever{dot}com

Send paper mail to Dr Patrick B Warren, Unilever R&D Port Sunlight,
Quarry Road East, Bebington, Wirral, CH63 3JW, UK.
