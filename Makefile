# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
# modifications copyright (c) 2009-2017 Unilever UK Central Resources
# Ltd (Registered in England & Wales, Company No 29140; Registered
# Office: Unilever House, Blackfriars, London, EC4P 4BQ, UK).
# Additional modifications copyright (c) 2020-2021 Patrick B Warren
# <patrick.warren@stfc.ac.uk> and STFC.

# SunlightDPD is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# SunlightDPD is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with SunlightDPD.  If not, see <http://www.gnu.org/licenses/>.

F90 = gfortran
# FFLAGS = -Wall -fbounds-check
# FFLAGS = -g -Wall -fcheck=all
FFLAGS = -O2
LIBFLAGS =  -L/usr/lib -llapack -lfftw3
INCFLAGS = -I/usr/include/
PDFL = pdflatex -synctex=1

DRIVERS = driver1 driver2 driver3 driver4

default: oz.so

all : $(DRIVERS) fftw_test oz.so

docs : oz_doc.pdf

oz_doc.pdf : oz_doc.tex gofr.png sofk.png
	$(PDFL) oz_doc.tex
	$(PDFL) oz_doc.tex
	$(PDFL) oz_doc.tex

drivers : $(DRIVERS)

oz.so : oz_mod.f90
	f2py3 --overwrite-signature $< -m oz -h oz.pyf
	f2py3 -c $< oz.pyf $(INCFLAGS) $(LIBFLAGS)

driver1 : driver1.o oz_mod.o
	$(F90) -o $@ $^ $(LIBFLAGS)

driver2 : driver2.o oz_mod.o
	$(F90) -o $@ $^ $(LIBFLAGS)

driver3 : driver3.o oz_mod.o
	$(F90) -o $@ $^ $(LIBFLAGS)

driver4 : driver4.o oz_mod.o
	$(F90) -o $@ $^ $(LIBFLAGS)

fftw_test : fftw_test.o
	$(F90) -o $@ $^ $(LIBFLAGS)

%.o : %.f90
	$(F90) -c $(FFLAGS) $(INCFLAGS) $<

driver1.o : driver1.f90 oz_mod.o
driver2.o : driver2.f90 oz_mod.o
driver3.o : driver3.f90 oz_mod.o
driver4.o : driver4.f90 oz_mod.o

clean:
	rm -f *~
	rm -f *.pyf *.mod *.o *.so
	rm -f oz_doc.aux oz_doc.toc oz_doc.log
	rm -f $(DRIVERS) fftw_test
