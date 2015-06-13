# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
# modifications copyright (c) 2009-2013 Unilever UK Central Resources
# Ltd (Registered in England & Wales, Company No 29140; Registered
# Office: Unilever House, Blackfriars, London, EC4P 4BQ, UK).

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
FFLAGS = -Wall -fbounds-check
LIBFLAGS =  -L/usr/lib -llapack -lfftw3
INCFLAGS = -I/usr/include/

default: oz.so

all : driver fftw_test

oz.so : oz_mod.f90
	f2py3 -c $< -m oz $(INCFLAGS) $(LIBFLAGS)

driver : driver.o oz_mod.o
	$(F90) -o $@ $^ $(LIBFLAGS)

fftw_test : fftw_test.o
	$(F90) -o $@ $^ $(LIBFLAGS)

%.o : %.f90
	$(F90) -c $(FFLAGS) $(INCFLAGS) $<

driver.o : driver.f90 oz_mod.o

clean:
	rm -f *~
	rm -f *.mod *.o *.so
	rm -f driver fftw_test
