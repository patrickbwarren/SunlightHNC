#!/usr/bin/python3

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  
# Additional modifications copyright (c) 2009-2016 Unilever UK Central
# Resources Ltd (Registered in England & Wales, Company No 29140;
# Registered Office: Unilever House, Blackfriars, London, EC4P 4BQ,
# UK).

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

# The results can be directly compared with Fig 5 of Wijmans et al
# [J. Chem. Phys. v114, 7644 (2001)].  The chemical potential
# difference muex_1 - muex_0 is exactly equal to the free energy
# change to convert a particle of species 0 into a particle of species
# 1.  Note that the calculation is done as a two component problem
# with a mole fraction of the second species set to zero.

from oz import wizard as w

w.ncomp = 2
w.initialise()

rho = 3.0
A = 25.0
Amax = 50.0

w.arep[0,0] = w.arep[1,1] = A
w.rho[0] = rho
w.rho[1] = 0.0

npt = 41

for i in range(npt):

    dA = (Amax - A) * i / (npt - 1.0)

    w.arep[0,1] = A + dA
    w.dpd_potential(1)
    w.hnc_solve()

    print("%f\t%f\t%g" % (dA, w.muex[1]-w.muex[0], w.error))

