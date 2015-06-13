#!/usr/bin/python3

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  
# Additional modifications copyright (c) 2009-2015 Unilever UK Central
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

# A plot of the chemical potential difference against x reveals nearly
# linear behaviour, which can be confirmed in Monte-Carlo simulations
# (unreported).  This linear dependence can be represented by Delta F
# = chi (1 - 2x) where chi is the Flory chi-parameter and is the basic
# reason why DPD is so good at representing fluid mixtures which fit
# Flory-Huggins (regular solution) theory.

from oz import wizard as w

w.ncomp = 2
w.initialise()

rho = 3.0
A = 25.0
dA = 5.0

w.arep[0,0] = w.arep[1,1] = A
w.arep[0,1] = A + dA
w.dpd_potential(1)

n = 41

for i in range(n):

    x = i / (n - 1.0)

    w.rho[0] = (1.0 - x) * rho
    w.rho[1] = x * rho
    w.hnc_solve()

    print("%f\t%f\t%g" % (x, w.muex[1]-w.muex[0], w.error))

