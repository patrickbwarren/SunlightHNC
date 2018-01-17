#!/usr/bin/env python3

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  
# Additional modifications copyright (c) 2009-2017 Unilever UK Central
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

# Compute log_10 activity in infinite dilution limit for a binary DPD
# fluid and compare to 0.144 * Delta A

from math import log
from oz import wizard as w

w.ncomp = 2 # two components
w.initialise()

w.rho[0] = 3.0
w.rho[1] = 0.0 # zero density (still works!)

# Warm up HNC by ramping the repulsion amplitude

Amin = 25.0
Amax = 106.5
npt = 50

for i in range(npt):
    A = Amin + (Amax-Amin)*i/(npt-1.0)
    w.arep[0,0] = w.arep[0,1] = w.arep[1,1] = A
    w.dpd_potential(1)
    w.hnc_solve()
    print("%f\t%f\t%f\t%g" % (A, w.muex[0], w.muex[1], w.error))

# Now ramp the extra repulsion amplitude

dAmin = -5.0
dAmax = 20.0
npt = 11

x = [0.0 for i in range(npt)]
y = [0.0 for i in range(npt)]

for i in range(npt):
    x[i] = dA = dAmin + (dAmax-dAmin)*i/(npt-1.0)
    w.arep[0,1] = A + dA
    w.dpd_potential(1)
    w.hnc_solve()
    y[i] = (w.muex[1] - w.muex[0]) / log(10.0) # ie log_10(gamma)
    print("%f\t%f\t%f\t%g" % (dA, w.muex[0], w.muex[1], w.error))

xx = [x[0], x[-1]] # first and last elements
yy = [0.144*x for x in xx] # straight line slope 0.144

import matplotlib.pyplot as plt

fig, ax = plt.subplots()

plt.plot(x, y, label='HNC')
plt.plot(xx, yy, '--', label='0.144 $\\Delta a$')
plt.xlabel('$\\Delta a$')
plt.ylabel('$log_{10}\\gamma$')
plt.legend(loc='lower right')

ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')

plt.show()
