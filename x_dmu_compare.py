#!/usr/bin/env python3

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

# A plot of the chemical potential difference against x reveals nearly
# linear behaviour, which can be confirmed in Monte-Carlo simulations
# (unreported).  This linear dependence can be represented by Delta F
# = chi (1 - 2x) where chi is the Flory chi-parameter and is the basic
# reason why DPD is so good at representing fluid mixtures which fit
# Flory-Huggins (regular solution) theory.

# The following table contains unpublished Monte-Carlo calculations of
# the chemical potential difference, using trial particle identity
# swaps.

#        x     d(mu)  std-error

data = [[0.0,  1.702,  0.004],
        [0.1,  1.281,  0.003],
        [0.2,  0.915,  0.003],
        [0.3,  0.597,  0.004],
        [0.4,  0.282,  0.004],
        [0.5, -0.002,  0.003],
        [0.6, -0.290,  0.003],
        [0.7, -0.596,  0.003],
        [0.8, -0.920,  0.004],
        [0.9, -1.284,  0.003],
        [1.0, -1.702,  0.004]]

xdata = list(data[i][0] for i in range(len(data)))
ydata = list(data[i][1] for i in range(len(data)))
edata = list(data[i][2] for i in range(len(data)))

from oz import wizard as w

w.ncomp = 2
w.initialise()

rho = 3.0
A = 25.0
dA = 5.0

w.arep[0,0] = w.arep[1,1] = A
w.arep[0,1] = A + dA
w.dpd_potential()

n = 41

xarr = []
yarr = []

for i in range(n):
    x = i / (n - 1.0)
    w.rho[0] = (1.0 - x) * rho
    w.rho[1] = x * rho
    w.hnc_solve()
    xarr.append(x)
    yarr.append(w.muex[1]-w.muex[0])
    print("%f\t%f\t%g" % (x, w.muex[1]-w.muex[0], w.error))

import matplotlib.pyplot as plt

# plt.errorbar(xdata, ydata, yerr=edata, fmt='ro', label='Monte-Carlo')
# error bars are smaller than symbols so just plot data

plt.plot(xdata, ydata, 'ro', label='Monte-Carlo')
plt.plot(xarr, yarr, label='HNC')
plt.xlabel('$x$')
plt.ylabel('$\\Delta\\mu$')
plt.legend(loc='upper right')

plt.show()
