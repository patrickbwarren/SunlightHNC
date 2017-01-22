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

# The results of this calculation can be directly compared with Fig 4
# of the Groot and Warren [J. Chem. Phys. v107, 4423 (1997)].  The
# data from that figure is coded below.

#       rho   (p-rho)/(A*rho^2)

data = [[0.0, 0.0379935086163],
        [1.5, 0.0751786298043],
        [2.5, 0.0886823425022],
        [3.0, 0.0924251622846],
        [3.5, 0.0946639891655],
        [4.0, 0.0965259421847],
        [5.0, 0.0987451548125],
        [6.0, 0.0998358473824],
        [7.0, 0.1005510671090],
        [8.0,  0.102017933031]]

xdata = list(data[i][0] for i in range(len(data)))
ydata = list(data[i][1] for i in range(len(data)))

from oz import wizard as w

w.initialise()
w.arep[0,0] = A = 25.0
w.dpd_potential()

npt = 41
rhomax = 10.0

x = []
y = []

for i in range(npt):
    w.rho[0] = rho = rhomax * (i + 1.0) / npt
    w.hnc_solve()
    x.append(rho)
    y.append((w.press-rho)/(A*rho*rho))
    print("%f\t%g\t%g" % (rho, (w.press-rho)/(A*rho*rho), w.error))

import matplotlib.pyplot as plt

plt.plot(xdata, ydata, 'ro', label='Groot & Warren (1997)')
plt.plot(x, y, label='HNC')
plt.xlabel('$\\rho$')
plt.ylabel('$(p-\\rho)/A\\rho^2$')
plt.legend(loc='lower right')

plt.show()
