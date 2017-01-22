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

# The results for the compressibility factor beta P / rho for hard
# spheres as a function of packing fraction eta can be directly
# compared with Fig 4.2 of Hansen and McDonald.  Calculations are done
# for both HNC and MSA.  For hard spheres the MSA is equivalent to
# Percus-Yevick which is solvable, so a comparison is also given with
# the exact results in that case.

from oz import wizard as w

# High resolution calculations are required for this to be accurate!

w.ncomp = 1
w.ng = 65536
w.deltar = 0.001

w.initialise()

w.sigma = 1.0

w.hs_potential()

# Calculate the compressibility-route pressure by integrating the
# compressibility along an isotherm, as in examples.py.  Note that you
# also need a pretty fine-grained integration step for this to work.

npt = 150
rhohi = 0.43 * 6.0 / (w.pi * w.sigma**3)
drho = rhohi / npt

x = []
y = []
yy = []
yyy = []
y1 = []
y2 = []
y3 = []
y4 = []

hnc_p_xc = hnc_prev = 0.0
msa_p_xc = msa_prev = 0.0

for i in range(npt):
    rho = drho * (i + 1.0)
    eta = w.pi * rho * w.sigma**3 / 6.0
    carnahan_starling = (1 + eta + eta**2 - eta**3) / (1 - eta)**3
    percus_yevick_virial = (1 + 2*eta + 3*eta**2) / (1 - eta)**2
    percus_yevick_comp = (1 + eta + eta**2) / (1 - eta)**3
    x.append(eta)
    y.append(carnahan_starling)
    yy.append(percus_yevick_virial)
    yyy.append(percus_yevick_comp)
    w.rho[0] = rho
    w.hnc_solve()
    hnc_p_xc = hnc_p_xc + 0.5*drho*(hnc_prev + w.comp_xc)
    hnc_prev = w.comp_xc
    y1.append(w.press/rho)
    y2.append(1.0 + hnc_p_xc/rho)
    print("%f\t%f\t%f\t%g\tHNC" % (eta, w.press/rho, 1.0+hnc_p_xc/rho, w.error))
    w.msa_solve()
    msa_p_xc = msa_p_xc + 0.5*drho*(msa_prev + w.comp_xc)
    msa_prev = w.comp_xc
    y3.append(w.press/rho)
    y4.append(1.0 + msa_p_xc/rho)
    print("%f\t%f\t%f\t%g\tMSA" % (eta, w.press/rho, 1.0+msa_p_xc/rho, w.error))

import matplotlib.pyplot as plt

plt.plot(x, y, 'k-', label='Carnahan-Starling')
plt.plot(x, yy, 'k--', label='PY(v) exact')
plt.plot(x, yyy, 'k-.', label='PY(c) exact')
plt.plot(x, y1, 'b-', label='HNC(v)')
plt.plot(x, y2, 'b--', label='HNC(c)')
plt.plot(x, y3, 'r-', label='MSA(v)')
plt.plot(x, y4, 'r--', label='MSA(c)')

plt.xlabel('$\\eta$')
plt.ylabel('$\\beta P/\\rho$')

plt.legend(loc='upper left')

plt.show()
