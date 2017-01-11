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

# The results can be directly compared with Fig 10.2 of Hansen and
# McDonald.  This figure is based on Rasaiah et al, J. Chem. Phys. 56,
# 248 (1972).  Here is data from part of Table II from that paper :

#         c_s/M      -U/NkT (MC)      contact (MC)     p/(rho kT)

data = [[0.00911,  0.1029, 0.0013,  0.0044, 0.0007,  0.9701, 0.0008],
        [0.10376,  0.2739, 0.0014,  0.0359, 0.0011,  0.9445, 0.0012],
        [0.42502,  0.4341, 0.0017,  0.1217, 0.0045,  0.9774, 0.0046],
        [1.0001,   0.5516, 0.0016,  0.2777, 0.0045,  1.094,  0.005],
        [1.9676,   0.6511, 0.0020,  0.5625, 0.0088,  1.346,  0.009]]

# The MC results are reported as a value +/- error

# The pressure satisfies p/(rho kT) = 1 + U/(3NkT) + contact

# Parameters are sigma = 0.425nm, lB = 0.71nm, hence lB/sigma = 1.68.
# Molar concentrations can be converted to ion concentrations by rho =
# 2 * c_s/M * 10^3 * N_A = 1.204 * c_s/M nm^(-3).  The factor two
# arises because each molecule of salt contributes two ions (an anion
# and a cation).  Thus rho d^3 ranges from about 5e-4 to 0.2 in the
# above table.

# Calculations are done for both HNC and MSA, but only HNC printed

import math as m

xdat = list(m.sqrt(0.425**3 * 1.204 * data[i][0]) for i in range(len(data)))
nundat = list(data[i][1] for i in range(len(data)))
ctcdat = list(data[i][3] for i in range(len(data)))
compdat = list(data[i][5] for i in range(len(data)))

from oz import wizard as w

w.ng = 4096
w.ncomp = 2
w.deltar = 0.01

w.initialise()

w.lb = 0.71 / 0.425
w.sigma = 1.0
w.kappa = -1.0

w.soft_rpm_potential(0)

w.write_params();

npt = 41

rholo = 0.0005 
rhohi = 0.2
lrlo = m.log(rholo)
lrhi = m.log(rhohi)

x = []
y1 = []
y2 = []
y3 = []
y4 = []

for i in range(npt):
    lr = lrlo + (lrhi - lrlo) * i / (npt - 1.0)
    rho = m.exp(lr)
    x.append(m.sqrt(rho))
    w.rho[0] = rho / 2.0
    w.rho[1] = rho / 2.0
    w.hnc_solve()
    comp = 1.0 + w.un / 3.0 + w.cf_gc
    y1.append(-w.un)
    y2.append(comp)
    print("%f\t%f\t%f\t%g" % (m.sqrt(rho), w.un, comp, w.error))
    w.msa_solve()
    comp = 1.0 + w.un / 3.0 + w.cf_gc
    y3.append(-w.un)
    y4.append(comp)

import matplotlib.pyplot as plt

plt.plot(xdat, nundat, 'ro', label='Rasaiah et al (1972): -U/V')
plt.plot(xdat, compdat, 'co', label='Rasaiah et al (1972): p')
plt.plot(x, y1, 'b-', label='HNC')
plt.plot(x, y2, 'b-')
plt.plot(x, y3, 'b--', label='MSA')
plt.plot(x, y4, 'b--')
plt.xlabel('$(\\rho\\sigma^3)^{1/2}$')
plt.ylabel('$-U_n$ and $p/\\rho k_\\mathrm{B}T$')
plt.legend(loc='upper left')

plt.show()
