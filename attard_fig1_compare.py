#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This file includes unicode characters

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Copyright (c) 2009-2019 Unilever UK Central Resources Ltd
# (Registered in England & Wales, Company No 29140; Registered
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

# Reproduce Fig 1a in Attard, Phys. Rev. E 48, 3604-21 (1993).

# This figure plots the total correlation functions for the RPM model
# of an electrolyte at 0.5 M, 2 M, and 5 M, with ion diameter 0.5 nm,
# in water with ε_r = 78.5 ant T = 300 K.

import argparse
import numpy as np
import math as m
import matplotlib.pyplot as plt

from oz import wizard as w
from numpy import pi as π

w.ncomp = 2
w.ng = 2**14 # or 2**16
w.deltar = 1e-3 # or 5e-4

w.initialise()

# fundamental quantities in SI units 

e = 1.6e-19 # charge on electron
ε0 = 8.854e-12 # permittivity of vacuum
kB = 1.38e-23 # Boltzmann constant
NA = 6.02e23 # Avogadro constant

Å = 1e-10 # in m, used for conversions

# solvent (water) parameters

T = 300 # temperature
εr = 78.5 # relative permittivity

# Bjerrum length (in m)

lb = e**2 / (4*π*εr*ε0*kB*T)

# solute (electrolyte) parameters

d = 5 * Å # ion diameter
concs = [0.5, 2, 5] # Molar units

print('Bjerrum length = %g Å = %g d' % (lb/Å, lb/d))

# We use the diameter d as a base length unit (assuming w.sigma = 1)

w.lb = lb / d

w.rpm_potential()

# plt.xkcd()

fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.set_xlim([5, 20]) # horizontal axis will be in Å
ax1.set_xticks([5, 10, 15, 20])
ax1.set_xlabel('r (Å)')
ax1.set_ylim([-0.5, 0.5])
ax1.set_yticks(np.around((np.arange(0, 11)*0.1-0.5), decimals=1).tolist())
ax1.set_title('h(r)')

ax2.set_xlim([5, 50]) # horizontal axis will be in Å
ax2.set_xticks(list(range(5, 55, 5)))
ax2.set_xlabel('r (Å)')
ax2.set_ylim([-8, 1])
ax2.set_yticks(list(range(-8, 2, 1)))
ax2.set_title('log10 h(r)')

imin = int(1.0 / w.deltar)

styles = ['-', '--', ':']

for i, (conc, style) in enumerate(zip(concs, styles)):

    # Conversion molar concentration to molecular density
    # ρ = (c/M) * 10^3 * N_A where the 10^3 comes
    # from the number of litres in 1 m^3.

    ρd3 = conc * 1e3 * NA * d**3

    w.rho[0] = ρd3
    w.rho[1] = ρd3

    w.hnc_solve()

    if w.return_code: exit()

    s = str(w.closure_name, 'utf-8').strip()
    print('conc = %5.1f M \tρ*d^3 = %8.5f \tϕ = %8.5f \t%s error = %g' %
          (conc, ρd3, π*ρd3/6, s, w.error))

    imax = int(4.0 / w.deltar)
    
    ax1.plot(w.r[imin:imax] * d / Å, w.hr[imin:imax, 0, 0], 'k'+style)
    ax1.plot(w.r[imin:imax] * d / Å, w.hr[imin:imax, 0, 1], 'k'+style)

    # In Attard's figure the line for 0.5 M is surely h_{+-}

    imax = int(10.0 / w.deltar)

    if i == 0:
        h = w.hr[imin:imax, 0, 1]
    else:
        h = w.hr[imin:imax, 0, 0]

    h[h < 0] = 1e-20 # chop off negative regions

    ax2.plot(w.r[imin:imax] * d / Å, np.log10(h), 'k'+style)

    # The commented out variant plots r h for both ++ and +- functions

    # for j, color in zip([0, 1], ['r', 'b']):
    #     h = w.r[imin:imax] * w.hr[imin:imax, 0, j]
    #     h[h < 0] = 1e-20
    #     plt.plot(w.r[imin:imax] * d / Å, np.log10(h), color+style)

plt.show()
