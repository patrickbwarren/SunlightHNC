#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This file includes unicode characters

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Copyright (c) 2026 Patrick B Warren <patrick.warren@stfc.ac.uk>.

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

# Reproduce Fig 2 in Warren et al., J. Chem. Phys. 138, 204907 (2013).
# Monte-Carlo data from that figure is in urpm_eos.ods.

# This figure plots the excess pressure and energy for the ultrasoft
# restricted primitive model (URPM) as a function of density for two
# values of the Bjerrum length lB.

import oz_aux
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from oz import wizard
from numpy import pi as π

parser = argparse.ArgumentParser(description='figure 2 in JCP 138, 204907 (2013)')
parser.add_argument('--dpi', default=72, type=int, help='resolution (dpi) for image output, default (for pdf) 72')
parser.add_argument('-o', '--output', help='output figure to, eg, pdf file')
args = parser.parse_args()

MC = pd.read_excel('urpm_eos.ods')
MC['pex'] = MC.p - MC.rho

lB_vals = MC.lB.unique()

half = 0.5 * np.ones(2)
grid = oz_aux.Grid(wizard, ncomp=2, ng=4096, deltar=0.01)

data = []
for lB in lB_vals:
    urpm = oz_aux.ultrasoft_restricted_primitive_model(grid, lB)
    for ρz in np.geomspace(1e-3, 10, 21):
        soln = oz_aux.solve(urpm, half*ρz, 'HNC')
        data.append([x+0 for x in [lB, ρz, soln.pex, soln.uex, soln.error]])
HNC = pd.DataFrame(data, columns=['lB', 'rho', 'pex', 'e', 'conv'])

lw, ms = 1.2, 4
gen_lw, line_lw = 1.2, 1.2
tick_fs, label_fs, legend_fs = 14, 14, 14

fig, ax = plt.subplots(figsize=(6, 4), dpi=args.dpi)

colors = [f'tab:{c}' for c in ['blue', 'red']]

ρ = np.array([1e-3, 10])

for lB, color in zip(lB_vals, colors):
    κ = np.sqrt(4*π*lB*ρ)
    ax.loglog(ρ, κ**3/(24*π), ':', lw=lw, c=color)
    cut = HNC.lB == lB
    ax.loglog(HNC[cut].rho, -HNC[cut].pex, '-', lw=lw, c=color, label=f'$l_B={lB}$')
    ax.loglog(HNC[cut].rho, -HNC[cut].e/3, '-', lw=lw, c=color)
    cut = MC.lB == lB
    ax.loglog(MC[cut].rho, -MC[cut].pex, 'o', ms=ms, c=color)
    ax.loglog(MC[cut].rho, -MC[cut].e/3, 's', ms=ms, c=color)

ax.legend(loc='upper left', frameon=False, fontsize=legend_fs, labelspacing=0.5)

ax.set_xlim(5e-3, 2)
xticks = [0.01, 0.1, 1]
xlabels = ['0.01', '0.1', '1']
ax.set_xticks(xticks, labels=xlabels)
ax.set_xlabel(r'$\rho_z$', fontsize=label_fs)

ax.set_ylim(1e-4, 10)
yticks = [1e-4, 1e-3, 0.01, 0.1, 1, 10]
ylabels = ['$10^{-4}$', '$10^{-3}$', '0.01', '0.1', '1', '10']
ax.set_yticks(yticks, labels=ylabels)
yaxlabels = [r'$-p^{\mathrm{ex}}$', r'$-\langle U\rangle/3V$']
ax.set_ylabel(' ; '.join(yaxlabels), fontsize=label_fs)

ax.minorticks_off()
ax.tick_params(direction='in', width=gen_lw, length=5, top=True, right=True, labelsize=tick_fs)
for spine in ax.spines:
    ax.spines[spine].set_linewidth(gen_lw)

plt.tight_layout()

if args.output:
    plt.savefig(args.output, bbox_inches='tight', pad_inches=0.05)
    print('Figure saved to', args.output)
else:
    plt.show()
