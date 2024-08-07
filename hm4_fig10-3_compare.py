#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This file includes unicode characters

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Copyright (c) 2009-2019 Unilever UK Central Resources Ltd
# (Registered in England & Wales, Company No 29140; Registered Office:
# Unilever House, Blackfriars, London, EC4P 4BQ, UK).  Additional
# modifications copyright (c) 2020-2024 Patrick B Warren
# <patrick.warren@stfc.ac.uk> and STFC.

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

# The results can be directly compared with Fig 10.3 of Hansen and
# McDonald (4th edition), based on Rasaiah et al., J. Chem. Phys. 56,
# 248 (1972).  Here is data from part of Table II from that paper :

#         c_s/M      -U/NkT (MC)      contact (MC)     p/(rho kT)

data = [[0.00911,  0.1029, 0.0013,  0.0044, 0.0007,  0.9701, 0.0008],
        [0.10376,  0.2739, 0.0014,  0.0359, 0.0011,  0.9445, 0.0012],
        [0.42502,  0.4341, 0.0017,  0.1217, 0.0045,  0.9774, 0.0046],
        [1.0001,   0.5516, 0.0016,  0.2777, 0.0045,  1.094,  0.005],
        [1.9676,   0.6511, 0.0020,  0.5625, 0.0088,  1.346,  0.009]]

# The MC results are reported as a value +/- error

# The pressure satisfies p/(rho kT) = 1 + U/(3NkT) + contact

# Parameters are d = 0.425 nm, lB = 0.71 nm, hence lB / d = 1.68.
# Molar concentrations can be converted to ion concentrations by ρ =
# 2 * c_s/M * 10^3 * N_A = 1.204 * c_s/M nm^(-3).  The factor two
# arises because each molecule of salt contributes two ions (an anion
# and a cation).  Thus ρd³ ranges from about 5e-4 to 0.2 in the
# above table.

# Calculations are done for both HNC, MSA and EXP.  Note that EXP
# requires the hard sphere reference state to be computed.  This in
# turn appears to require a cold start each time since the EXP
# solution for the RPM appears to deviate too far from the hard sphere
# solution to be a good starting point.

import oz_aux
import argparse
import numpy as np
from numpy import sqrt, exp, log10
from numpy import log as ln
from oz import wizard

d = 0.425
lb = 0.71

ndata = len(data)
sqrdat = np.array([sqrt(d**3 * 1.204 * data[i][0]) for i in range(ndata)])
nundat = np.array([data[i][1] for i in range(ndata)])
ctcdat = np.array([data[i][3] for i in range(ndata)])
cmpdat = np.array([data[i][5] for i in range(ndata)])

parser = argparse.ArgumentParser(description='Reproduce Fig 10-4 from Hansen and McDonald (4th edn)')
parser.add_argument('--npt', default=41, type=int, help='number of points, default 41')
parser.add_argument('--hnc', action='store_true', help='include HNC prediction')
parser.add_argument('--msa', action='store_true', help='include MSA prediction')
parser.add_argument('--exp', action='store_true', help='include EXP prediction')
parser.add_argument('--all', action='store_true', help='include all predictions')
args = parser.parse_args()

if args.all:
    args.hnc = args.msa = args.exp = True    
elif not args.hnc and not args.msa and not args.exp:
    args.hnc = args.msa = True

grid = oz_aux.Grid(wizard, ncomp=2, ng=8192, deltar=0.01)

model = oz_aux.restricted_primitive_model(grid, lb/d)

x = sqrt(np.logspace(log10(5e-4), log10(0.2), args.npt))
y = np.empty((args.npt, 6))

u = np.array([0.5, 0.5])

if args.hnc:
    for i in range(args.npt):
        rho = x[i]**2
        soln = oz_aux.solve(model, rho*u, 'HNC')
        y[i, 0] = -soln.uex/rho
        y[i, 1] = soln.press/rho
        print("%i\t%f\t%f\t%f\t%g\t%g\tHNC" %
              (i, sqrt(rho), soln.uex/rho, soln.wizard.cf_gc, soln.wizard.cf_xc, soln.error))
        
if args.msa:
    for i in range(args.npt):
        rho = x[i]**2
        soln = oz_aux.solve(model, rho*u, 'MSA')
        y[i, 2] = -soln.uex/rho
        y[i, 3] = soln.press/rho
        print("%i\t%f\t%f\t%f\t%g\tMSA" %
              (i, sqrt(rho), soln.uex/rho, soln.wizard.cf_gc, soln.error))

if args.exp: # temporary fix here as details not yet implemented in oz_aux
    for i in range(args.npt):
        rho = x[i]**2
        wizard.lb = lb / d
        wizard.sigma = 1.0
        wizard.kappa = -1.0
        wizard.rho = rho * u
        wizard.cold_start = True
        wizard.hs_potential()
        wizard.msa_solve()
        wizard.save_reference()
        wizard.rpm_potential()
        wizard.msa_solve()
        wizard.exp_refine()
        y[i, 4] = -wizard.uex/rho
        y[i, 5] = wizard.press/rho
        print("%i\t%f\t%f\t%f\t%g\tEXP" %
              (i, sqrt(rho), wizard.uex/rho, wizard.cf_gc, wizard.error))

import matplotlib.pyplot as plt

plt.plot(sqrdat, nundat, 'ro', label=r"Rasaiah et al. '72: Energy")
plt.plot(sqrdat, cmpdat, 'co', label=r"Rasaiah et al. '72: Pressure")
if args.hnc:
    plt.plot(x, y[:, 0], 'k-', label='HNC')
    plt.plot(x, y[:, 1], 'k-')
if args.msa:
    plt.plot(x, y[:, 2], 'b--', label='MSA (virial)')
    plt.plot(x, y[:, 3], 'b--')
if args.exp:
    plt.plot(x, y[:, 4], 'r--', label='EXP')
    plt.plot(x, y[:, 5], 'r--')
plt.xlabel(r'$(\rho\sigma^3)^{1/2}$')
plt.ylabel(r'$-\beta U^{ex}/N$ and $\beta P/\rho$')
plt.legend(loc='upper left')

plt.show()
