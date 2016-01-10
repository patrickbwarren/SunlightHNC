#!/usr/bin/python3

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Copyright (c) 2009-2016 Unilever UK Central Resources Ltd
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

# Scan the screening behaviour of various URPM models, using command
# line parameters as a control.

# Typical use starts off with a broad (--lo, --hi) interval where the
# numbers are -log10(rho_z).  Then the interval which contains
# the Kirkwood crossover between pure exponential and oscillatory
# decay can be narrowed by a factor 10 each time.  For example :
#   python urpm_scan.py --lo=1.0 --hi=2.0
#   python urpm_scan.py --lo=1.5 --hi=1.6
#   python urpm_scan.py --lo=1.52 --hi=1.53
# Inspecting the output of the last narrows the interval to (1.524,
# 1.526) - this for the default parameters.

import argparse
import math as m
import matplotlib.pyplot as plt
from oz import wizard as w

parser = argparse.ArgumentParser(description='URPM h(r) scanner')

parser.add_argument('--ncomp', action='store', default=2, type=int, help='number of components (species) (default 2)')
parser.add_argument('--ng', action='store', default=4096, type=int, help='number of grid points (default 4096)')
parser.add_argument('--deltar', action='store', default=0.01, type=float, help='grid spacing (default 0.01)')
parser.add_argument('--alpha', action='store', default=0.2, type=float, help='Picard mixing fraction (default 0.2)')
parser.add_argument('--npic', action='store', default=6, type=int, help='number of Picard steps (default 6)')

parser.add_argument('--lb', action='store', default=1.0, type=float, help='Bjerrum length (default 1.0)')
parser.add_argument('--sigma', action='store', default=1.0, type=float, help='Gaussian charge size (default 1.0)')
parser.add_argument('--R', action='store', default=0.0, type=float, help='Grootian charge size (defaults to sigma*sqrt(15/2))')
parser.add_argument('--rc', action='store', default=1.0, type=float, help='DPD length scale (default 1.0)')
parser.add_argument('--A', action='store', default=0.0, type=float, help='DPD repulsion amplitude (default 0.0)')
parser.add_argument('--z1', action='store', default=1, type=int, help='valency of positive ions (default +1)')
parser.add_argument('--z2', action='store', default=-1, type=int, help='valency of negative ions (default -1)')
parser.add_argument('--type', action='store', default=1, type=int, help='charge type (1=Gaussian, 2=Bessel, 3=Groot, 4=Mexican)')

parser.add_argument('--rho', action='store', default=3.0, type=float, help='total density if ncomp = 3 (default 3.0)')
parser.add_argument('--npt', action='store', default=11, type=int, help='number of rho points (default 11)')
parser.add_argument('--lo', action='store', default=-2.0, type=float, help='log10(rho_min) (default -2.0)')
parser.add_argument('--hi', action='store', default=-1.0, type=float, help='log10(rho_max) (default -1.0)')

args = parser.parse_args()

w.ng = args.ng
w.ncomp = args.ncomp
w.deltar = args.deltar
w.alpha = args.alpha
w.npic = args.npic

w.initialise()

w.lb = args.lb
w.sigma = args.sigma

if (args.R > 0): w.rgroot = args.R
else: w.rgroot = args.sigma * m.sqrt(7.5)

w.rc = args.rc
w.arep[:,:] = args.A
w.z[0] = args.z1
w.z[1] = args.z2

w.dpd_potential(args.type)

n = args.npt

off = -2
eps = 1e-20

plt.figure(1)

for i in range(n):

    log10rho = args.lo + (args.hi - args.lo) * i / (n - 1.0)
    rhoz = 10.0**(log10rho)
    w.rho[0] = rhoz / (args.z1 * (args.z1 - args.z2))
    w.rho[1] = rhoz / (args.z2 * (args.z2 - args.z1))
    if (w.ncomp > 2): w.rho[2] = args.rho - w.rho[0] - w.rho[1]

    if (i == 0): w.write_params()

    w.hnc_solve()

    if (w.ncomp == 2): print('%i\t%g\t%g\t%g' % (i, log10rho, rhoz, w.error))
    else: print('%i\t%g\t%g\t%g\t%g' % (i, log10rho, rhoz, args.rho, w.error))

    plt.plot(w.r[:], list(map(lambda x, y: i*off + m.log10(eps + m.fabs(x*y)), w.hr[:, 0, 0], w.r[:])), label="%g" % log10rho)

plt.legend(loc='upper right')
plt.xlabel('$r$')
plt.ylabel('$\log_{10} r\,h_{++}$')

plt.show()
