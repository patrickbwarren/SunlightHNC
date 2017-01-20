#!/usr/bin/env python3

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

######################################################################

# Given lb, the code finds the density so that the pressure matches
# the target value (which may be zero).  The code is used to generate
# phase coexistence points on the right hand (high density) side of
# the (rho, lb) phase diagram.

# The output is lb, density, pressure, chemical potential.
# The pressure should match the target value.
# The chemical potential is only correct if the HNC solution is used.

import argparse
import math as m
import scipy.optimize
from oz import wizard as w

def urpm_press(rhoz):
    w.rho[0] = rhoz / 2.0
    w.rho[1] = w.rho[0]
    if args.rpa: w.rpa_solve()
    else: w.hnc_solve()
    if args.exp: w.exp_refine()
    if (args.verbose):
        w.write_params()
        w.write_thermodynamics()
    p = w.press - args.targp
    return p

parser = argparse.ArgumentParser(description='find density for target pressure for URPM')

parser.add_argument('--ng', action='store', default=4096, type=int, help='number of grid points (default 4096)')
parser.add_argument('--deltar', action='store', default=0.01, type=float, help='grid spacing (default 0.01)')
parser.add_argument('--alpha', action='store', default=0.2, type=float, help='Picard mixing fraction (default 0.2)')
parser.add_argument('--npic', action='store', default=6, type=int, help='number of Picard steps (default 6)')

parser.add_argument('--lb', action='store', default=100.0, type=float, help='Bjerrum length (default 100.0)')
parser.add_argument('--sigma', action='store', default=1.0, type=float, help='like charge size (default 1.0)')
parser.add_argument('--rhoz', action='store', default=0.1, type=float, help='total charge density (default 0.1)')
parser.add_argument('--targp', action='store', default=0.0, type=float, help='target pressure (default 0.0)')

parser.add_argument('--rpa', action='store_true', help='use RPA (default HNC)')
parser.add_argument('--exp', action='store_true', help='use EXP refinement')

parser.add_argument('--verbose', action='store_true', help='more output')

args = parser.parse_args()

w.ng = args.ng
w.ncomp = 2
w.deltar = args.deltar
w.alpha = args.alpha
w.npic = args.npic
w.initialise()

w.lb = args.lb
w.sigma = args.sigma
w.z[0] = 1.0
w.z[1] = -1.0
w.dpd_potential()

rho1 = scipy.optimize.fsolve(urpm_press, args.rhoz)
p = urpm_press(rho1)

print("%g\t%g\t%g\t%g" % (w.lb, 2.0*w.rho[0], w.press, m.log(w.rho[0]) + w.muex[0]))

