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

# The code solves the Wertheim association model for the softened URPM
# as reference fluid. The code is used to generate phase coexistence
# points on the left hand (low density) side of the (rho, lb) phase
# diagram.

# It writes out lb, density, pressure, chemical potential, sigmap 
# where sigmap is either selected from the command line (--min not set) or
# the value at the minimum of the free energy f (--min set).

# The minimisation is performed in two ways using scipy optimise
# package.  The first is a direct attempt to minimise f(sigmap).  The
# second is an attempt to find the root df/d(sigmap) = 0.  Results can
# be found using verbose.  The final answer is pragmatically chosen as
# the average.  Currently it appears a minimum in f wrt sigmap appears
# only at low enough densities.  Otherwise it appears f is an
# increasing function of sigmap.

# Note that in calculating df/d(rho) (by finite differencing), the
# dependence on sigmap is ignored, which is correct at the minimum
# since df/d(rho) = df/d(rho) + df/d(sigmap)*d(sigmap)/d(rho) and the
# second term vanishes at the minimum.

import argparse
import math as m
import scipy.optimize
from oz import wizard as w

# Function to solve for the HNC + association free energy

def solve_wertheim(lb, rhoz, sigma, sigmap):
    w.lb = lb
    w.sigma = sigma
    w.sigmap = sigmap
    w.urpm_potential()
    w.rho[0] = rhoz / 2.0
    w.rho[1] = w.rho[0]
    w.hnc_solve()

    if (sigmap == sigma):
        x = 1.0
        rhomon = rhoz
        rhodim = 0.0
        fvass = 0.0
    else:

        # Calculate the Wertheim integral (this is no longer done in oz mod)

        d12 = 0.0

        for i in range(w.ng-1):
            g12 = 1.0 + w.hr[i, 0, 1]
            du12 = w.lb * (m.erfc(0.5*w.r[i]/w.sigma) - m.erfc(0.5*w.r[i]/w.sigmap)) / w.r[i]
            d12 = d12 + w.fourpi * w.deltar * (m.exp(-du12) - 1.0) * g12 * w.r[i]**2

        x = (m.sqrt(1.0 + 2.0*d12*rhoz) - 1.0) / (d12*rhoz)
        rhomon = x * rhoz
        rhodim = d12 * rhomon**2 / 4.0
        fvass = rhoz * m.log(x) + rhodim

    fvtot = rhoz * (m.log(rhoz) - 1.0) + w.fvex + fvass

    if (args.verbose > 1):
        print('wertheim solution at lb = ', w.lb, 'rhoz = ', 2*w.rho[0], 'sigmap = ', w.sigmap)
        print('wertheim: D+- = ', d12)
        print('wertheim: x = ', x)
        print('wertheim: rhoz = ', rhoz)
        print('wertheim: rhomon = ', rhomon)
        print('wertheim: rhodim = ', rhodim)
        print('wertheim: rhomon + 2*rhodim = ', rhomon + 2.0*rhodim)
        print('wertheim: rhodim / (0.5*rhomon)^2 = ', rhodim / (0.5*rhomon)**2)
        print('wertheim: fvass = ', fvass, 'fvass/rhoz = ', fvass/rhoz)
        print('wertheim: fvref = ', w.fvex, 'fvref/rhoz = ', w.fvex/rhoz)
        print('wertheim: fvtot = ', fvtot, 'fvtot/rhoz = ', fvtot/rhoz)

    return fvtot

# Auxiliary function returns free energy for given sigmap, used in minimisation.

def fv(sigmap):
    ans = solve_wertheim(args.lb, args.rhoz, args.sigma, sigmap)
    return ans

# Auxiliary function returns derivative by finite differencing of free
# energy with respect to sigmap, used in root finding (ie the minimum
# is df/dsigmap = 0).

def dfv(sigmap):
    fvtotp = solve_wertheim(args.lb, args.rhoz, args.sigma, sigmap*(1.0+args.wobble))
    fvtotm = solve_wertheim(args.lb, args.rhoz, args.sigma, sigmap*(1.0-args.wobble))
    ans = (fvtotp - fvtotm) / (2.0*sigmap*args.wobble)
    return ans

parser = argparse.ArgumentParser(description='Softened URPM Wertheim solver')

parser.add_argument('--ng', action='store', default=4096, type=int, help='number of grid points (default 4096)')
parser.add_argument('--deltar', action='store', default=0.01, type=float, help='grid spacing (default 0.01)')
parser.add_argument('--alpha', action='store', default=0.2, type=float, help='Picard mixing fraction (default 0.2)')
parser.add_argument('--npic', action='store', default=6, type=int, help='number of Picard steps (default 6)')

parser.add_argument('--lb', action='store', default=1.0, type=float, help='Bjerrum length (default 1.0)')
parser.add_argument('--sigma', action='store', default=1.0, type=float, help='like charge size (default 1.0)')
parser.add_argument('--sigmap', action='store', default=1.5, type=float, help='unlike charge size (default 1.5)')
parser.add_argument('--rhoz', action='store', default=0.1, type=float, help='total charge density (default 0.1)')
parser.add_argument('--wobble', action='store', default=0.02, type=float, help='wobble amount (default 0.02)')

parser.add_argument('--minimize', action='store_true', help='attempt to minimize wrt sigmap')
parser.add_argument('--verbose', action='store', default=0, type=int, help='increasing amounts of output')

args = parser.parse_args()

w.ng = args.ng
w.ncomp = 2
w.deltar = args.deltar
w.alpha = args.alpha
w.npic = args.npic

w.initialise()

if (args.minimize):
    sigmapmin1 = scipy.optimize.fmin(fv, args.sigmap, disp=0)
    if (args.verbose > 0):
        print("scipy.optimize.fmin: minimising f(sigmap) gives sigmap = %g" % (sigmapmin1))
    
    sigmapmin2 = scipy.optimize.fsolve(dfv, args.sigmap)
    if (args.verbose > 0):
        print("scipy.optimize.fsolve: df/d(sigmap)=0 gives     sigmap = %g" % (sigmapmin2))

    sigmap = 0.5*(sigmapmin1+sigmapmin2)
    if (args.verbose > 0):
        print("averaged value gives                            sigmap = %g" % (sigmap))

else:
    sigmap = args.sigmap

args.verbose = args.verbose + 1;
fvtot = solve_wertheim(args.lb, args.rhoz, args.sigma, sigmap)
args.verbose = args.verbose - 1;

if (args.verbose > 0):
    w.write_params()
    w.write_thermodynamics()

fvtotp = solve_wertheim(args.lb, args.rhoz*(1.0+args.wobble), args.sigma, sigmap)
fvtotm = solve_wertheim(args.lb, args.rhoz*(1.0-args.wobble), args.sigma, sigmap)

dfvtot = (fvtotp - fvtotm) / (2.0*args.rhoz*args.wobble)

p = args.rhoz * dfvtot - fvtot

print("%g\t%e\t%g\t%g\t%g" % (args.lb, args.rhoz, p, dfvtot, sigmap))
