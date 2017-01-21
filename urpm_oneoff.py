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

import argparse
import math as m
import numpy as np
from oz import wizard as w

parser = argparse.ArgumentParser(description='URPM one off calculator')

parser.add_argument('--ncomp', action='store', default=2, type=int, help='number of components (species) (default 2)')
parser.add_argument('--ng', action='store', default=4096, type=int, help='number of grid points (default 4096)')
parser.add_argument('--deltar', action='store', default=0.01, type=float, help='grid spacing (default 0.01)')
parser.add_argument('--alpha', action='store', default=0.2, type=float, help='Picard mixing fraction (default 0.2)')
parser.add_argument('--npic', action='store', default=6, type=int, help='number of Picard steps (default 6)')

parser.add_argument('--lb', action='store', default=1.0, type=float, help='Bjerrum length (default 1.0)')
parser.add_argument('--sigma', action='store', default=1.0, type=float, help='Gaussian/Bessel charge size (default 1.0)')
parser.add_argument('--lambda', action='store', default=1.0, dest='lbda', type=float, help='Slater charge size (default 1.0)')
parser.add_argument('--R', action='store', default=1.0, dest='rgroot', type=float, help='Groot charge size (defaults 1.0)')
parser.add_argument('--rc', action='store', default=1.0, type=float, help='DPD length scale (default 1.0)')
parser.add_argument('--A', action='store', default=0.0, type=float, help='DPD repulsion amplitude (default 0.0)')
parser.add_argument('--z1', action='store', default=1, type=int, help='valency of positive ions (default +1)')
parser.add_argument('--z2', action='store', default=-1, type=int, help='valency of negative ions (default -1)')
parser.add_argument('--type', action='store', default=1, type=int, help='charge type (1=Gaussian, 2=Bessel, 3=Groot, 4=Slater)')
parser.add_argument('--case', action='store', default=1, type=int, help='Slater method (1=exact, 2=good, 3=bad)')

parser.add_argument('--rho', action='store', default=3.0, type=float, help='total density if ncomp = 3 (default 3.0)')
parser.add_argument('--rhoz', action='store', default=0.1, type=float, help='total charge density (default 0.1)')
parser.add_argument('--grcut', action='store', default=15.0, type=float, help='r cut off for g(r) plots (default 15.0)')
parser.add_argument('--skcut', action='store', default=15.0, type=float, help='k cut off for S(k) plots (default 15.0)')

parser.add_argument('--rpa', action='store_true', help='use RPA (default HNC)')
parser.add_argument('--exp', action='store_true', help='use EXP refinement')

args = parser.parse_args()

w.ng = args.ng
w.ncomp = args.ncomp
w.deltar = args.deltar
w.alpha = args.alpha
w.npic = args.npic

w.initialise()

w.lb = args.lb
w.sigma = args.sigma
w.lbda = args.lbda
w.rgroot = args.rgroot

w.rc = args.rc
w.arep[:,:] = args.A
w.z[0] = args.z1
w.z[1] = args.z2

# potential type = 4 (exact), or potential type = 5 (approximate) with
# beta = 1/lambda, or beta=5/(8lambda)

if args.type < 4:
    w.dpd_potential(args.type)
else:
    if args.case == 1:
        w.dpd_potential(w.dpd_slater_exact_charges)
    else:
        if args.case == 2: w.beta = 5 / (8*w.lbda)
        else: w.beta = 1 / w.lbda
        w.dpd_potential(w.dpd_slater_approx_charges)

# The calculation here solves rhoz = z1^2*rho1 + z2^2*rho2, z1*rho1 + z2*rho2 = 0

w.rho[0] = args.rhoz / (args.z1 * (args.z1 - args.z2))
w.rho[1] = args.rhoz / (args.z2 * (args.z2 - args.z1))

if (w.ncomp > 2): w.rho[2] = args.rho - w.rho[0] - w.rho[1]

eps = 1e-20

w.write_params()

if (args.rpa or args.exp): w.rpa_solve()
else: w.hnc_solve()

if args.exp: w.exp_refine()

if w.return_code: exit()

w.write_thermodynamics()

# code plots log10(r h(r)) versus r

import matplotlib.pyplot as plt

plt.figure(1)

plt.subplot(2, 2, 1)

plt.title('%s solution, error = %0.1g' % (str(w.closure_name, 'utf-8'), w.error))

plt.plot(w.r[:], 
         list(map(lambda x, y: m.log10(eps + m.fabs(x*y)), w.hr[:, 0, 0], w.r[:])), 
         label="$+\!+$")

plt.plot(w.r[:], 
         list(map(lambda x, y: m.log10(eps + m.fabs(x*y)), w.hr[:, 0, 1], w.r[:])), 
         label="$+ -$")

plt.plot(w.r[:], 
         list(map(lambda x, y: m.log10(eps + m.fabs(x*y)), w.hr[:, 1, 1], w.r[:])), 
         label=" $- -$")

if (w.ncomp == 3):

    plt.plot(w.r[:], 
             list(map(lambda x, y: m.log10(eps + m.fabs(x*y)), w.hr[:, 0, 2], w.r[:])), 
             label="$+0$")

    plt.plot(w.r[:], 
             list(map(lambda x, y: m.log10(eps + m.fabs(x*y)), w.hr[:, 1, 2], w.r[:])), 
             label="$-\,0$")

    plt.plot(w.r[:], 
             list(map(lambda x, y: m.log10(eps + m.fabs(x*y)), w.hr[:, 2, 2], w.r[:])), 
             label="$0\,0$")

plt.legend(loc='upper right')
plt.xlabel('$r$')
plt.ylabel('$\log_{10}|rh|$')

plt.subplot(2, 2, 2)

imax = int(args.grcut / w.deltar)

plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 0, 0], label="$+\!+$")
plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 0, 1], label="$+ -$")
plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 1, 1], label=" $- -$")

if (w.ncomp == 3):
    plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 0, 2], label="$+0$")
    plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 1, 2], label="$-\,0$")
    plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 2, 2], label="$0\,0$")

plt.legend(loc='upper right')
plt.xlabel('$r$')
plt.ylabel('$g(r)$')

# density-density structure factor

ddsf = np.sum(np.sum(w.sk, axis=2), axis=1) / np.sum(w.rho)

# charge-charge structure factor (notice how elegant this is :-)

ccsf = np.dot(np.dot(w.z, w.sk), w.z) / np.sum(w.rho)

plt.subplot(2, 2, 3)

jmax = int(args.skcut / w.deltak)
plt.plot(w.k[:jmax], ddsf[:jmax], label='$S_{NN}$')
plt.plot(w.k[:jmax], ccsf[:jmax], label='$S_{ZZ}$')
plt.legend(loc='lower right')
plt.xlabel('$k$')
plt.ylabel('$S(k)$')

plt.show()
