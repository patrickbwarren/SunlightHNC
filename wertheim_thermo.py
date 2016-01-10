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

import argparse
import math as m
import numpy as np
import matplotlib.pyplot as plt
from oz import wizard as w

# Function calculates the excess free energy per particle by coupling
# constant integration along an isochore.

def fnex1(nlb):
    fnex_xc = prev = 0.0
    dlb = w.lb / nlb
    n = int(w.lb/dlb + 0.5)
    w.cold_start = 1
    for i in range(n):
        w.lb = dlb * (i + 1.0)
        w.soft_urpm_potential(0)
        w.hnc_solve()
        curr = w.un_xc / w.lb
        fnex_xc = fnex_xc + 0.5*dlb*(prev + curr)
        prev = curr
#        print('fnex1: w.lb = ', w.lb
    return w.un_mf + fnex_xc

# Function calculates the excess free energy per particle by
# integrating the chemical potential along an isotherm.  Note that the
# calculation relies on rho1 = rho2 = rho/2.

def fnex2(nrho):
    fvex = prev = 0.0
    drho = 2.0*w.rho[0]/nrho
    n = int(2.0*w.rho[0]/drho + 0.5)
    w.cold_start = 1
    for i in range(n):
        w.rho[0] = drho * (i + 1.0) / 2.0
        w.rho[1] = w.rho[0]
        w.hnc_solve()
        curr = w.muex[0]
        fvex = fvex + 0.5*drho*(prev + curr)
        prev = curr
#        print('fnex2: w.rho[0]+w.rho[1] = ', w.rho[0] + w.rho[1]
    return fvex / (2.0*w.rho[0])

parser = argparse.ArgumentParser(description='softened URPM one off calculator for Wertheim thermodynamics')

parser.add_argument('--ng', action='store', default=4096, type=int, help='number of grid points (default 4096)')
parser.add_argument('--deltar', action='store', default=0.01, type=float, help='grid spacing (default 0.01)')
parser.add_argument('--alpha', action='store', default=0.2, type=float, help='Picard mixing fraction (default 0.2)')
parser.add_argument('--npic', action='store', default=6, type=int, help='number of Picard steps (default 6)')

parser.add_argument('--lb', action='store', default=1.0, type=float, help='Bjerrum length (default 1.0)')
parser.add_argument('--sigma', action='store', default=1.0, type=float, help='like charge size (default 1.0)')
parser.add_argument('--sigmap', action='store', default=1.5, type=float, help='unlike charge size (default 1.5)')
parser.add_argument('--rhoz', action='store', default=0.1, type=float, help='total charge density (default 0.1)')

parser.add_argument('--rpa', action='store_true', help='use RPA (default HNC)')
parser.add_argument('--exp', action='store_true', help='use EXP (default HNC)')
parser.add_argument('--ushort', action='store_true', help='use U_short in potential')
parser.add_argument('--test', action='store', default=0, type=int, help='select a test to perform (see code for details)')

args = parser.parse_args()

w.ng = args.ng
w.ncomp = 2
w.deltar = args.deltar
w.alpha = args.alpha
w.npic = args.npic

w.initialise()

w.lb = args.lb
w.sigma = args.sigma
w.sigmap = args.sigmap

if (args.ushort): w.soft_urpm_potential(1)
else: w.soft_urpm_potential(0)

w.rho[0] = args.rhoz / 2.0
w.rho[1] = w.rho[0]

eps = 1e-20

# w.write_params()

if (args.rpa):
    w.rpa_solve()
    print('*** RPA solved')
elif (args.exp):
    w.exp_solve()
    print('*** EXP solved')
else:
    w.hnc_solve()
    print('*** HNC solved, error = ', w.error)

# w.write_thermodynamics()

if (args.test == 1):
    print('fnex (direct)       =', w.fnex)
    print('fnex (energy route) =', fnex1(100))
    print('fnex (mu route)     =', fnex2(100))

rhotot = args.rhoz
x = (m.sqrt(1.0 + 2.0*w.d12*rhotot) - 1.0) / (w.d12*rhotot)
rhomon = x * rhotot
rhodim = w.d12 * rhomon**2 / 4.0

print('Wertheim D+- = ', w.d12)
print('Wertheim x = ', x)
print('Total density = ', rhotot)
print('Monomer density = ', rhomon)
print('Dimer density = ', rhodim)
print('Mass balance monomer + 2*dimer = ', rhomon + 2.0*rhodim)
print('Association constant dimer / (0.5 monomer)^2 = ', rhodim / (0.5*rhomon)**2)

fvass = rhotot * m.log(x) + rhodim
fvtot = w.fvex + fvass

print('Wertheim association free energy density = ', fvass)
print('Reference fluid free energy density = ', w.fvex)
print('Total free energy density = ', fvtot)
print('Total free energy per particle = ', fvtot / rhotot)

