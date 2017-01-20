#!/usr/bin/env python3

# The code writes lb, density, pressure, chemical potential.
# The chemical potential only makes sense for the HNC solution.

import argparse
import math as m
from oz import wizard as w

parser = argparse.ArgumentParser(description='URPM thermodynamics generator')

parser.add_argument('--ng', action='store', default=4096, type=int, help='number of grid points (default 4096)')
parser.add_argument('--deltar', action='store', default=0.01, type=float, help='grid spacing (default 0.01)')
parser.add_argument('--alpha', action='store', default=0.2, type=float, help='Picard mixing fraction (default 0.2)')
parser.add_argument('--npic', action='store', default=6, type=int, help='number of Picard steps (default 6)')

parser.add_argument('--lb', action='store', default=1.0, type=float, help='Bjerrum length (default 1.0)')
parser.add_argument('--sigma', action='store', default=1.0, type=float, help='like charge size (default 1.0)')
parser.add_argument('--rhoz', action='store', default=0.1, type=float, help='total charge density (default 0.1)')

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

w.rho[0] = args.rhoz / 2.0
w.rho[1] = w.rho[0]

if args.rpa: w.rpa_solve()
else: w.hnc_solve()

if args.exp: w.exp_refine()

if (args.verbose):
    w.write_params()
    w.write_thermodynamics()

print("%g\t%g\t%g\t%g" % (w.lb, 2.0*w.rho[0], w.press, m.log(w.rho[0]) + w.muex[0]))

