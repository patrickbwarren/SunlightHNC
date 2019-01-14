#!/usr/bin/env python3

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

# Note on diameter setting in --diam
# first three entries are HS diameters z = +1, -1, 0
# [1,1,1] is equivalent to default, [1,1,0.8] is solvent of smaller spheres
# next three entries are the pairwise (+-), (+0), (-0) hard core distances
# [1,1,0] would be RPM + gas of points (not much interest)
# [1,1,0,1,0.9] is RPM + excluded volume points (Asakura-Oosawa model)
# or equivalently --diam="[1,1,0,1,(1+0.8)/2]" at q = 0.8 (size ratio)
# can also set --sigma=1 with this as the charged spheres have unit diameter

# The following reproduce Fig 1 insets in Coupette et al., PRL 121, 075501 (2018)
# For this lB = 0.7 nm, sigma = 0.3 nm, [salt] = 1 M, and [solvent] = 10 M and 40 M.
# Note the use of computed values for lB/sigma and rho*sigma^3.  Also 1 M = 0.602 molecules per nm^3
# For the salt, NaCl --> Na+ and Cl-, and rhoz = [Na+] + [Cl-] = 2 [NaCl], hence the factor 2 in --rhoz

# python3 rpm_oneoff.py --ncomp=3 --lb=0.7/0.3 --rhoz=2*0.602*0.3^3 --add --rho=10*0.602*0.3^3 --show --all --tail --only
# python3 rpm_oneoff.py --ncomp=3 --lb=0.7/0.3 --rhoz=2*0.602*0.3^3 --add --rho=40*0.602*0.3^3 --show --all --tail --only

# Add --diam='[0.25/0.3,0.3373/0.3,1]' to reproduce the size-asymmetric model shown in Fig S1.
# Note though that the + and - are the wrong way around in this figure.

# Replace "--show --all --tail --only" by "--dump", to get to the raw data, then filter on 'H' to extract h_ij(r)

import argparse
import numpy as np
from oz import wizard as w

parser = argparse.ArgumentParser(description='RPM one off calculator')

parser.add_argument('--ncomp', action='store', default=2, type=int, help='number of components (species) (default 2)')
parser.add_argument('--ng', action='store', default='65536', help='number of grid points (default 2^16 = 65536)')
parser.add_argument('--deltar', action='store', default=1e-3, type=float, help='grid spacing (default 1e-3)')
parser.add_argument('--alpha', action='store', default=0.2, type=float, help='Picard mixing fraction (default 0.2)')
parser.add_argument('--npic', action='store', default=6, type=int, help='number of Picard steps (default 6)')
parser.add_argument('--nps', action='store', default=6, type=int, help='length of history array (default 6)')
parser.add_argument('--maxsteps', action='store', default=100, type=int, help='number of iterations (default 100)')

parser.add_argument('--diam', action='store', default='[1]', help='hard core diameters (default [1])')
parser.add_argument('--sigma', action='store', default=0.0, type=float, help='inner core diameter (default min diam)')
parser.add_argument('--rho', action='store', default='0.5', help='total hard sphere density (default 0.5)')
parser.add_argument('--rhoz', action='store', default='0.1', help='total ion density (default 0.1)')
parser.add_argument('--lb', action='store', default='1.0', help='Bjerrum length (default 1.0)')
parser.add_argument('--kappa', action='store', default=-1.0, type=float, help='softening parameter (default off)')
parser.add_argument('--ushort', action='store_true', help='use U_short in potential')

parser.add_argument('--grcut', action='store', default=15.0, type=float, help='r cut off for g(r) plots (default 15.0)')
parser.add_argument('--skcut', action='store', default=15.0, type=float, help='k cut off for S(k) plots (default 15.0)')

parser.add_argument('--msa', action='store_true', help='use MSA instead of HNC')
parser.add_argument('--exp', action='store_true', help='use EXP (implies MSA)')
parser.add_argument('--nwarm', action='store', default=1, type=int, help='number of intermediate warm-up steps')

parser.add_argument('--dump', action='store_true', help='write out g(r)')
parser.add_argument('--show', action='store_true', help='plot results')

parser.add_argument('--eps', action='store', default=1e-20, type=float, help='floor for log tails (default 1e-20)')
parser.add_argument('--only', action='store_true', help='plot only pair functions')
parser.add_argument('--tail', action='store_true', help='plot tails of pair functions')
parser.add_argument('--all', action='store_true', help='include solvent pair functions')
parser.add_argument('--add', action='store_true', help='add solute rather than make up as total')

parser.add_argument('--verbose', action='store_true', help='more output')

args = parser.parse_args()

w.ncomp = args.ncomp
w.ng = eval(args.ng.replace('^', '**')) # catch 2^10 etc
w.deltar = args.deltar
w.alpha = args.alpha
w.npic = args.npic
w.nps = args.nps
w.maxsteps = args.maxsteps
w.verbose = args.verbose

w.initialise()

w.sigma = args.sigma
w.kappa = args.kappa

# Now construct the hard core diameters

diam = eval(args.diam)

if not isinstance(diam, list):
    diam = [diam]

if len(diam) == 1: diam.append(diam[0])

# These are for the charged species

w.diam[0, 0] = diam[0]
w.diam[0, 1] = 0.5*(diam[0] + diam[1])
w.diam[1, 1] = diam[1]

if w.ncomp == 3: # solvent spheres
    if len(diam) == 2: diam.append(0.5*(diam[0] + diam[1]))
    w.diam[0, 2] = 0.5*(diam[0] + diam[2])
    w.diam[1, 2] = 0.5*(diam[1] + diam[2])
    w.diam[2, 2] = diam[2]

# Now over-write pairwise diameters for non-additivity

if len(diam) == 4: w.diam[0, 1] = diam[3]
if len(diam) == 5: w.diam[0, 2] = w.diam[1, 2] = diam[4]
if len(diam) == 6: w.diam[1, 2] = diam[5]

rhoz = eval(args.rhoz.replace('^', '**')) # total charged species density
rho = eval(args.rho.replace('^', '**')) # total density, or solvent density with --add

w.rho[0] = 0.5 * rhoz
w.rho[1] = 0.5 * rhoz

if w.ncomp == 3:
    if args.add:
        w.rho[2] = rho
    else:
        w.rho[2] = rho - rhoz

if args.exp:
    w.hs_potential()
    w.msa_solve()
    w.save_reference()
    s = str(w.closure_name, 'utf-8').strip()
    print('rho = %g \thard spheres \t%s error = %g' % (np.sum(w.rho), s, w.error))
    args.msa = True

for i in range(args.nwarm):
    w.lb = (i + 1.0) / args.nwarm * eval(args.lb)
    w.rpm_potential(args.ushort)
    if w.verbose: w.write_params()
    if args.msa:
        w.msa_solve()
    else:
        w.hnc_solve()
    s = str(w.closure_name, 'utf-8').strip()
    print('rhoz = %g \trho = %g \tlb = %g \tkappa = %g \t%s error = %g' %
          (w.rho[0]+w.rho[1], np.sum(w.rho), w.lb, w.kappa, s, w.error))

if args.exp:
    w.exp_refine()
    s = str(w.closure_name, 'utf-8').strip()
    print('rhoz = %g \trho = %g \tlb = %g \tkappa = %g \t%s error = %g' %
          (w.rho[0]+w.rho[1], np.sum(w.rho), w.lb, w.kappa, s, w.error))

if w.return_code: exit()

if not args.dump:
    w.write_params()
    w.write_thermodynamics()

# density-density and charge-charge structure factor
# (notice how elegant this is :-)

snn = np.sum(np.sum(w.sk, axis=2), axis=1) / np.sum(w.rho)
szz = np.dot(np.dot(w.z, w.sk), w.z) / np.sum(w.rho * w.z**2)

if args.dump:

    for i in range(w.ng-1): 
        print("%g\t%g\t%g\t%g\tL" % (w.r[i], w.ulong[i, 0],
                                     w.ulong[i, 1], w.ulong[i, 2]))

    for i in range(w.ng-1):
        print("%g\t%g\t%g\tC" % (w.r[i],
                                     0.5*(w.c[i, 0, 0]-w.ulong[i, 0]
                                          +w.c[i, 1, 0]-w.ulong[i, 1]),
                                     0.5*(w.c[i, 0, 0]-w.ulong[i, 0]
                                          -w.c[i, 1, 0]+w.ulong[i, 1])))

    for i in range(w.ng-1):
        if w.ncomp == 2:
            print("%g\t%g\t%g\t%g\tH" % (w.r[i], w.hr[i, 0, 0],
                                         w.hr[i, 0, 1], w.hr[i, 1, 1]))
        else:
            print("%g\t%g\t%g\t%g\t%g\t%g\t%g\tH" % (w.r[i], w.hr[i, 0, 0], w.hr[i, 0, 1], w.hr[i, 1, 1],
                                                     w.hr[i, 0, 2], w.hr[i, 1, 2], w.hr[i, 2, 2]))
            
    for i in range(w.ng-1):
        print("%g\t%g\t%g\tS" % (w.k[i], snn[i], szz[i]))

elif args.show:

    import math as m
    import matplotlib.pyplot as plt

    def plot_gr():
        "Plot pair correlation functions"
        imax = int(args.grcut / w.deltar)
        plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 0, 0], label="$g_{11}(r)$")
        plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 0, 1], label="$g_{12}(r)$")
        plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 1, 1], label="$g_{22}(r)$")
        if args.all:
            plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 0, 2], label="$g_{13}(r)$")
            plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 1, 2], label="$g_{23}(r)$")
            plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 2, 2], label="$g_{33}(r)$")
        plt.legend(loc='lower right')
    
    def plot_sk():
        "Plot structure factors"
        jmax = int(args.skcut / w.deltak)
        plt.plot(w.k[:jmax], snn[:jmax], label='$S_{NN}(k)$')
        plt.plot(w.k[:jmax], szz[:jmax], label='$S_{ZZ}(k)$')
        plt.legend(loc='lower right')
    
    def plot_hk():
        "Plot fourier transformed pair functions"
        jmax = int(args.skcut*3 / w.deltak)
        plt.plot(w.k[0:jmax],w.ek[0:jmax, 0]+w.ck[0:jmax, 0], label="${\\tilde h}_{11}(k)$")
        plt.plot(w.k[0:jmax],w.ek[0:jmax, 1]+w.ck[0:jmax, 1], label="${\\tilde h}_{12}(k)$")
        plt.legend(loc='lower right')

    def plot_cr():
        "Plot direct correlation functions"
        imax = int(args.grcut / w.deltar)
        plt.plot(w.r[0:imax],
                 0.5*(w.c[0:imax, 0, 0]-w.ulong[0:imax, 0]+w.c[0:imax, 1, 0]
                      -w.ulong[0:imax, 1]),
                 label="$[c_{11}+c_{12}]/2$")
        plt.plot(w.r[0:imax],
                 0.5*(w.c[0:imax, 0, 0]-w.ulong[0:imax, 0]-w.c[0:imax, 1, 0]
                      +w.ulong[0:imax, 1]),
                 label="$[c_{11}-c_{12}]/2$")
        plt.legend(loc='lower right')

    def plot_rhrtail():
        "Plot log r|h| versus r, to show tails"
        sumrhoz2 = np.sum(w.rho * w.z**2) # ionic strength
        if sumrhoz2 > 0 and w.lb > 0:
            lD = 1.0 / m.sqrt(4*w.pi*w.lb*sumrhoz2) # Debye length
        else:
            lD = 0
        print("ionic strength sum rho z^2 = %g" % sumrhoz2)
        print("Debye length lD = %g" % lD)
        plt.plot(w.r[:],
                 list(map(lambda x, y: m.log10(args.eps + m.fabs(x*y)), w.hr[:, 0, 0], w.r[:])),
                 label="$r|h_{11}|$")
        plt.plot(w.r[:],
                 list(map(lambda x, y: m.log10(args.eps + m.fabs(x*y)), w.hr[:, 0, 1], w.r[:])),
                 label="$r|h_{12}|$")
        plt.plot(w.r[:],
                 list(map(lambda x, y: m.log10(args.eps + m.fabs(x*y)), w.hr[:, 1, 1], w.r[:])),
                 label=" $r|h_{22}|$")
        if lD > 0:
            plt.plot(w.r[:],
                     list(map(lambda x: m.log10(args.eps + m.exp(-x/lD)), w.r[:])),
                     label=" $e^{-\kappa r}$", linestyle='--', color='black')
        if args.all:
            plt.plot(w.r[:],
                     list(map(lambda x, y: m.log10(args.eps + m.fabs(x*y)), w.hr[:, 0, 2], w.r[:])),
                     label="$r|h_{13}|$")
            plt.plot(w.r[:],
                     list(map(lambda x, y: m.log10(args.eps + m.fabs(x*y)), w.hr[:, 1, 2], w.r[:])),
                     label="$r|h_{23}|$")
            plt.plot(w.r[:],
                     list(map(lambda x, y: m.log10(args.eps + m.fabs(x*y)), w.hr[:, 2, 2], w.r[:])),
                     label="$r|h_{33}|$")
        plt.legend(loc='upper right')

    # Here's where the plotting starts

    plt.figure(1)

    if args.only: # only do one plot depending on other arguments
        
        plt.title('%s solution, error = %0.1g' % (str(w.closure_name, 'utf-8'), w.error))
        if args.tail: # plot log10(r h(r)) versus r
            plot_rhrtail()
        else: # plot pair functions
            plot_gr()

    else: # do four subplots

        plt.subplot(2, 2, 1) # pair functions
        plt.title('%s solution, error = %0.1g' % (str(w.closure_name, 'utf-8'), w.error))
        plot_gr()
        plt.subplot(2, 2, 2) # structure factors
        plot_sk()
        plt.subplot(2, 2, 3) # direct correlation functions
        plot_cr()
        plt.subplot(2, 2, 4)
        if args.tail: # plot log10(r h(r)) versus r
            plot_rhrtail()
        else: # plot total correlation functions in reciprocal space
            plot_hk()
        
    plt.show()

