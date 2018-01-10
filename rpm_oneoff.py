#!/usr/bin/env python3

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Copyright (c) 2009-2017 Unilever UK Central Resources Ltd
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
import numpy as np
import matplotlib.pyplot as plt
from oz import wizard as w

parser = argparse.ArgumentParser(description='RPM one off calculator')

parser.add_argument('--ng', action='store', default=4096, type=int, help='number of grid points (default 4096)')
parser.add_argument('--deltar', action='store', default=0.01, type=float, help='grid spacing (default 0.01)')
parser.add_argument('--alpha', action='store', default=0.2, type=float, help='Picard mixing fraction (default 0.2)')
parser.add_argument('--npic', action='store', default=6, type=int, help='number of Picard steps (default 6)')

parser.add_argument('--sigma', action='store', default=1.0, type=float, help='hard core diameter (default 1.0)')
parser.add_argument('--rho', action='store', default=0.1, type=float, help='total charge density (default 0.1)')
parser.add_argument('--lb', action='store', default=1.0, type=float, help='Bjerrum length (default 1.0)')
parser.add_argument('--kappa', action='store', default=-1.0, type=float, help='softening parameter (default off)')
parser.add_argument('--ushort', action='store_true', help='use U_short in potential')

parser.add_argument('--grcut', action='store', default=15.0, type=float, help='r cut off for g(r) plots (default 15.0)')
parser.add_argument('--skcut', action='store', default=15.0, type=float, help='k cut off for S(k) plots (default 15.0)')

parser.add_argument('--msa', action='store_true', help='use MSA instead of HNC')
parser.add_argument('--exp', action='store_true', help='use EXP (implies MSA)')
parser.add_argument('--nwarm', action='store', default=1, type=int, help='number of intermediate warm-up steps')

parser.add_argument('--dump', action='store_true', help='write out g(r)')
parser.add_argument('--show', action='store_true', help='plot results')

parser.add_argument('--verbose', action='store_true', help='more output')

args = parser.parse_args()

w.ncomp = 2
w.ng = args.ng
w.deltar = args.deltar
w.alpha = args.alpha
w.npic = args.npic

w.initialise()

# print(w.c.shape)

w.sigma = args.sigma
w.kappa = args.kappa

w.rho[0] = 0.5 * args.rho
w.rho[1] = 0.5 * args.rho

if args.verbose:
    w.write_params()
    w.verbose = True

if args.exp:    
    w.hs_potential()
    w.msa_solve()
    w.save_reference()
    w.backup()
    s = str(w.closure_name, 'utf-8').strip()
    print('rho = %g \thard spheres \t%s error = %g' % (np.sum(w.rho), s, w.error))
    args.msa = True

for i in range(args.nwarm):
    w.lb = (i + 1.0) / args.nwarm * args.lb
    w.rpm_potential(args.ushort)
    if args.msa:
        if args.exp: w.restore()
        w.msa_solve()
    else:
        w.hnc_solve()
    s = str(w.closure_name, 'utf-8').strip()
    print('rho = %g \tlb = %g \tkappa = %g \t%s error = %g' % (np.sum(w.rho), w.lb, w.kappa, s, w.error))
        
if args.exp:
    w.backup()
    w.exp_refine()
    s = str(w.closure_name, 'utf-8').strip()
    print('rho = %g \tlb = %g \tkappa = %g \t%s error = %g' % (np.sum(w.rho), w.lb, w.kappa, s, w.error))

if w.return_code: exit()

if not args.dump: w.write_thermodynamics()

# density-density and charge-charge structure factor
# (notice how elegant this is :-)

snn = np.sum(np.sum(w.sk, axis=2), axis=1) / np.sum(w.rho)
szz = np.dot(np.dot(w.z, w.sk), w.z) / np.sum(w.rho)

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
        print("%g\t%g\t%g\t%g\tH" % (w.r[i], w.hr[i, 0, 0],
                                     w.hr[i, 0, 1], w.hr[i, 1, 1]))

    for i in range(w.ng-1):
        print("%g\t%g\t%g\tS" % (w.k[i], snn[i], szz[i]))

elif args.show:

    plt.figure(1)

    plt.subplot(2, 2, 1)

    plt.title('%s solution, error = %0.1g' % (str(w.closure_name, 'utf-8'),
                                              w.error))

    imax = int(args.grcut / w.deltar)
    plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 0, 0], label="$g_{11}(r)$")
    plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 0, 1], label="$g_{12}(r)$")
    plt.legend(loc='lower right')
    
    plt.subplot(2, 2, 2)

    jmax = int(args.skcut / w.deltak)
    plt.plot(w.k[:jmax], snn[:jmax], label='$S_{NN}(k)$')
    plt.plot(w.k[:jmax], szz[:jmax], label='$S_{ZZ}(k)$')
    plt.legend(loc='lower right')

    plt.subplot(2, 2, 3)
    
    plt.plot(w.r[0:imax],
             0.5*(w.c[0:imax, 0, 0]-w.ulong[0:imax, 0]+w.c[0:imax, 1, 0]
                  -w.ulong[0:imax, 1]),
             label="$[c_{11}+c_{12}]/2$")
    plt.plot(w.r[0:imax],
             0.5*(w.c[0:imax, 0, 0]-w.ulong[0:imax, 0]-w.c[0:imax, 1, 0]
                  +w.ulong[0:imax, 1]),
             label="$[c_{11}-c_{12}]/2$")
    plt.legend(loc='lower right')
    
    plt.subplot(2, 2, 4)

    jmax = int(args.skcut*3 / w.deltak)
    plt.plot(w.k[0:jmax],w.ek[0:jmax, 0]+w.ck[0:jmax, 0],
             label="${\\tilde h}_{11}(k)$")
    plt.plot(w.k[0:jmax],w.ek[0:jmax, 1]+w.ck[0:jmax, 1],
             label="${\\tilde h}_{12}(k)$")
    plt.plot(w.k[0:jmax],
             1.0 + w.rho[0]*(w.ek[0:jmax, 0]+w.ck[0:jmax, 0]
                             -w.ek[0:jmax, 1]-w.ck[0:jmax, 1]), label="$\\ne 0$ !") 
#    plt.plot(w.k[0:jmax],w.ck[0:jmax, 0],label="${\\tilde c}_{11}$")
#    plt.plot(w.k[0:jmax],w.ck[0:jmax, 1],label="${\\tilde c}_{12}$")
    plt.legend(loc='lower right')
    
    plt.show()
