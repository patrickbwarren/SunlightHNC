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
import matplotlib.pyplot as plt
from oz import wizard as w

parser = argparse.ArgumentParser(description='RPM one off calculator')

parser.add_argument('--ng', action='store', default=65536, type=int, help='number of grid points (default 4096)')
parser.add_argument('--deltar', action='store', default=0.001, type=float, help='grid spacing (default 0.01)')
parser.add_argument('--alpha', action='store', default=0.2, type=float, help='Picard mixing fraction (default 0.2)')
parser.add_argument('--npic', action='store', default=6, type=int, help='number of Picard steps (default 6)')

parser.add_argument('--sigma', action='store', default=1.0, type=float, help='hard core diameter (default 1.0)')

parser.add_argument('--eta', action='store', default=0.3, type=float, help='packing fraction (default 0.3)')
parser.add_argument('--grcut', action='store', default=15.0, type=float, help='r cut off for g(r) plots (default 15.0)')
parser.add_argument('--skcut', action='store', default=15.0, type=float, help='k cut off for S(k) plots (default 15.0)')

parser.add_argument('--msa', action='store_true', help='use MSA (default HNC)')

parser.add_argument('--dump', action='store_true', help='write out g(r)')

parser.add_argument('--verbose', action='store_true', help='more output')

args = parser.parse_args()

w.ncomp = 1
w.ng = args.ng
w.deltar = args.deltar
w.alpha = args.alpha
w.npic = args.npic

w.initialise()

# print(w.c.shape)

w.sigma = args.sigma

w.hs_potential()

w.rho[0] = 6.0 * args.eta / (w.pi * w.sigma**3)

eps = 1e-20

if args.verbose:
    w.write_params()
    w.verbose = 1

if args.msa: w.msa_solve()
else: w.hnc_solve()

if not args.dump:
    if args.msa: print('*** MSA solved, error = ', w.error)
    else: print('*** HNC solved, error = ', w.error)
    w.write_thermodynamics()

if args.dump:

    for i in range(w.ng-1):
        print("%g\t%g\tC" % (w.r[i], w.c[i, 0, 0]))

    for i in range(w.ng-1):
        print("%g\t%g\tH" % (w.r[i], w.hr[i, 0, 0]))

    for i in range(w.ng-1):
        print("%g\t%g\tS" % (w.k[i], w.sk[i, 0, 0]))

else:

    plt.figure(1)

    plt.subplot(2, 2, 1)

    if (args.msa): plt.title('MSA solution, error = %0.1g' % w.error)
    else: plt.title('HNC solution, error = %0.1g' % w.error)

    imax = int(args.grcut / w.deltar)
    plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 0, 0], label="$g(r)$")
    plt.legend(loc='lower right')
    
    plt.subplot(2, 2, 2)

    jmax = int(args.skcut / w.deltak)
    plt.plot(w.k[:jmax], w.sk[:jmax, 0, 0], label='$S(k)$')
    plt.legend(loc='lower right')

    plt.subplot(2, 2, 3)
    
    plt.plot(w.r[0:imax], w.c[0:imax, 0, 0], label="$c(r)$")
    plt.legend(loc='lower right')
    
    plt.subplot(2, 2, 4)

    jmax = int(args.skcut*3 / w.deltak)
    plt.plot(w.k[0:jmax], w.ek[0:jmax, 0]+w.ck[0:jmax, 0], label="${\\tilde h}(k)$")
    plt.legend(loc='lower right')
    
    plt.show()

