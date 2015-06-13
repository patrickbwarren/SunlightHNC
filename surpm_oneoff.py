#!/usr/bin/python3

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Copyright (c) 2009-2014 Unilever UK Central Resources Ltd
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

parser = argparse.ArgumentParser(description='softened URPM one off calculator for structure and thermodynamics')

parser.add_argument('--ng', action='store', default=4096, type=int, help='number of grid points (default 4096)')
parser.add_argument('--deltar', action='store', default=0.01, type=float, help='grid spacing (default 0.01)')
parser.add_argument('--alpha', action='store', default=0.2, type=float, help='Picard mixing fraction (default 0.2)')
parser.add_argument('--npic', action='store', default=6, type=int, help='number of Picard steps (default 6)')

parser.add_argument('--lb', action='store', default=1.0, type=float, help='Bjerrum length (default 1.0)')
parser.add_argument('--sigma', action='store', default=1.0, type=float, help='like charge size (default 1.0)')
parser.add_argument('--sigmap', action='store', default=1.5, type=float, help='unlike charge size (default 1.5)')
parser.add_argument('--rhoz', action='store', default=0.1, type=float, help='total charge density (default 0.1)')

parser.add_argument('--grcut', action='store', default=15.0, type=float, help='r cut off for g(r) plots (default 15.0)')
parser.add_argument('--skcut', action='store', default=15.0, type=float, help='k cut off for S(k) plots (default 15.0)')

parser.add_argument('--rpa', action='store_true', help='use RPA (default HNC)')
parser.add_argument('--exp', action='store_true', help='use EXP (default HNC)')
parser.add_argument('--ushort', action='store_true', help='use U_short in potential')

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
w.rho[1] = args.rhoz / 2.0

eps = 1e-20

w.write_params()

if (args.rpa):
    w.rpa_solve()
    print('*** RPA solved')
elif (args.exp):
    w.exp_solve()
    print('*** EXP solved')
else:
    w.hnc_solve()
    print('*** HNC solved, error = ', w.error)

w.write_thermodynamics()

# code plots log10(r h(r)) versus r

plt.figure(1)

plt.subplot(2, 2, 1)

if (args.rpa): plt.title('RPA solution')
elif (args.exp): plt.title('EXP solution')
else: plt.title('HNC solution, error = %0.1g' % w.error)

plt.plot(w.r[:], 
         list(map(lambda x, y: m.log10(eps + m.fabs(x*y)), w.hr[:, 0, 0], w.r[:])), 
         label="$h_{++}$")

plt.plot(w.r[:], 
         list(map(lambda x, y: m.log10(eps + m.fabs(x*y)), w.hr[:, 0, 1], w.r[:])), 
         label="$h_{+-}$")

plt.plot(w.r[:], 
         list(map(lambda x, y: m.log10(eps + m.fabs(x*y)), w.hr[:, 1, 1], w.r[:])), 
         label="$h_{--}$")

plt.legend(loc='upper right')
plt.xlabel('$r$')
plt.ylabel('$\log_{10} r\,h$')

plt.subplot(2, 2, 2)
#plt.figure(2)

imax = int(args.grcut / w.deltar)

plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 0, 0], label="$g_{++}$")
plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 0, 1], label="$g_{+-}$")
plt.plot(w.r[0:imax], 1.0 + w.hr[0:imax, 1, 1], label="$g_{--}$")

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
