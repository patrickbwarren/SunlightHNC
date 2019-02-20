#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This file includes unicode characters

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

# By default (--ncomp=2) this is for RPM without solvent
# Run with --ncomp=3 for solvent primitive model

# The following correspond to the Fig 1 insets in Coupette et al., PRL 121, 075501 (2018)
# For this lB = 0.7 nm, sigma = 0.3 nm, [salt] = 1 M, and [solvent] = 10 M and 40 M.
# Note the use of computed values for lB/sigma and rho*sigma^3.  Also 1 M = 0.602 molecules per nm^3
# For the salt, NaCl --> Na+ and Cl-, and rhoz = [Na+] + [Cl-] = 2 [NaCl], hence the factor 2 in --rhoz

# python3 rpm_explorer.py --ncomp=3 --lb=0.7/0.3 --rhoz=2*0.602*0.3^3 --rhos=10*0.602*0.3^3
# python3 rpm_explorer.py --ncomp=3 --lb=0.7/0.3 --rhoz=2*0.602*0.3^3 --rhos=40*0.602*0.3^3

# Add --diam='[0.25/0.3,0.3373/0.3,1]' to reproduce the size-asymmetric model shown in Fig S1.

import argparse
import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

from oz import wizard as w
from numpy import pi as π

parser = argparse.ArgumentParser(description='Interactive RPM explorer')

parser.add_argument('--ncomp', action='store', default=2, type=int, help='number of components (species) (default 2)')

parser.add_argument('--ng', action='store', default='16384', help='number of grid points (default 2^14 = 16384)')
parser.add_argument('--deltar', action='store', default=1e-3, type=float, help='grid spacing (default 1e-3)')

parser.add_argument('--alpha', action='store', default=0.2, type=float, help='Picard mixing fraction (default 0.2)')
parser.add_argument('--npic', action='store', default=6, type=int, help='number of Picard steps (default 6)')
parser.add_argument('--nps', action='store', default=6, type=int, help='length of history array (default 6)')
parser.add_argument('--maxsteps', action='store', default=100, type=int, help='number of iterations (default 100)')

parser.add_argument('--diam', action='store', default='[1]', help='hard core diameters (default [1])')
parser.add_argument('--sigma', action='store', default=0.0, type=float, help='inner core diameter (default min diam)')
parser.add_argument('--rhoz', action='store', default='0.1', help='total ion density (default 0.1)')
parser.add_argument('--rhos', action='store', default='0.4', help='added solvent density (default 0.4)')
parser.add_argument('--lb', action='store', default='1.0', help='Bjerrum length (default 1.0)')
parser.add_argument('--swap', action='store_true', help='swap representations of h(r)')
parser.add_argument('--neg', action='store_true', help='invert sign of h(r)')

parser.add_argument('--rmax', action='store', default=15.0, type=float, help='maximum radial distance')

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

w.lb = lb_init = eval(args.lb)

# Now construct the hard core diameters

diam = eval(args.diam)

if not isinstance(diam, list):
    diam = [diam]

if len(diam) == 1: diam.append(diam[0])

w.diam[0, 0] = diam[0]
w.diam[0, 1] = 0.5*(diam[0] + diam[1])
w.diam[1, 1] = diam[1]

if w.ncomp == 3: # solvent spheres
    if len(diam) == 2: diam.append(0.5*(diam[0] + diam[1]))
    w.diam[0, 2] = 0.5*(diam[0] + diam[2])
    w.diam[1, 2] = 0.5*(diam[1] + diam[2])
    w.diam[2, 2] = diam[2]

# Over-write pairwise diameters for non-additivity

if len(diam) == 4: w.diam[0, 1] = diam[3]
if len(diam) == 5: w.diam[0, 2] = w.diam[1, 2] = diam[4]
if len(diam) == 6: w.diam[1, 2] = diam[5]

w.rpm_potential()

rhoz_init = eval(args.rhoz.replace('^', '**')) # total charged species density
rhos_init = eval(args.rhos.replace('^', '**')) # added solvent density

def solve(rhoz, rhos):
    """solve the structure at the given densities"""
    w.rho[0] = rhoz/2
    w.rho[1] = rhoz/2
    if w.ncomp == 3:
        w.rho[2] = rhos
    w.hnc_solve()
    if w.return_code: exit()

def selector(individual):
    """return a list according to  (hnn, hzz) and (h00, h01) representations"""
    if individual:
        return [[1, 0, 'g', 'h00'], [0, 1, 'b', 'h01']]
    else:
        return [[0.5, 0.5, 'k', '(h00+h01)/2'], [0.5, -0.5, 'r', '(h00-h01)/2']]

def update(val):
    """update state point from sliders, solve, and replot"""
    w.lb = 1 / tstar_slider.val
    w.rpm_potential()
    rhoz = 10**rhoz_slider.val
    rhos = 10**rhos_slider.val if rhos_slider else rhos_init
    solve(rhoz, rhos)
    replot()

def get_ann_txt():
    """get a string for annotating the plot"""
    rhoz = w.rho[0] + w.rho[1]
    if w.ncomp == 2:
        msg = 'T* = %5.2f  ρz = %8.3f  HNC err = %0.1g' % (1/w.lb, rhoz, w.error)
    else:
        rhos = w.rho[2]
        msg = 'T* = %5.2f  ρz = %8.3f  ρs = %8.3f  HNC err = %0.1g' % (1/w.lb, rhoz, rhos, w.error)
    return msg

def replot():
    """replot the lines and re-annotate"""
    for i, (p, q, color, text) in enumerate(selector(args.swap)):
        r = w.r[imin:imax]
        sgn = 1.0 if args.neg else -1.0 # double negative intentional
        rh = sgn * r * (p*w.hr[imin:imax, 0, 0] + q*w.hr[imin:imax, 0, 1])
        rh[rh < 0] = 1e-20
        line[i].set_ydata(np.log10(rh))
    ann.set_text(get_ann_txt())
    fig.canvas.draw_idle()

# Set up the plot area

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.30)

ax.set_xlim([1, args.rmax])
ax.set_ylim([-12, 1])
ax.set_xticks(list(range(0, int(5+args.rmax+0.5), 5)))
ax.set_yticks(list(range(-12, 2, 1)))
ax.set_xlabel('r / σ', labelpad=-10)
ax.set_ylabel('log10 r h(r)')

# report on diameters

txt = 'diams : %0.2f %0.2f' % (w.diam[0, 0], w.diam[1, 1])

if w.ncomp == 3:
    txt = txt + ' %0.2f' % w.diam[2, 2]

txt = txt + ' excess : %0.2f' % (w.diam[0, 1] - 0.5*(w.diam[0, 0] + w.diam[1, 1]))

if w.ncomp == 3:
    txt = txt + ' %0.2f' % (w.diam[0, 2] - 0.5*(w.diam[0, 0] + w.diam[2, 2]))
    txt = txt + ' %0.2f' % (w.diam[1, 2] - 0.5*(w.diam[1, 1] + w.diam[2, 2]))

ax.annotate(txt, xy=(0.02, 1.08), xycoords='axes fraction')

# solve and make initial plot

imin = int(1.0 / w.deltar)
imax = int(args.rmax / w.deltar)

line = [None, None] # will contain the data for the lines
label = [None, None] # will contain the labels for the lines

solve(rhoz_init, rhos_init)
ann = ax.annotate(get_ann_txt(), xy=(0.02, 1.02), xycoords='axes fraction')

style = 'dashed' if args.neg else 'solid'

for i, (p, q, color, text) in enumerate(selector(args.swap)):
    r = w.r[imin:imax]
    sgn = 1.0 if args.neg else -1.0 # double negative intentional
    rh = sgn * r * (p*w.hr[imin:imax, 0, 0] + q*w.hr[imin:imax, 0, 1])
    rh[rh < 0] = 1e-20
    line[i], = ax.plot(r, np.log10(rh), color+'-', linestyle=style)
    label[i] = ax.annotate(text, xy=(0.2+0.3*i, 0.92), color=color, xycoords='axes fraction')

# Set up sliders for lB, rho_z, and rho_s if required

ax_tstar = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor='PaleTurquoise')
tstar_slider = Slider(ax_tstar, 'T*', 0.2, 2.0, valinit=1/lb_init, valstep=0.01)
tstar_slider.on_changed(update)

ax_rhoz = plt.axes([0.25, 0.10, 0.65, 0.03], facecolor='PaleTurquoise')
rhoz_slider = Slider(ax_rhoz, 'ρ_z', -3, 0, valinit=m.log10(rhoz_init))
rhoz_slider.on_changed(update)

if w.ncomp == 3:
    ax_rhos = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor='PaleTurquoise')
    rhos_slider = Slider(ax_rhos, 'ρ_s', -3, 0, valinit=m.log10(rhos_init))
    rhos_slider.on_changed(update)
else:
    rhos_slider = None

# Set up buttons

def flip_sign(event):
    """invert signs of h"""
    args.neg = not args.neg
    style = 'dashed' if args.neg else 'solid'
    [ line[i].set_linestyle(style) for i in [0, 1] ]
    flip_button.label.set_text('-ve' if args.neg else '+ve')
    replot()

ax_flip = plt.axes([0.05, 0.25, 0.1, 0.03])
flip_button = Button(ax_flip, '-ve' if args.neg else '+ve', color='PaleTurquoise', hovercolor='0.975')
flip_button.on_clicked(flip_sign)

def swap(event):
    """swap between (hnn, hzz) and (h00, h01) representations"""
    args.swap = not args.swap
    for i, (p, q, color, text) in enumerate(selector(args.swap)):
        line[i].set_color(color)
        label[i].set_color(color)
        label[i].set_text(text)
    replot()

ax_swap = plt.axes([0.05, 0.20, 0.1, 0.03])
swap_button = Button(ax_swap, 'swap', color='PaleTurquoise', hovercolor='0.975')
swap_button.on_clicked(swap)

def reset(event):
    """reset all slider positions"""
    tstar_slider.reset()
    rhoz_slider.reset()
    if rhos_slider:
        rhos_slider.reset()

ax_reset = plt.axes([0.05, 0.15, 0.1, 0.03])
reset_button = Button(ax_reset, 'reset', color='PaleTurquoise', hovercolor='0.975')
reset_button.on_clicked(reset)

def dump(event):
    """write state point (T*, rho_z, kappa, [rho_s]) to std out"""
    rhoz = w.rho[0] + w.rho[1]
    kappa = m.sqrt(4*π*rhoz*w.lb)
    if w.ncomp == 2:
        print('%g\t%g\t%g' % (1/w.lb, rhoz, kappa))
    else:
        print('%g\t%g\t%g\t%g' % (1/w.lb, rhoz, kappa, w.rho[2]))

ax_dump = plt.axes([0.05, 0.10, 0.1, 0.03])
dump_button = Button(ax_dump, 'dump', color='PaleTurquoise', hovercolor='0.975')
dump_button.on_clicked(dump)

def quit(event):
    exit(0)

ax_quit = plt.axes([0.05, 0.05, 0.1, 0.03])
quit_button = Button(ax_quit, 'quit', color='PaleTurquoise', hovercolor='0.975')
quit_button.on_clicked(quit)

plt.show()

