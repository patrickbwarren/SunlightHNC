#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This file includes unicode characters like π = 3.14159

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Main script copyright (c) 2009-2019 Unilever UK Central Resources Ltd
# (Registered in England & Wales, Company No 29140; Registered
# Office: Unilever House, Blackfriars, London, EC4P 4BQ, UK).

# ZoomPan was adapted from
# https://stackoverflow.com/questions/11551049/matplotlib-plot-zooming-with-scroll-wheel
# copyright (c) remains with the original authors.

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

# By default this is for RPM without solvent
# Run with --solvated for solvent primitive model

# MOUSE AND KEYBOARD CONTROLS

# pan and zoom in/out with mouse wheel work in the main plot window
# control + mouse wheel zooms in/out horizontal axis

# sliders can be adjusted with mouse (move pointer onto slider)
# click to jump to a given position
# scroll with mouse wheel for medium adjustment
# shift + mouse wheel for fine adjustment
# control + mouse wheel for coarse adjustment
# 'v' with the pointer on a slider to request entry of a specific value (in terminal window)

# The following correspond to the Fig 1 insets in Coupette et al., PRL 121, 075501 (2018)
# For this lB = 0.7 nm, sigma = 0.3 nm, [salt] = 1 M, and [solvent] = 10 M and 40 M.
# Note the use of computed values for T* = sigma/lB and rho*sigma^3.
# Also 1 M = 0.602 molecules per nm^3
# For the salt, rhoz = [Na+] + [Cl-] = 2 [NaCl], hence the factor 2 in --rhoz

# python3 rpm_explorer.py --solvated --tstar=0.3/0.7 --rhoz=2*0.602*0.3^3 --rho=10*0.602*0.3^3
# python3 rpm_explorer.py --solvated --tstar=0.3/0.7 --rhoz=2*0.602*0.3^3 --rho=40*0.602*0.3^3

# Add --diam='[0.25/0.3,0.3373/0.3,1]' to reproduce the size-asymmetric model shown in Fig S1.

# Note that the densities are defined such that ρ+ = ρ- = ρz/2 , ρt = ρz + ρs (if solvated)

import argparse
import math as m
import numpy as np
import matplotlib.pyplot as plt

from numpy import pi as π
from oz import wizard as w
from matplotlib.widgets import Slider, Button, RadioButtons

# What might be shown in the plot:

selectors = [[['np.outer([0.5, 0.5, 0.0], [0.5, 0.5, 0.0])', 'k', '(h00+2h01+h11)/4'],
              ['np.outer([0.5, -0.5, 0.0], [0.5, -0.5, 0.0])', 'r', '(h00-2h01+h11)/4']],
             [['np.array([[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.0]])', 'g', '(h00+h11)/2'],
              ['np.array([[0.0, 0.5, 0.0], [0.5, 0.0, 0.0], [0.0, 0.0, 0.0]])', 'b', '(h01+h10)/2']],
             [['np.outer(x, x)', 'k', 'h_ρρ'],
              ['np.outer(z, z)', 'y', 'h_zz']]]

parser = argparse.ArgumentParser(description='Interactive RPM explorer')

parser.add_argument('--ng', action='store', default='16384', help='number of grid points (default 2^14 = 16384)')
parser.add_argument('--deltar', action='store', default=1e-3, type=float, help='grid spacing (default 1e-3)')

parser.add_argument('--alpha', action='store', default=0.2, type=float, help='Picard mixing fraction (default 0.2)')
parser.add_argument('--npic', action='store', default=6, type=int, help='number of Picard steps (default 6)')
parser.add_argument('--nps', action='store', default=6, type=int, help='length of history array (default 6)')
parser.add_argument('--maxsteps', action='store', default=100, type=int, help='number of iterations (default 100)')

parser.add_argument('--diam', action='store', default='1', help='hard core diameters (default 1)')
parser.add_argument('--sigma', action='store', default=0.0, type=float, help='inner core diameter (default min diam)')
parser.add_argument('--rhoz', action='store', default='0.1', help='total ion density (default 0.1)')
parser.add_argument('--rhos', action='store', default='0.4', help='added solvent density (default 0.4)')
parser.add_argument('--rhot', action='store', default=None, help='added solvent density (default computed)')
parser.add_argument('--tstar', action='store', default='1.0', help='reduced temperature (default 1.0)')
parser.add_argument('-s', '--solvated', action='store_true', help='for solvated primitive models')

parser.add_argument('--rmax', action='store', default=15.0, type=float, help='maximum radial distance (default 15)')
parser.add_argument('--floor', action='store', default=1e-20, type=float, help='floor for r h(r) (default 1e-20)')

parser.add_argument('--verbose', action='store_true', help='more output')

args = parser.parse_args()

args.show = 0
args.choice = ['both', 'both'] # which signs of h(r) to show

w.ncomp = 3 if args.solvated else 2
w.ng = eval(args.ng.replace('^', '**')) # catch 2^10 etc
w.deltar = args.deltar
w.alpha = args.alpha
w.npic = args.npic
w.nps = args.nps
w.maxsteps = args.maxsteps
w.verbose = args.verbose

w.initialise()

# if the user sets tstar to a string (eg 'infinity') this is caught here

try:
    tstar_init = eval(args.tstar)
except NameError:
    tstar_init = 0

w.lb = 1/tstar_init if tstar_init else 0.0

# Now construct the hard core diameters

diam = eval(f'[{args.diam}]') # wrap into a list

if len(diam) == 1: diam.append(diam[0])

w.diam[0, 0] = diam[0]
w.diam[0, 1] = 0.5*(diam[0] + diam[1])
w.diam[1, 1] = diam[1]

if args.solvated:
    if len(diam) == 2: diam.append(0.5*(diam[0] + diam[1]))
    w.diam[0, 2] = 0.5*(diam[0] + diam[2])
    w.diam[1, 2] = 0.5*(diam[1] + diam[2])
    w.diam[2, 2] = diam[2]

# Over-write pairwise diameters for non-additivity

if len(diam) == 4: w.diam[0, 1] = diam[3]
if len(diam) == 5: w.diam[0, 2] = w.diam[1, 2] = diam[4]
if len(diam) == 6: w.diam[1, 2] = diam[5]

w.sigma = args.sigma

w.rpm_potential()

def to_val(s):
    '''return a value from a string, replacing exponentiation symbol'''
    return eval(s.replace('^', '**'))

rhoz_init = to_val(args.rhoz) # total charged species density
rhos_init = to_val(f'{args.rhot} - {args.rhoz}') if args.rhot is not None else to_val(args.rhos) # solvent density

def solve(rhoz, rhos):
    """solve the structure at the given densities"""
    w.rho[0] = rhoz/2 # single point of truth where the densities are set for the HNC solver
    w.rho[1] = rhoz/2
    if args.solvated:
        w.rho[2] = rhos
    w.hnc_solve()
    if w.return_code: exit()

def update(val):
    """update state point from sliders, solve, and replot"""
    if tstar_slider:
        w.lb = 1 / tstar_slider.val
    w.sigma = args.sigma
    w.rpm_potential()
    rhoz = 10**rhoz_slider.val
    rhos = 10**rhos_slider.val if rhos_slider else rhos_init
    solve(rhoz, rhos)
    replot()

def update_both(val):
    '''call this if should update both sliders to match rhoz'''
#    rhoz, rhos = [10**slider.val for slider in [rhoz_slider, rhos_slider]] 
#    rhot = np.sum(w.rho)
#    new_rhos = rhot - rhoz
    rhos = np.sum(w.rho) - 10**rhoz_slider.val
    # print('rhoz, rhos, rhot, new_rhos, new_rhos+rhoz =', rhoz, rhos, rhot, new_rhos, new_rhos+rhoz)
    rhos_slider.set_val(m.log10(rhos if rhos > 0 else rhos_init))
    # print('slider vals =', [10**slider.val for slider in [rhoz_slider, rhos_slider]])
    update(val)

def get_ann_txt():
    """get a string for annotating the plot"""
    rhoz = w.rho[0] + w.rho[1]
    tstar = '%5.3f' % (1/w.lb) if w.lb else '∞'
    if args.solvated:
        rhos = w.rho[2]
        msg = 'T, ρ = %s, %8.4f + %8.4f = %8.4f  [err %0.1g]' % (tstar, rhos, rhoz, rhos+rhoz, w.error)
    else:
        msg = 'T, ρ = %s, %8.4f  [err %0.1g]' % (tstar, rhoz, w.error)
    return msg

def replot():
    """replot the lines and re-annotate"""
    for i, (wgts, color, text) in zip([0, 1], selectors[args.show]):
        x = w.rho / np.sum(w.rho) # mole fractions
        z = 0.5 * w.z # for charge-charge correlation
        wgt = eval(wgts) # this covers all cases
        # print(wgts) ; print(wgt)
        r = w.r[imin:imax]
        h = wgt[0, 0]*w.hr[imin:imax, 0, 0] + wgt[0, 1]*w.hr[imin:imax, 0, 1] \
            + wgt[1, 0]*w.hr[imin:imax, 1, 0] + wgt[1, 1]*w.hr[imin:imax, 1, 1]
        if args.solvated:
            h = h + wgt[0, 2]*w.hr[imin:imax, 0, 2] + wgt[1, 2]*w.hr[imin:imax, 1, 2] \
                + wgt[2, 2]*w.hr[imin:imax, 2, 2]
        rh_pos = r * h
        rh_neg = - rh_pos
        rh_pos[rh_pos < args.floor] = args.floor
        rh_neg[rh_neg < args.floor] = args.floor
        if line[2*i]: # if line[i] is not None then reset the y data.
            line[2*i].set_ydata(np.log10(rh_neg))
            line[2*i+1].set_ydata(np.log10(rh_pos))
        else: # plotting for the first time
            line[2*i], = ax.plot(r, np.log10(rh_neg), color+'-')
            line[2*i+1], = ax.plot(r, np.log10(rh_pos), color+'--')
            label[i] = ax.annotate(text, xy=(0.2+0.4*i, 0.92), color=color, xycoords='axes fraction')
        choice = args.choice[i]
        line[2*i].set_linestyle('None' if choice == 'none' or choice == '-ve' else 'solid')
        line[2*i+1].set_linestyle('None' if choice == 'none' or choice == '+ve' else 'dashed')
    ann.set_text(get_ann_txt())
    fig.canvas.draw_idle()

# Set up the plot area

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.30)

ax.set_xlim([1, args.rmax])
ax.set_ylim([-12, 1])
ax.set_xlabel('r / σ')
ax.set_ylabel('log10[- r h(r)]')

# report on diameters

txt = 'diams : %0.2f %0.2f' % (w.diam[0, 0], w.diam[1, 1])

if args.solvated:
    txt = txt + ' %0.2f' % w.diam[2, 2]

txt = txt + ' excess : %0.2f' % (w.diam[0, 1] - 0.5*(w.diam[0, 0] + w.diam[1, 1]))

if args.solvated:
    txt = txt + ' %0.2f' % (w.diam[0, 2] - 0.5*(w.diam[0, 0] + w.diam[2, 2]))
    txt = txt + ' %0.2f' % (w.diam[1, 2] - 0.5*(w.diam[1, 1] + w.diam[2, 2]))

ax.annotate(txt, xy=(0.02, 1.08), xycoords='axes fraction')

# solve and make initial plot

imin = int(1.0 / w.deltar)
imax = int(args.rmax / w.deltar)

line = [None, None, None, None] # will contain the data for the lines
label = [None, None] # will contain the labels for the lines

ann = ax.annotate(get_ann_txt(), xy=(0.02, 1.02), xycoords='axes fraction')

solve(rhoz_init, rhos_init)
replot()

# Set up sliders for lB, rho_z, and rho_s if required

back_color = 'powderblue'

if tstar_init > 0:
    ax_tstar = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor=back_color)
    tstar_slider = Slider(ax_tstar, 'T*', 0.0, max(2.0, tstar_init), valinit=tstar_init, valstep=0.01, valfmt='%5.3f')
    tstar_slider.on_changed(update)
else:
    tstar_slider = None

ax_rhoz = plt.axes([0.25, 0.10 if tstar_slider else 0.05, 0.65, 0.03], facecolor=back_color)
rhoz_slider = Slider(ax_rhoz, 'ρ_z', -3, 0, valinit=m.log10(rhoz_init), valfmt='%5.3f')
rhoz_slider.on_changed(update_both if args.solvated and args.rhot is not None else update)

if args.solvated:
    ax_rho = plt.axes([0.25, 0.15 if tstar_slider else 0.10, 0.65, 0.03], facecolor=back_color)
    rhos_slider = Slider(ax_rho, 'ρ_t' if args.solvated else 'ρ_s', -3, 0, valinit=m.log10(rhos_init), valfmt='%5.3f')
    rhos_slider.on_changed(update)
else:
    rhos_slider = None

# Set up buttons

def radio1(val):
    """Select between showing both, +ve, -ve or none for first h(r)"""
    args.choice[0] = val
    replot()

def radio2(val):
    """Select between showing both, +ve, -ve or none for second h(r)"""
    args.choice[1] = val
    replot()

radio = [radio1, radio2]

ax_choice = [None, None]
ax_choice[0] = plt.axes([0.05, 0.42, 0.1, 0.15])
ax_choice[1] = plt.axes([0.05, 0.25, 0.1, 0.15])

choice = [None, None]

for i in [0, 1]:
    choice[i] = RadioButtons(ax_choice[i], ('none', '+ve', '-ve', 'both'), active=3)
    choice[i].on_clicked(radio[i])

for i, (wgts, color, text) in zip([0, 1], selectors[args.show]):
    [ label.set_color(color) for label in choice[i].labels ]

def advance(event):
    """advance between (hnn, hzz) and (h00, h01) representations"""
    args.show = (args.show + 1) % len(selectors) # advance through the selections
    for i, (wgts, color, text) in zip([0, 1], selectors[args.show]):
        [ line[2*i+j].set_color(color) for j in [0, 1] ]
        label[i].set_color(color)
        label[i].set_text(text)
        [ label.set_color(color) for label in choice[i].labels ]
    replot()

ax_advance = plt.axes([0.05, 0.20, 0.1, 0.03])
advance_button = Button(ax_advance, 'cycle', color=back_color,  hovercolor='0.975')
advance_button.on_clicked(advance)

def reset(event):
    """reset all slider positions and plot area"""
    if tstar_slider:
        tstar_slider.reset()
    rhoz_slider.reset()
    if rhos_slider:
        rhos_slider.reset()
    ax.set_xlim([1, args.rmax])
    ax.set_ylim([-12, 1])
    ax.figure.canvas.draw()

ax_reset = plt.axes([0.05, 0.15, 0.1, 0.03])
reset_button = Button(ax_reset, 'reset', color=back_color, hovercolor='0.975')
reset_button.on_clicked(reset)

def dump(event):
    """write state point (T*, rho_z, kappa, [rho_s]) to std out"""
    rhoz = w.rho[0] + w.rho[1]
    kappa = m.sqrt(4*π*w.lb*rhoz)
    rhos = w.rho[2] if args.solvated else 0
    tstar = '%g' % (1/w.lb) if w.lb else '∞'
    print('%s\t%g\t%g\t%g\t%g' % (tstar, rhos, rhoz, rhoz+rhoz, kappa))

ax_dump = plt.axes([0.05, 0.10, 0.1, 0.03])
dump_button = Button(ax_dump, 'dump', color=back_color, hovercolor='0.975')
dump_button.on_clicked(dump)

def quit(event):
    exit(0)

ax_quit = plt.axes([0.05, 0.05, 0.1, 0.03])
quit_button = Button(ax_quit, 'quit', color=back_color, hovercolor='0.975')
quit_button.on_clicked(quit)

# ZoomPan was adapted from
# https://stackoverflow.com/questions/11551049/matplotlib-plot-zooming-with-scroll-wheel
# copyright (c) remains with the original authors.

class ZoomPan:

    def __init__(self, ax):
        self.ax = ax
        self.press = False
        self.cur_xlim = None
        self.cur_ylim = None
        self.xpress = None
        self.ypress = None
        self.control_down = False

    def factory(self, base_scale=1.1):

        def zoom(event):
            if event.inaxes != self.ax: return
            cur_xlim = self.ax.get_xlim()
            cur_ylim = self.ax.get_ylim()
            xdata, ydata = event.xdata, event.ydata
            if event.button == 'down':
                scale_factor = 1 / base_scale
            elif event.button == 'up':
                scale_factor = base_scale
            new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
            if self.control_down:
                new_height = (cur_ylim[1] - cur_ylim[0])
            else:
                new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor
            relx = (cur_xlim[1] - xdata)/(cur_xlim[1] - cur_xlim[0])
            rely = (cur_ylim[1] - ydata)/(cur_ylim[1] - cur_ylim[0])
            self.ax.set_xlim([xdata - new_width * (1-relx), xdata + new_width * (relx)])
            self.ax.set_ylim([ydata - new_height * (1-rely), ydata + new_height * (rely)])
            self.ax.figure.canvas.draw()

        def button_down(event):
            if event.inaxes != self.ax: return
            self.cur_xlim = self.ax.get_xlim()
            self.cur_ylim = self.ax.get_ylim()
            self.press = True
            self.xpress, self.ypress = event.xdata, event.ydata

        def button_up(event):
            self.press = False
            self.ax.figure.canvas.draw()

        def mouse_move(event):
            if self.press is False: return
            if event.inaxes != self.ax: return
            dx = event.xdata - self.xpress
            dy = event.ydata - self.ypress
            self.cur_xlim -= dx
            self.cur_ylim -= dy
            self.ax.set_xlim(self.cur_xlim)
            self.ax.set_ylim(self.cur_ylim)
            self.ax.figure.canvas.draw()

        def key_down(event):
            if event.key == 'control':
                self.control_down = True

        def key_up(event):
            if event.key == 'control':
                self.control_down = False

        fig = self.ax.get_figure()
        fig.canvas.mpl_connect('button_press_event', button_down)
        fig.canvas.mpl_connect('button_release_event', button_up)
        fig.canvas.mpl_connect('motion_notify_event', mouse_move)
        fig.canvas.mpl_connect('key_press_event', key_down)
        fig.canvas.mpl_connect('key_release_event', key_up)
        fig.canvas.mpl_connect('scroll_event', zoom)

        return zoom, button_down, button_up, mouse_move

zp = ZoomPan(ax).factory()

class SliderScroll:

    def __init__(self, ax):
        self.ax = ax
        self.shift_down = False
        self.control_down = False

    def factory(self, sliders, superfine_scale=0.001, fine_scale=0.01, coarse_scale=0.1):

        def scroll(event):
            if event.inaxes in sliders:
                slider = sliders[event.inaxes]
                saltus = superfine_scale if self.shift_down else coarse_scale if self.control_down else fine_scale
                val = slider.val
                if event.button == 'down':
                    val = val - saltus
                elif event.button == 'up':
                    val = val + saltus
                slider.set_val(val)
                update(val)

        def key_down(event):
            if event.key == 'shift':
                self.shift_down = True
            elif event.key == 'control':
                self.control_down = True

        def key_up(event):
            if event.key == 'shift':
                self.shift_down = False
            elif event.key == 'control':
                self.control_down = False
            elif event.key == 'v':
                if event.inaxes in sliders:
                    slider = sliders[event.inaxes]
                    s = input('enter a value for %s\n' % slider_name[slider])
                    try:
                        val = float(s)
                        print('%s set to %g' % (slider_name[slider], val))
                        slider.set_val(m.log10(val) if log_slider[slider] else val)
                        update(val)
                    except ValueError:
                        print('invalid number ', s)

        fig = self.ax.get_figure()
        fig.canvas.mpl_connect('scroll_event', scroll)
        fig.canvas.mpl_connect('key_press_event', key_down)
        fig.canvas.mpl_connect('key_release_event', key_up)

        return scroll, key_up, key_down

sliders = {ax_rhoz: rhoz_slider}
slider_name = {rhoz_slider: 'rho_z'}
log_slider = {rhoz_slider: True}

if tstar_slider:
    sliders[ax_tstar] = tstar_slider
    slider_name[tstar_slider] = 'T*'
    log_slider[tstar_slider] = False

if rhos_slider:
    sliders[ax_rho] = rhos_slider
    slider_name[rhos_slider] = 'rho_t' if args.rhot is not None else 'rho_s'
    log_slider[rhos_slider] = True

ss = SliderScroll(ax).factory(sliders)

plt.show()
