#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Wrapper for functionality in oz_mod.f90

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.  The main intent is to hide the 'wizard' interface.

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
# modifications copyright (c) 2009-2017 Unilever UK Central Resources
# Ltd (Registered in England & Wales, Company No 29140; Registered
# Office: Unilever House, Blackfriars, London, EC4P 4BQ, UK).  Later
# modifications copyright (c) 2020-2024 Patrick B Warren
# <patrick.warren@stfc.ac.uk> and STFC.

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

import numpy as np
from copy import copy
from numpy import pi as π
from numpy import sin, sqrt, exp
from scipy.special import erf, erfc

def F90str(s): # convert FORTRAN string to python
    return str(s, 'utf-8').strip()

class Grid: # instantiate one of these first

    def __init__(self, wizard, ncomp=1, ng=16384, deltar=0.01): # instantiate and initialise
        wizard.ng = ng
        wizard.ncomp = ncomp
        wizard.deltar = deltar
        wizard.initialise() # all the arrays in the FORTRAN sector
        self.ng = wizard.ng
        self.ncomp = wizard.ncomp
        self.r = wizard.r # copy these for use outside
        self.k = wizard.k
        self.deltar = wizard.deltar
        self.deltak = wizard.deltak
        self.version = F90str(wizard.version)
        self.wizard = wizard # give access to underlying FORTRAN sector

class Model: # build a model on top of the grid

    def __init__(self, grid, name):
        self.name = name # give it a name, other properties may be attached later
        self.wizard = grid.wizard # give access to underlying FORTRAN sector

    def write_params(self):
        self.wizard.write_params()

class Solution: # intended for internal use only; assumes direct access to wizard

    def __init__(self, wizard): # copy some things that we might need (add more here !)
        self.error = wizard.error
        self.deficit = wizard.deficit
        self.aex = wizard.aex
        self.press = wizard.press
        self.uex = wizard.uex
        self.muex = wizard.muex
        self.hc = wizard.hc
        self.hr = wizard.hr
        self.closure = F90str(wizard.closure_name)
        self.nsteps = wizard.nsteps
        self.wizard = wizard # give access to underlying FORTRAN sector

    def write_params(self):
        self.wizard.write_params()

    def write_thermodynamics(self):
        self.wizard.write_thermodynamics()

def additive_primitive_model(grid, lb, diam, z):
    model = Model(grid, 'additive PM')
    model.lb = copy(lb)
    model.z = np.copy(z)
    model.diam = np.copy(diam)
    w = model.wizard
    r, k, ncomp = w.r, w.k, w.ncomp
    model.dd = np.empty((ncomp, ncomp)) # a new model property
    for j in range(ncomp):
        for i in range(j+1):
            ij = i + j*(j+1)//2
            model.dd[i, j] = 0.5*(diam[i] + diam[j])
            w.dd[ij] = model.dd[i, j]
    σ = np.min(w.dd) # this defines the minimum hard core
    model.sigma = σ # stash this as a property of the model
    cut = round(σ / w.deltar)
    for j in range(ncomp):
        for i in range(j+1):
            ij = i + j*(j+1)//2
            zzlb = z[i] * z[j] * lb
            w.ulong[:, ij] = zzlb / r
            w.ulong[:cut, ij] = zzlb / σ # cut-off inside min hard core (see docs)
            w.dulong[:, ij] = - zzlb / r**2
            w.dulong[:cut, ij] = 0.0 # -- ditto --
            w.ulongk[:, ij] = 4*π*zzlb*sin(k*σ) / (σ * k**3)
    w.ushort[:, :] = 0.0
    w.dushort[:, :] = 0.0
    w.expnegus[:, :] = 1.0
    for ij in range(w.nfnc):
        cut = round(w.dd[ij] / w.deltar)
        w.expnegus[:cut, ij] = 0.0
    w.u0[:] = z**2 * lb / σ
    w.tp[:] = 0.0
    w.tu[:] = 0.0
    w.tl[:] = 0.0
    return model

def restricted_primitive_model(grid, lb):
    diam = np.array([1.0, 1.0])
    z = np.array([1.0, -1.0])
    model = additive_primitive_model(grid, lb, diam, z)
    model.name = 'RPM'
    return model

def soften_rpm(model, kappa, ushort=False): # to be called after RPM
    lb = model.lb
    σ = model.sigma
    w = model.wizard
    model.kappa = copy(kappa)
    r, k, κ = w.r, w.k, kappa
    if ushort: # modify Ushort, keeping Ulong unchanged
        w.ushort[:, 1] = lb*erfc(κ*r) / r # note this flattens as r --> 0, unlike the bare Coulomb law
        w.dushort[:, 1] = - lb*erfc(κ*r) / r**2 - 2*κ*lb*exp(-κ**2*r**2) / (sqrt(π)*r)
        cut = round(w.dd[1] / w.deltar)
        w.expnegus[cut:, 1] = exp(-w.ushort[cut:, 1])
    else: # modify Ulong, keeping Ushort unchanged
        w.ulong[:, 1]  = - lb*erf(κ*r) / r
        w.dulong[:, 1] = lb*erf(κ*r) / r**2 - 2*κ*lb*exp(-κ**2*r**2) / (sqrt(π)*r)
        w.ulongk[:, 1] = - 4*π*lb*exp(-k**2/(4*κ**2)) / k**2
        w.tl[1] = π*lb/κ**2 - 2/3*π*lb*σ**2 # off SYM condition (see docs)
    w.tp[1] = π*lb * ( σ*exp(-κ**2*σ**2)/(κ*sqrt(π)) + (1/(2*κ**2) - 1/3*σ**2) * erfc(κ*σ) )
    w.tu[1] = π*lb * ( σ*exp(-κ**2*σ**2)/(κ*sqrt(π)) + (1/(2*κ**2) - σ**2) * erfc(κ*σ) )
    what = 'Ushort' if ushort else 'Ulong'
    model.name = f'softened RPM ({what} changed)'
    return model

def solve(model, rho, closure='HNC', npic=6, alpha=0.2, maxsteps=100, cold_start=None):
    model.wizard.rho = rho
    model.wizard.npic = npic
    model.wizard.alpha = alpha
    model.wizard.maxsteps = maxsteps
    if cold_start is not None: # avoid setting this if not required
        model.wizard.cold_start = cold_start
    if 'hnc' in closure.lower():
        model.wizard.hnc_solve()
    elif 'msa' in closure.lower():
        model.wizard.msa_solve()
    else:
        raise ValueError(f'unrecognised closure {closure} in solve; use HNC or MSA')
    return Solution(model.wizard) # return a solution object

def hnc_solve(*args, **kwargs): # for backwards compatibility ; deprecated
    return solve(*args, closure='HNC', **kwargs)
