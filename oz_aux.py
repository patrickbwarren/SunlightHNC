#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Wrapper for functionality in oz_mod.f90

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
# modifications copyright (c) 2009-2017 Unilever UK Central Resources
# Ltd (Registered in England & Wales, Company No 29140; Registered
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

import numpy as np
from numpy import pi as π
from numpy import sin, sqrt, exp
from scipy.special import erf, erfc

def F90str(s): # convert FORTRAN string to python
    return str(s, 'utf-8').strip()

class Solution:

    def __init__(self, wizard): # copy everything over that we might need
        self.error = wizard.error
        self.deficit = wizard.deficit
        self.aex = wizard.aex
        self.press = wizard.press
        self.uex = wizard.uex
        self.muex = wizard.muex
        self.hc = wizard.hc
        self.closure = F90str(wizard.closure_name)

    def write_params(self):
        self.wizard.write_params()

    def write_thermodynamics():
        self.wizard.write_thermodynamics()

class OZSolver:

    def __init__(self, wizard, ncomp=1, ng=16384, deltar=0.01): # instantiate and initialise
        wizard.ncomp = ncomp
        wizard.ng = ng
        wizard.deltar = deltar
        wizard.initialise()
        self.version = F90str(wizard.version)
        self.wizard = wizard # gives access from the instantiated object

    def additive_primitive_model(self, lb, diam, z):
        w = self.wizard
        r, k = w.r, w.k
        for j in range(w.ncomp):
            for i in range(j+1):
                ij = i + j*(j+1)//2
                w.dd[ij] = 0.5*(diam[i] + diam[j])
        sigma = np.min(w.dd) # this defines the minimum hard core
        cut = round(sigma / w.deltar)
        for j in range(w.ncomp):
            for i in range(j+1):
                ij = i + j*(j+1)//2
                zzlb = z[i] * z[j] * lb
                w.ulong[:, ij] = zzlb / r
                w.ulong[:cut, ij] = zzlb / sigma # cut-off inside min hard core (see docs)
                w.dulong[:, ij] = - zzlb / r**2
                w.dulong[:cut, ij] = 0.0 # -- ditto --
                w.ulongk[:, ij] = 4*π*zzlb*sin(k*sigma) / (sigma * k**3)
        w.ushort[:, :] = 0.0
        w.dushort[:, :] = 0.0
        w.expnegus[:, :] = 1.0
        for ij in range(w.nfnc):
            cut = round(w.dd[ij] / w.deltar)
            w.expnegus[:cut, ij] = 0.0
        w.u0[:] = z**2 * lb / sigma
        w.tp[:] = 0.0
        w.tu[:] = 0.0
        w.tl[:] = 0.0
        self.model = 'additive PM'

    def restricted_primitive_model(self, lb):
        diam = np.array([1.0, 1.0])
        z = np.array([1.0, -1.0])
        self.additive_primitive_model(lb, diam, z)
        self.model = 'RPM'

    def soften_rpm(self, lb, kappa, ushort=False): # to be called after RPM
        w = self.wizard
        σ = np.min(w.dd)
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
        self.model = f'softened RPM ({what} changed)'

    def set_verbosity(self, verbosity):
        self.verbosity = verbosity
        self.wizard.verbosity = (verbosity > 1)

    def hnc_solve(self, rho, alpha=0.3, npic=6, cold_start=False): # solve and return a solution object
        self.wizard.rho = rho
        self.wizard.npic = npic
        self.wizard.alpha = alpha
        self.wizard.cold_start = cold_start
        self.wizard.hnc_solve()
        return Solution(self.wizard)
