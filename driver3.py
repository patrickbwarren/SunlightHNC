#!/usr/bin/env python3

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
# modifications copyright (c) 2009-2016 Unilever UK Central Resources
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
import matplotlib.pyplot as plt
from oz import wizard as w

w.ng = 4096
w.ncomp = 3

w.initialise()

w.lb = 200.0
w.arep[:,:] = 25.0
w.arep[0, 1] = 30.0
w.arep[0, 2] = 35.0
w.arep[1, 2] = 20.0
w.z[0] = 1
w.z[1] = -1

w.dpd_potential()

rho = 3.0
xc = 0.2

w.rho[0] = 0.5 * rho * xc
w.rho[1] = 0.5 * rho * xc
w.rho[2] = rho * (1 - xc)

w.write_params()

w.verbose = 1
w.hnc_solve()
w.write_thermodynamics()

# density-density structure factor

ddsf = np.sum(np.sum(w.sk, axis=2), axis=1) / np.sum(w.rho)

# charge-charge structure factor (notice how elegant this is :-)

ccsf = np.dot(np.dot(w.z, w.sk), w.z)

plt.figure(1)
imax = int(5.0 / w.deltar)
plt.plot(w.r[:imax], 1.0 + w.hr[:imax, 0, 0], label='$g_{11}$')
plt.plot(w.r[:imax], 1.0 + w.hr[:imax, 0, 1], label='$g_{12} = g_{21}$')
plt.plot(w.r[:imax], 1.0 + w.hr[:imax, 1, 1], label='$g_{22}$')
plt.plot(w.r[:imax], 1.0 + w.hr[:imax, 1, 2], label='$g_{13} = g_{31}$')
plt.plot(w.r[:imax], 1.0 + w.hr[:imax, 2, 2], label='$g_{33}$')
plt.legend(loc='lower right')
plt.xlabel('$r$')
plt.ylabel('$g(r)$')

plt.figure(2)
jmax = int(25.0 / w.deltak)
plt.plot(w.k[:jmax], ddsf[:jmax], label='$S_{NN}$')
plt.plot(w.k[:jmax], ccsf[:jmax], label='$S_{ZZ}$')
plt.legend(loc='lower right')
plt.xlabel('$k$')
plt.ylabel('$S(k)$')

plt.show()

