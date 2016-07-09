#!/usr/bin/python3

# This file is part of SunlightDPD - a home for open source software
# related to the dissipative particle dynamics (DPD) simulation
# method.

# Based on an original code copyright (c) 2007 Lucian Anton.
# Modifications copyright (c) 2008, 2009 Andrey Vlasov.  
# Additional modifications copyright (c) 2009-2016 Unilever UK Central
# Resources Ltd (Registered in England & Wales, Company No 29140;
# Registered Office: Unilever House, Blackfriars, London, EC4P 4BQ,
# UK).

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

# The results of this calculation can be directly compared with Fig 4
# of the Groot and Warren [J. Chem. Phys. v107, 4423 (1997)].  Here is
# the data from that figure.  Extract this by:

#  gawk '/##/ && NF==3 { print $2, $3 }' gw_p_compare.py > temp.dat

#     rho (p-rho)/(A*rho^2)

##    0.0  0.0379935086163
##    1.5  0.0751786298043
##    2.5  0.0886823425022
##    3.0  0.0924251622846
##    3.5  0.0946639891655
##    4.0  0.0965259421847
##    5.0  0.0987451548125
##    6.0  0.0998358473824
##    7.0  0.100551067109
##    8.0  0.102017933031

from oz import wizard as w

w.initialise()
w.arep[0,0] = A = 25.0
w.dpd_potential(1)

npt = 41
rhomax = 10.0

for i in range(npt):
    w.rho[0] = rho = rhomax * (i + 1.0) / npt
    w.hnc_solve()
    print("%f\t%g\t%g" % (rho, (w.press-rho)/(A*rho*rho), w.error))

