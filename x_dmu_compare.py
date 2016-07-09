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

# A plot of the chemical potential difference against x reveals nearly
# linear behaviour, which can be confirmed in Monte-Carlo simulations
# (unreported).  This linear dependence can be represented by Delta F
# = chi (1 - 2x) where chi is the Flory chi-parameter and is the basic
# reason why DPD is so good at representing fluid mixtures which fit
# Flory-Huggins (regular solution) theory.

# The following table contains Monte-Carlo calculations of the
# chemical potential difference, using trial particle identity swaps.
# Extract the data by:

#  gawk '/##/ && NF==4 { print $2, $3, $4 }' x_dmu_compare.py > temp.dat

#    x     d(mu)     std-error

##  0.0   1.70174    0.00407275
##  0.1   1.28057    0.00260028
##  0.2   0.914987   0.00333303
##  0.3   0.597154   0.00386779
##  0.4   0.281795   0.00428741
##  0.5  -0.00181835 0.00325142
##  0.6  -0.290226   0.00341724
##  0.7  -0.595639   0.00328827
##  0.8  -0.919721   0.00357011
##  0.9  -1.28446    0.00297603
##  1.0  -1.70174    0.00407275

from oz import wizard as w

w.ncomp = 2
w.initialise()

rho = 3.0
A = 25.0
dA = 5.0

w.arep[0,0] = w.arep[1,1] = A
w.arep[0,1] = A + dA
w.dpd_potential(1)

n = 41

for i in range(n):

    x = i / (n - 1.0)

    w.rho[0] = (1.0 - x) * rho
    w.rho[1] = x * rho
    w.hnc_solve()

    print("%f\t%f\t%g" % (x, w.muex[1]-w.muex[0], w.error))

