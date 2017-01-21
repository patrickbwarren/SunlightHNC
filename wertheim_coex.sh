#!/bin/bash

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

######################################################################

# Solution method for phase coexistence using URPM for right hand
# (high density) side and Wertheim with softened URPM for left hand
# (low density) side.  First off: note HNC solution exists at high
# densities but not low densities, this is why use softened URPM
# solved by HNC as reference.  Second note that the left hand phase
# boundary is at a low density, so the coexistence pressure is nearly
# zero.  Hence the zero pressure isochor of the HNC forms a very good
# first approximation to the right hand (high density) phase boundary.
# This can be used to determine a first guess at the coexistence
# chemical potential, which can be used to fix the low density phase
# boundary, as in the following scheme:

# (1) Choose a value for lB. (2) Calculate the HNC density and
# chemical potential for the URPM for which the pressure vanishes -
# this will form the first estimate of the high density phase
# boundary.  (3) Determine the density of the Wertheim + softened URPM
# fluid (solved by HNC) which has the same chemical potential - this
# can be done by hand, by interval halving, and forms the first
# estimate of the low density phase boundary.  (4) Record the pressure
# of the low density phase (which will be small and positive) and
# recalculate the high density phase boundary with this target
# pressure.  Iterate to self consistency.  In practice one iteration
# of this scheme seems sufficient.

# Note for reference the expected phase coexistence at lB = 120
# between rhoz = 0.015 and 0.05, from MC simulations.

# The following give some good estimates of the coexistence solutions,
# at various values of lB.  Smaller values of lB may not have a
# solution.  The initial guess for sigmap may need fine tuning as lB
# varies.  

# The results are :

# lB    density         pressure        chem_pot        sigmap
# 50	2.600000e-04	0.000172844	-17.3088	3.86101
# 50	0.021681	0.000172844	-17.2971
# 80	1.920000e-04	0.000135821	-25.2931	4.49317
# 80	0.0518256	0.000135821	-25.2926
# 100	1.590000e-04	0.000113434	-30.7495	4.87625
# 100	0.0689763	0.000114403	-30.7538
# 120	1.365000e-04	9.45216e-05	-36.2619	5.18208
# 120	0.0843975	9.45e-05	-36.2613
# 150	1.170000e-04	7.38239e-05	-44.5716	5.4927
# 150	0.105023	7.38239e-05	-44.5724

# The rows where sigmap is reported correspond to state points on the
# low density phase boundary (Wertheim + softened URPM, HNC).  The
# rows where sigmap is not reported correspond to state points on the
# high density phase boundary (URPM, HNC).  Obviously the aim is that
# the two state points should have the same pressure and chemical
# potential (to within the target numerical accuracy).

/usr/bin/env python3 wertheim_solver.py --lb=50 --rhoz=2.6e-4 --sigmap=4.5 --min
/usr/bin/env python3 urpm_targp.py --lb=50 --targp=0.000172844

/usr/bin/env python3 wertheim_solver.py --lb=80 --rhoz=1.92e-4 --sigmap=5 --min
/usr/bin/env python3 urpm_targp.py --lb=80 --targp=0.000135821

/usr/bin/env python3 wertheim_solver.py --lb=100 --rhoz=1.59e-4 --sigmap=5 --min
/usr/bin/env python3 urpm_targp.py --lb=100 --targp=0.000114403

/usr/bin/env python3 wertheim_solver.py --lb=120 --rhoz=1.365e-4 --sigmap=5 --min
/usr/bin/env python3 urpm_targp.py --lb=120 --targp=9.45e-5

/usr/bin/env python3 wertheim_solver.py --lb=150 --rhoz=1.17e-4 --sigmap=5 --min
/usr/bin/env python3 urpm_targp.py --lb=150 --targp=7.38239e-05
