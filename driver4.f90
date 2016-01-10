! This file is part of SunlightDPD - a home for open source software
! related to the dissipative particle dynamics (DPD) simulation
! method.

! Based on an original code copyright (c) 2007 Lucian Anton.
! Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
! modifications copyright (c) 2009-2016 Unilever UK Central Resources
! Ltd (Registered in England & Wales, Company No 29140; Registered
! Office: Unilever House, Blackfriars, London, EC4P 4BQ, UK).

! SunlightDPD is free software: you can redistribute it and/or
! modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 3 of the
! License, or (at your option) any later version.

! SunlightDPD is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with SunlightDPD.  If not, see <http://www.gnu.org/licenses/>.

program driver4
  use wizard
  implicit none
  double precision :: rhotot, mfcharge

  verbose = 1

  ng = 4096
  ncomp = 4

  call initialise

  lb = 200.0
  arep = 25.0

  arep(1, 2) = 32.0
  arep(1, 3) = 22.0
  arep(1, 4) = 15.0
  arep(2, 3) = 30.0
  arep(2, 4) = 35.0
  arep(3, 4) = 20.0
  
  z(2) = 1
  z(3) = -1
  z(4) = -1

  call dpd_potential(1)

  rhotot = 3.0d0
  mfcharge = 0.2d0

  rho(1) = rhotot * (1.0d0 - mfcharge)
  rho(2) = 0.5d0 * rhotot * mfcharge
  rho(3) = rho(2) * 0.2
  rho(4) = rho(2) - rho(3)

  call write_params

  call hnc_solve

  if (error .gt. 1.0d-10) &
       & print *, 'Warning, did not converge to 1e-10'

  call write_thermodynamics

end program driver4
