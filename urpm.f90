! This file is part of SunlightDPD - a home for open source software
! related to the dissipative particle dynamics (DPD) simulation
! method.

! Based on an original code copyright (c) 2007 Lucian Anton.
! Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
! modifications copyright (c) 2009-2013 Unilever UK Central Resources
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

program driver
  use wizard
  implicit none
  integer :: j
  double precision :: rhotot

  verbose = 1

  ng = 4096
  ncomp = 2

  call initialise

  lb = 20
  sigma = 0.5
  arep = 25.0
  z(1) = 1
  z(2) = -1

  call dpd_potential(1)

  rhotot = 3.0d0

  rho(1) = 0.5d0 * rhotot
  rho(2) = rho(1)

  call write_params

  call hnc_solve

  if (error .gt. 1.0d-10) &
       & print *, 'Warning, did not converge to 1e-10'

  open (unit=11, file='urpm_hij.dat')
  do j = 1, ng-1
     write (11, '(4(f15.8,2x))') r(j), hr(j,1,1), hr(j,1,2), hr(j,2,2)
  end do
  close (11)
  print *, 'pair correlation functions written to urpm_hij.dat'

  call write_thermodynamics

end program driver
