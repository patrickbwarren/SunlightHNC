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

program driver2
  use wizard
  implicit none
  integer :: j
  real(kind=dp) :: rhotot
  real(kind=dp), allocatable :: ccsf(:), ddsf(:)

  verbose = .true.

  ng = 4096
  ncomp = 2

  call initialise

  lb = 120.0_dp
  arep = 25.0_dp
  arep(1, 2) = 30.0_dp
  z(1) = 1.0_dp
  z(2) = -1.0_dp

  call dpd_potential

  rhotot = 3.0_dp

  rho(1) = 0.5_dp * rhotot
  rho(2) = rho(1)

  call write_params

  call hnc_solve

  if (error .gt. 1.0E-10_dp) &
       & print *, 'Warning, did not converge to 1e-10'

  allocate(ccsf(ng-1))
  allocate(ddsf(ng-1))

  ddsf = sum(sum(sk, dim=3), dim=2) / sum(rho)
  ccsf = sk(:,1,1) - sk(:,1,2) - sk(:,2,1) + sk(:,2,2)

  open (unit=10, file='driver2_sf.dat')
  do j = 1, ng-1
     write (10, '(3(f15.8,2x))') k(j), ddsf(j), ccsf(j)
  end do
  close (10)
  print *, 'structure factors written to driver2_sf.dat'

  open (unit=11, file='driver2_hij.dat')
  do j = 1, ng-1
     write (11, '(4(f15.8,2x))') r(j), hr(j,1,1), hr(j,1,2), &
          & hr(j,2,2)
  end do
  close (11)
  print *, 'pair correlation functions written to driver2_gij.dat'

  call write_thermodynamics

end program driver2
