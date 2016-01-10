! This file is part of SunlightDPD - a home for open source software
! related to the dissipative particle dynamics (DPD) simulation
! method.

! Copyright (c) 2009-2016 Unilever UK Central Resources Ltd
! (Registered in England & Wales, Company No 29140; Registered Office:
! Unilever House, Blackfriars, London, EC4P 4BQ, UK).

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

! Test code to comprehend how FFTW works.

program fftw_test
implicit none
include "fftw3.f"
integer, parameter :: n = 47
integer :: i, j, k
integer*8 :: planA, planB
double precision :: a(n), b(n), c(n), d(n)
double precision, parameter :: pi=3.141592653589793d0

do i=1,n-1
   a(i) = 1.0d0 / sqrt(dble(i))
end do

call dfftw_plan_r2r_1d(planA, n-1, a, b, FFTW_RODFT00, FFTW_ESTIMATE)
call dfftw_plan_r2r_1d(planB, n-1, b, c, FFTW_RODFT00, FFTW_ESTIMATE)

call dfftw_execute(planA)
b = b / sqrt(2.0d0*n)

call dfftw_execute(planB)
c = c / sqrt(2.0d0*n)

do k = 1, n-1
   d(k) = 0.0d0
   do j = 1, n-1
      d(k) = d(k) + a(j) * sin(j*k*pi/n)
   enddo
   d(k) = d(k) * sqrt(2.0d0/n)
enddo

print *, "n = ", n

print *, "fftw_rodft00 = ", fftw_rodft00
print *, "fftw_estimate = ",fftw_estimate
 
print *, "address of planA = ", planA
print *, "address of planB = ", planB

print *, '  a       b         c      d'
print *, '  A  FFT(FFT(A))  FFT(A)  FT(A)'  

do i = 1, n-1
   print *, i, a(i), c(i), b(i), d(i)
end do
 
print *, "||a|| = ", sqrt(sum(a(:)**2))
print *, "||b|| = ", sqrt(sum(b(:)**2))
print *, "||c|| = ", sqrt(sum(c(:)**2))
print *, "||d|| = ", sqrt(sum(d(:)**2))

print *, "||a-c|| = ", sqrt(sum((a(:)-c(:))**2))
print *, "||b-d|| = ", sqrt(sum((b(:)-d(:))**2))

call dfftw_destroy_plan(planA)
call dfftw_destroy_plan(planB)

end

