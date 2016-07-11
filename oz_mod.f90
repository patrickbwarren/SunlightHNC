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

! ===================================================================

! The main purpose of this module is to provide integral equation
! closures for the multicomponent Ornstein-Zernike equations.  The
! closures currently provided are RPA, EXP and HNC.

! Better documentation (which reproduces the below) is found in the
! accompanying LaTeX document.

! Note on Fourier transforms, based on "A fast solver for the
! Ornstein–Zernike equations",  C. T. Kelley and B. Montgomery Pettitt,
! Journal of Computational Physics, Volume 197, Issue 2, 1 July 2004,
! Pages 491-501 [K&M-P] (check for typos in the paper!).

! Let h(k) be the 3d Fourier transform of h(r), with
!   h(k) = int d^3r exp(-ik.r) h(r)
! then one can show
!   h(k) = (4 pi / k) int_0^infty dr sin(kr) r h(r)
! which defines the forward Fourier-Bessel transform.  Likewise, if
!   h(r) = int d^3k/(2pi)^3 exp(ik.r) h(k)
! then one can show
!   h(r) = (1 / 2 pi^2 r) int_0^infty dk sin(kr) k h(k)
! which defines the backward Fourier-Bessel transform.
! The Ornstein-Zernike eqn in these terms is h(k) = c(k) + rho h(k) c(k)
! where rho is the density (with suitable generalisation to multicomponent
! systems).  This confirms the relevance of the standard choice of
! normalisation of the 3d Fourier transform pair which puts 1/(2 pi)^3 into
! the back transform.

! The discrete version of the forward Fourier-Bessel transform is
!   h_j = (2 pi delta / k_j)  *  2 sum_i=1^(n-1) r_i h_i sin(pi i j / n)
! where j = 1..n-1.  In this delta = L / n is the spacing in r, and
! r_i = i*delta.  The spacing in k is delta_k = pi / L = pi / (n delta),
! and k_j = j*delta_k.  The quantity 2 sum_i=1^(n-1) r_i h_i sin(pi i j / n)
! is computed by calling the FFTW routine RODFT00 on r_i h_i, with
! length n-1 -- see fftw_test.f90.
! Note that both the arrays in real and reciprocal space are of length n-1.

! The discrete version of the backward Fourier-Bessel transform
!   h(r) = (1 / 2 pi^2 r) int_0^infty dk sin(kr) k h(k)
! is
!   h_i = (delta_k / (2 pi)^2 r_i)  *  2 sum_j=1^(n-1) k_j h_j sin(pi i j / n)
! The second factor is computed by calling RODFT00 on k_j h_j, with length n-1.
! With reference to K&M-P, i and j here are i-1 and j-1 in this paper,
! and n = N-1.  K&M-P also suggest that c and e at r = 0 are evaluated by
! simple linear extrapolation,
!   c_0 = 2 c_1 - c_2,  e_0 = 2 e_1 - e_2
! and c and e at r = L corresponding to i = n are zero,
!   c_n = e_n = 0.
! Note that the indirect correlation function 'e' is called gamma by
! Vrbka et al [JCP v131, 154109 (2009)], and 'b' in Hansen and McDonald.

! Species pair index (i, j) is mapped to function index ij using the
! following rule:
!
!  if (i <= j) then ij = i + j*(j-1)/2
!  else ij = j + i*(i-1)/2
!
! This effectively labels the upper triangular entries in a symmetric
! matrix as (i labels the row, j labels the columm)
!
!     | 1  2  3  4 ...
!  ---------------------
!   1 | 1  2  4  7 ...
!   2 |    3  5  8 ...
!   3 |       6  9 ...
!   4 |         10 ...
!  ...|            ...
!
! A common idiom used below to walk through the function index is
!
!  do j = 1, ncomp
!    do i = 1, j
!      ij = i + j*(j-1)/2
!      .....
!

module wizard

  implicit none

  include "fftw3.f"

  double precision, parameter :: &
       & pi = 3.141592653589793d0, &
       & twopi = 2.0d0 * pi, &
       & fourpi = 4.0d0 * pi

  integer :: &
       & verbose = 0,     & ! how much info to generate
       & cold_start = 1,  & ! if the solver needs initialising
       & start_type = 3,  & ! how to initialise in a cold start
       & auto_fns = 1,    & ! whether to calculate things at end
       & model_type = 0,  & ! which potential was (last) chosen
       & istep,           & ! current position in iterative solver
       & ng = 4096,       & ! grid size
       & ncomp = 1,       & ! number of chemical components
       & nfnc = 0,        & ! number of functions, = ncomp (ncomp + 1) / 2
       & nps = 6,         & ! number of previous states used in Ng method
       & npic = 6,        & ! number of Picard steps
       & maxsteps = 100     ! max number of steps to take for convergence
  
  integer*8 :: plan  ! FFTW plan for fast discrete sine transforms

  double precision :: &
       & deltar = 0.01,   & ! real space grid spacing
       & deltak,          & ! reciprocal space grid spacing (computed)
       & error,           & ! difference between current and previous solns
       & alpha = 0.2,     & ! Picard method, fraction of new solution
       & tol = 1.0d-12,   & ! Error tolerance for claiming convergence
       & rc = 1.0,        & ! short-range DPD repulsion range
       & lb = 0.0,        & ! long-range Coulomb coupling length
       & sigma = 1.0,     & ! see below (*)
       & sigmap = 1.0,    & ! +- long-range Coulomb smearing length (URPM)
       & kappa = -1.0,    & ! +- long-range Coulomb smoothing parameter (RPM)
       & rgroot = 1.0,    & ! linear charge smearing range (Groot)
       & lambda = 1.0,    & ! Slater charge smearing range (exact)
       & beta = 1.0,      & ! Slater charge smearing range (approx)
       & cf_mf, cf_xc,    & ! the virial route pressure contributions ..
       & cf_gc, press,    & ! .. and the virial route pressure
       & comp, comp_xc,   & ! compressibility, and excess
       & fvex, fnex,      & ! excess free energy, density and per particle
       & un_mf, un_xc,    & ! energy per particle contributions
       & un, uv             ! energy per particle and density

  ! (*) sigma is used both for the long-range Coulomb smearing length
  ! for the soft potentials, and for hard core diameter for RPM models
  
  double precision, allocatable :: &
       & rho(:),            & ! density array
       & z(:),              & ! valence array
       & arep(:, :),        & ! repulsion amplitude array
       & diam(:),           & ! hard core diameter array
       & tp(:),             & ! mean field pressure contribution
       & tu(:),             & ! mean field energy contribution
       & tl(:),             & ! mean field long range potential
       & muex(:),           & ! chemical potential array
       & c(:, :, :),        & ! direct correlation functions (dcfs)
       & e(:, :, :),        & ! indirect correlation functions (icfs)
       & hr(:, :, :),       & ! total correlation functions (tcfs)
       & ck(:, :),          & ! transform of dcfs
       & ek(:, :),          & ! transform of icfs
       & hk(:, :),          & ! transform of total correlation functions
       & sk(:, :, :),       & ! partial structure factors
       & ushort(:, :),      & ! short range potential in real space
       & expnegus(:, :),    & ! exp(-ushort) (includes hard cores)
       & dushort(:, :),     & ! derivative of the same
       & ulong(:, :),       & ! long range potential in real space
       & dulong(:, :),      & ! derivative of the same
       & ulongk(:, :),      & ! long range potential in reciprocal space
       & r(:), k(:),        & ! r and k grids
       & fftwx(:), fftwy(:)   ! arrays for fast discrete sine transform

contains

  subroutine initialise
    implicit none
    integer i

    nfnc = ncomp * (ncomp + 1) / 2

    allocate(rho(ncomp))
    allocate(z(ncomp))
    allocate(muex(ncomp))
    allocate(arep(ncomp, ncomp))
    allocate(diam(nfnc))
    allocate(tp(nfnc))
    allocate(tu(nfnc))
    allocate(tl(nfnc))
    allocate(c(ng-1, nfnc, nps))
    allocate(e(ng-1, nfnc, nps))
    allocate(hr(ng-1, ncomp, ncomp))
    allocate(ck(ng-1, nfnc))
    allocate(ek(ng-1, nfnc))
    allocate(hk(ng-1, nfnc))
    allocate(sk(ng-1, ncomp, ncomp))
    allocate(ushort(ng-1, nfnc))
    allocate(expnegus(ng-1, nfnc))
    allocate(dushort(ng-1, nfnc))
    allocate(ulong(ng-1, nfnc))
    allocate(dulong(ng-1, nfnc))
    allocate(ulongk(ng-1, nfnc))
    allocate(r(ng-1))
    allocate(k(ng-1))
    allocate(fftwx(ng-1))
    allocate(fftwy(ng-1))

    ! Default values

    rho = 0.0
    arep = 0.0
    z = 0.0
    diam = 0.0

    ! Make grids

    deltak = pi / (dble(ng) * deltar)

    r = (/ (i*deltar, i=1, ng-1) /)
    k = (/ (i*deltak, i=1, ng-1) /)

    ! Make the FFTW plan

    call dfftw_plan_r2r_1d(plan, ng-1, fftwx, fftwy, &
         & FFTW_RODFT00, FFTW_ESTIMATE)

  end subroutine initialise

! Write out parameters for system and potentials.

  subroutine write_params
    implicit none
    integer :: i
    print *, '====================================================='
    print *, 'GRID DETAILS'
    print *, ' ng = ', ng, ' ncomp = ', ncomp, ' nfnc = ', nfnc, ' nps = ', nps
    print *, ' deltar = ', deltar, ' deltak = ', deltak
    print *, ' deltar*deltak*ng/pi = ', deltar*deltak/pi*dble(ng)
    print *, ' r(ng-1) = ', r(ng-1), ' k(ng-1) = ', k(ng-1)
    print *, 'POTENTIAL DETAILS (model type', model_type, ')'
    if (model_type.eq.0) then
       print *, 'No potential has been selected'
    else if (model_type.lt.10) then
       print *, 'DPD potential was selected, matrix A = '
       do i = 1, ncomp
          print *, ' ', arep(i, :)
       end do
       print *, ' valencies, z = ', z
       print *, ' rc = ', rc, ' lb = ', lb
       if (model_type.eq.1) then
          print *, ' Gaussian smearing, sigma = ', sigma
       end if
       if (model_type.eq.2) then
          print *, ' Bessel smearing, sigma = ', sigma
       end if
       if (model_type.eq.3) then
          print *, ' linear (Groot) smearing, R = ', rgroot
          print *, ' equivalent Gaussian sigma = ', &
               & sqrt(2.0d0/15.0d0) * rgroot
       end if
       if (model_type.eq.4) then
          print *, ' Slater smearing (exact), lambda = ', lambda
          print *, ' equivalent Gaussian sigma = ', lambda
       end if
       if (model_type.eq.5) then
          print *, ' Slater smearing (approx), beta = ', beta
          print *, ' equivalent lambda = ', 0.625d0 / beta, &
               & ' ( 1 / beta = ', 1.0d0 / beta, ' )'
          print *, ' equivalent Gaussian sigma = ', sqrt(0.5d0) / beta
       end if
    else if (model_type.lt.20) then
       if (model_type.eq.10) then
          print *, 'softened URPM potential was selected with ushort unused'
       else
          print *, 'softened URPM potential was selected with ushort used'
       end if
       print *, ' lb = ', lb, ' sigma = ', sigma, ' sigmap = ', sigmap
    else if (model_type.lt.30) then
       if (model_type.eq.20) then
          print *, 'softened RPM potential was selected with ushort unused'
       else
          print *, 'softened RPM potential was selected with ushort used'
       end if
       if (kappa.lt.0) then
          print *, ' lb = ', lb, ' sigma = ', sigma, ' kappa -> infinity'
       else
          print *, ' lb = ', lb, ' sigma = ', sigma, ' kappa = ', kappa
       end if
    else
       print *, 'Undefined potential'
    end if
    print *, 'SYSTEM DETAILS'
    print *, ' rho = ', rho
    if (sum(rho).gt.0.0) then
       print *, ' x = ', rho(:) / sum(rho)
    end if
    print *, ' sum(rho) = ', sum(rho)
    print *, ' sum(rho*z) = ', sum(rho(:)*z(:))
    print *, '====================================================='
  end subroutine write_params

! Build the potential arrays, with parameters rc and arep(:,:) for the
! short-range DPD repulsion, and lb, sigma and z(:) for the long-range
! Gaussian-smeared Coulomb part.  A factor beta = 1/kT is implicit in
! these definitions.  The parameter charge_type is Gaussian (1),
! Bessel (2), Groot (3), Slater (exact) (4), Slater (approximate) (5).

  subroutine dpd_potential(charge_type)
    implicit none
    integer, intent(in) :: charge_type
    integer :: i, j, ij, irc
    double precision :: aa(nfnc), zz(nfnc)
    double precision :: rootpi

    rootpi = sqrt(pi)

    ! Sort out some recoded potential parameters.  Also force the
    ! amplitude matrix to be symmetric, set by upper triangular
    ! entries.  This enforces the rule in the documentation and also
    ! simplifies the printing.

    do j = 1, ncomp
       do i = 1, j
          ij = i + j*(j-1)/2
          aa(ij) = arep(i, j)
          zz(ij) = z(i) * z(j)
          if (i.lt.j) then
             arep(j, i) = arep(i, j)
          end if
       end do
    end do

    irc = nint(rc/deltar)

    ! Leave out the amplitude, then the function can be re-used
    ! (see below)

    ushort(:,nfnc) = 0.0d0
    ushort(1:irc,nfnc) = 0.5d0 * (1.0d0 - r(1:irc)/rc)**2

    dushort(:,nfnc) = 0.0d0
    dushort(1:irc,nfnc) = - (1.0d0 - r(1:irc)/rc) / rc

    ! Gaussian charges

    if (charge_type .eq. 1) then

       ulong(:,nfnc) = (lb / r) * erf(0.5d0*r/sigma)

       ulongk(:,nfnc) = (fourpi * lb / k**2) * exp(-k**2*sigma**2)

       dulong(:,nfnc) = lb * exp(-0.25d0*r**2/sigma**2) &
            & / (rootpi * r * sigma) &
            & - lb * erf(0.5d0*r/sigma) / r**2

    end if

    ! Bessel charges

    if (charge_type .eq. 2) then

       ulong(:,nfnc) = (lb / r) * (1.0d0 - exp(-r/sigma))

       ulongk(:,nfnc) = (fourpi * lb / k**2) &
            & * 1.0d0 / (1.0d0 + k**2*sigma**2)

       dulong(:,nfnc) = - (lb / r**2) * (1.0d0 - exp(-r/sigma) &
            & * (1.0d0 + r / sigma))

    end if

    ! Linear charge smearing as in Groot [JCP v118, 11265 (2003)].
    ! Note we do not give the real space part here hence the
    ! thermodynamic calculations will be wrong.

    if (charge_type .eq. 3) then

       ulong(:,nfnc) = 0.0d0

       ulongk(:,nfnc) = (fourpi * lb / k**2) * 144.0d0 &
            & * (2.0d0 - 2.0d0*cos(k*rgroot) &
            &      - k*rgroot*sin(k*rgroot))**2 &
            &                  / (k**8 * rgroot**8)

       dulong(:,nfnc) = 0.0d0

    end if

    ! Slater charge smearing as in Gonzales-Melchor et al, [JCP v125,
    ! 224107 (2006)] with exact expression for interaction (here
    ! translated into reciprocal space).

    if (charge_type .eq. 4) then

       ulong(:,nfnc) = (lb / r) * (1.0d0 - exp(-2.0d0*r/lambda) &
            & * (1.0d0 + 1.375d0*r/lambda + 0.75d0*r**2/lambda**2 &
            &   + r**3/(6.0d0*lambda**3)) )

       ulongk(:,nfnc) = (fourpi * lb / k**2) * &
            & 1.0d0 / (1.0d0 + k**2*lambda**2/4.0d0)**4

       dulong(:,nfnc) = - (lb / r**2) * (1.0d0 - exp(-2.0d0*r/sigma) &
            & * (1.0d0 + 2.0d0*r/lambda + 2.0d0*r**2/lambda**2 &
            &     + 7.0d0*r**3/(6.0d0*lambda**3) + r**4/(3.0d0*lambda**4)) )

    end if

    ! Slater charge smearing as in Gonzales-Melchor et al, [JCP v125,
    ! 224107 (2006)] with original approximate expression (here
    ! translated into reciprocal space).

    if (charge_type .eq. 5) then

       ulong(:,nfnc) = (lb / r) * (1.0d0 - exp(-2*beta*r) * &
            & (1.0d0 + beta*r) )

       ulongk(:,nfnc) = (fourpi * lb / k**2) * &
            & 1.0d0 / (1.0d0 + k**2/(4.0d0*beta**2))**2

       dulong(:,nfnc) = - (lb / r**2) * (1.0d0 - exp(-2.0d0*beta*r) &
            & * (1.0d0 + 2.0d0*beta*r*(1.0d0 + beta*r)) )

    end if

    ! Generate the pair potentials by walking through the functions.
    ! In the final step we correctly normalise the final function.

    do i = 1, nfnc
       ushort(:,i)  = aa(i) * ushort(:,nfnc)
       dushort(:,i) = aa(i) * dushort(:,nfnc)
       ulong(:,i)   = zz(i) * ulong(:,nfnc)
       ulongk(:,i)  = zz(i) * ulongk(:,nfnc)
       dulong(:,i)  = zz(i) * dulong(:,nfnc)
    end do

    ! These individual species-pair contributions to the mean field
    ! compressibility factor and the mean-field internal energy per
    ! particle can be calculated analytically for the DPD potential.

    tp = pi * rc**3 * aa / 30.0
    tu = tp
    tl = 0.0d0

    ! Generate auxiliary function

    expnegus = exp(-ushort)

    ! Record the model type

    model_type = charge_type

  end subroutine dpd_potential

! Build the potential arrays for the softened URPM (Gaussian charges),
! with parameters lb, sigma and sigmap.  This expects ncomp = 2, and
! will set z(1) = 1, z(2) = -1.  The parameter (0 or 1) controls
! whether ushort is used or not.

  subroutine soft_urpm_potential(use_ushort)
    implicit none
    integer, intent(in) :: use_ushort
    double precision :: rootpi

    rootpi = sqrt(pi)

    if (ncomp.ne.2) then
       print *, 'oz_mod.f90: soft_urpm_potential: ncomp = ', ncomp, &
            & '(should be ncomp = 2)'
       stop
    end if

    z(1) = 1; z(2) = -1;

    ulong(:,1) = lb * erf(0.5d0*r/sigma) / r

    ulongk(:,1) = fourpi * lb * exp(-k**2*sigma**2) / k**2

    dulong(:,1) = lb * exp(-0.25d0*r**2/sigma**2) / (rootpi * r * sigma) &
         & - lb * erf(0.5d0*r/sigma) / r**2

    ulong(:,2) = - lb * erf(0.5d0*r/sigmap) / r

    ulongk(:,2) = - fourpi * lb * exp(-k**2*sigmap**2) / k**2

    dulong(:,2) = - lb * exp(-0.25d0*r**2/sigmap**2) / (rootpi * r * sigmap) &
         & + lb * erf(0.5d0*r/sigmap) / r**2

    ulong(:,3) = ulong(:,1)
    ulongk(:,3) = ulongk(:,1)
    dulong(:,3) = dulong(:,1)

    ushort(:,:) = 0.0d0
    dushort(:,:) = 0.0d0

    if (use_ushort.ne.0) then
       ushort(:,2) = ulong(:,2) + ulong(:,1)
       dushort(:,2) = dulong(:,2) + dulong(:,1)
       ulong(:,2) = - ulong(:,1)
       ulongk(:,2) = - ulongk(:,1)
       dulong(:,2) = - dulong(:,1)
    end if

    ! These individual species-pair contributions to the mean field
    ! compressibility factor and the mean-field internal energy per
    ! particle can be calculated analytically for the URPM potential.
    ! These are the same whether using ushort or not, as they are
    ! defined in terms of the total potential.

    tp(1) = 0.0d0
    tp(2) = 2*pi*lb*(sigmap**2 - sigma**2)
    tp(3) = 0.0d0

    tu = tp

    ! If not using ushort, we are off the symmetry point condition and
    ! the contribution of the long range part should be incorporated
    ! into the compressibility and chemical potential expressions.

    if (use_ushort.eq.0) then
       tl = 2.0*tp
    else
       tl = 0.0d0
    end if

    ! Generate auxiliary function

    expnegus = exp(-ushort)

    ! Record the model type

    model_type = 10 + use_ushort

  end subroutine soft_urpm_potential

! Build the potential arrays for the softened RPM (charged hard
! spheres) with parameters lb, sigma and kappa.  This expects ncomp =
! 2, and will set z(1) = 1, z(2) = -1, and hard core diameters to
! sigma.  The parameter (0 or 1) controls whether ushort is used or
! not.  Using kappa < 0 implies the pure RPM case (kappa -> infinity).

  subroutine soft_rpm_potential(use_ushort)
    implicit none
    integer, intent(in) :: use_ushort
    integer :: i, irc
    double precision :: rootpi

    rootpi = sqrt(pi)

    if (ncomp.ne.2) then
       print *, 'oz_mod.f90: soft_urpm_potential: ncomp = ', ncomp, &
            & '(should be ncomp = 2)'
       stop
    end if

    z(1) = 1; z(2) = -1;
    diam(1) = sigma
    diam(2) = sigma
    diam(3) = sigma

    ulong(:,1) = lb / r
    ulongk(:,1) = fourpi * lb / k**2
    dulong(:,1) = - lb / r**2

    if (kappa.gt.0.0) then
       ulong(:,2) = - lb * erf(kappa*r) / r
       ulongk(:,2) = - fourpi * lb * exp(-k**2/(4.0d0*kappa**2)) / k**2
       dulong(:,2) = - 2.0d0*kappa*lb * exp(-kappa**2*r**2) / (rootpi * r) &
            & + lb * erf(kappa*r) / r**2
    else
       ulong(:,2) = - lb / r
       ulongk(:,2) = - fourpi * lb / k**2
       dulong(:,2) = lb / r**2
    end if

    ulong(:,3) = ulong(:,1)
    ulongk(:,3) = ulongk(:,1)
    dulong(:,3) = dulong(:,1)

    ushort(:,:) = 0.0d0
    dushort(:,:) = 0.0d0

    if (use_ushort.ne.0) then
       ushort(:,2) = ulong(:,2) + ulong(:,1)
       dushort(:,2) = dulong(:,2) + dulong(:,1)
       ulong(:,2) = - ulong(:,1)
       ulongk(:,2) = - ulongk(:,1)
       dulong(:,2) = - dulong(:,1)
    end if

    ! These are the analytic contributions to the thermodynamics.

    tp = 0.0d0
    tu = 0.0d0
    tl = 0.0d0

    if (kappa.gt.0) then
       tp(2) = pi*lb * ( sigma * exp(-kappa**2*sigma**2) / (kappa*rootpi) &
            & + (1/(2.0d0*kappa**2) - sigma**2/3.0d0) * erfc(kappa*sigma) )
       tu(2) = pi*lb * ( sigma * exp(-kappa**2*sigma**2) / (kappa*rootpi) &
            & + (1/(2.0d0*kappa**2) - sigma**2) * erfc(kappa*sigma) )
       if (use_ushort.eq.0) then
          tl(2) = pi*lb / kappa**2
       end if
    end if

    ! Generate auxiliary function

    expnegus = exp(-ushort)

    ! Impose the hard core condition

    do i = 1, nfnc
       irc = nint(diam(i) / deltar)
       ushort(1:irc, i) = 0.0d0
       ulong(1:irc, i) = 0.0d0
       dushort(1:irc, i) = 0.0d0
       expnegus(1:irc, i) = 0.0d0
    end do

    ! Record the model type

    model_type = 20 + use_ushort

  end subroutine soft_rpm_potential

! The next routine solves the Ornstein-Zernicke equation to determine
! e = h - c, given c.  We re-partition the long range part of the
! potential so that the routine actually calculates c' = c + Ulong
! and e' = e - Ulong.  This is because the Fourier transform of
! Ulong can be computed in closed form.  Note h = e + c = e' + c'.

  subroutine oz_solve
    implicit none
    integer :: i1, i, j, ij, ik, irc
    integer :: perm(ncomp)
    double precision :: &
         & a(ncomp, ncomp), b(ncomp, ncomp), &
         & cmat(ncomp, ncomp), umat(ncomp, ncomp), rhomat(ncomp, ncomp), &
         & m0(ncomp, ncomp), unita(ncomp, ncomp)

    i1 = mod(istep-1, nps) + 1

    ! Forward transform the real space functions c, to the reciprocal
    ! space functions ck.

    do i=1, nfnc
       fftwx(1:ng-1) = r(1:ng-1) * c(1:ng-1, i, i1)
       call dfftw_execute(plan)
       ck(1:ng-1, i) =  (twopi * deltar) * fftwy(1:ng-1) / k(1:ng-1)
    end do

    if (ncomp .eq. 1) then

       ! In the one component case OZ inversion is straightforward.
       ! Note the implicit indexing on wavevector k.

       ek(:, 1) = ( ck(:, 1) - ulongk(:, 1) ) &
            & / ( 1.0d0 - rho(1) * (ck(:, 1) - ulongk(:, 1)) ) &
            & - ck(:, 1)

    else ! Multicomponent OZ inversion

       ! First set up a unit matrix, and the diagonal R matrix -- see
       ! the documentation for the math here.

       rhomat = 0.0d0
       unita = 0.0d0

       do i = 1, ncomp
          rhomat(i,i) = rho(i)
          unita(i,i) = 1.0d0
       end do

       ! Do the matrix calculations for each wavevector k.

       do ik = 1, ng-1

          ! Unpack the reciprocal space functions into matrices.

          do j = 1, ncomp
             do i = 1, j
                ij = i + j*(j-1)/2
                cmat(i, j) = ck(ik, ij)
                umat(i, j) = ulongk(ik, ij)
                if (i.lt.j) then
                   cmat(j, i) = cmat(i, j)
                   umat(j, i) = umat(i, j)
                end if
             end do
          end do

          ! The following matrices are constructed:
          !   M0 = (C - beta UL) . R
          !   A = I - (C - beta UL) . R
          !   B = (C - beta UL) . R . C - beta UL
          ! (note that beta = 1/kT = 1 is not written explicitly)

          m0 = matmul(cmat - umat, rhomat)
          a = unita - m0
          b = matmul(m0, cmat) - umat

          ! Solve the equation A.X = B so that
          ! X = [I - (C - beta U) . R]^(-1) . [(C - beta U) . R . C - beta U]

          call axeqb_reduce(a, ncomp, b, ncomp, perm, irc)

          if (irc.gt.0) then
             print *, 'oz_solve(oz_mod): axeqb_reduce returned irc = ', irc
             stop
          end if

          ! Now X(I, :) = B(PERM(I), :) is the new estimate for the
          ! reciprocal space functions ek.  They are built from the
          ! upper triangle of the matrix.

          do j = 1, ncomp
             do i = 1, j
                ij = i + j*(j-1)/2
                ek(ik, ij) = b(perm(i), j)
             end do
          end do

       end do ! loop over k vectors

    end if ! select single component or multicomponent case

    ! Do the Fourier back transforms

    do i = 1, nfnc
       fftwx(1:ng-1) = k(1:ng-1) * ek(1:ng-1, i)
       call dfftw_execute(plan)
       e(1:ng-1, i, i1) =  (deltak / twopi**2) * fftwy(1:ng-1) / r(1:ng-1)
    end do

  end subroutine oz_solve

! This routine solves an alternate version of the Ornstein-Zernicke
! equation to determine c and e from h.  In practice as always we
! actually calculate c' = c + Ulong and e' = e - Ulong.  Note h = e +
! c = e' + c'.

  subroutine oz_solve2
    implicit none
    integer :: i1, i, j, ij, ik, irc
    integer :: perm(ncomp)
    double precision :: &
         & a(ncomp, ncomp), b(ncomp, ncomp), &
         & h(ng-1, nfnc), hmat(ncomp, ncomp), &
         & rhomat(ncomp, ncomp), unita(ncomp, ncomp)

    i1 = mod(istep-1, nps) + 1

    ! This is called after the correlation function h_ij have been
    ! calculated (using RPA).  We therefore unpack the matrix
    ! functions and Fourier transform to reciprocal space

    do j = 1, ncomp
       do i = 1, j
          ij = i + j*(j-1)/2
          h(:, ij) = hr(:, i, j)
       end do
    end do

    do i=1, nfnc
       fftwx(1:ng-1) = r(1:ng-1) * h(1:ng-1, i)
       call dfftw_execute(plan)
       hk(1:ng-1, i) =  (twopi * deltar) * fftwy(1:ng-1) / k(1:ng-1)
    end do

    if (ncomp .eq. 1) then

       ! In the one component case the OZ solution is trivial.

       ck(:, 1) = hk(:, 1) / (1.0d0 + rho(1) * hk(:, 1)) &
            & + ulongk(:, 1)

    else ! Multicomponent OZ solution

       ! As above set up a unit matrix, and the diagonal R matrix

       rhomat = 0.0d0
       unita = 0.0d0

       do i = 1, ncomp
          rhomat(i,i) = rho(i)
          unita(i,i) = 1.0d0
       end do

       ! Do the matrix calculations for each wavevector k.

       do ik = 1, ng-1

          ! Convert the reciprocal space functions into matrices.

          do j = 1, ncomp
             do i = 1, j
                ij = i + j*(j-1)/2
                hmat(i, j) = hk(ik, ij)
                if (i.lt.j) then
                   hmat(j, i) = hmat(i, j)
                end if
             end do
          end do

          ! Construct A = I + H . R, and B = H

          a = unita + matmul(hmat, rhomat)
          b = hmat

          ! Solve A.X = B so that X = (I + H.R)^(-1) . H.

          call axeqb_reduce(a, ncomp, b, ncomp, perm, irc)

          ! Now compute C = (I + H.R)^(-1) . H + beta UL
          ! (map back to functions, and unravel the pivoting)

          do j = 1, ncomp
             do i = 1, j
                ij = i + j*(j-1)/2
                ck(ik, ij) = b(perm(i), j) + ulong(ik, ij)
             end do
          end do

       end do

    end if

    ! Do the Fourier back transforms and calculate e = h - c.  This
    ! means the results can be used in the structure and
    ! thermodynamics routines as though they had come from the HNC
    ! solution.

    do i = 1, nfnc

       fftwx(1:ng-1) = k(1:ng-1) * ck(1:ng-1, i)
       call dfftw_execute(plan)
       c(1:ng-1, i, i1) =  (deltak / twopi**2) * fftwy(1:ng-1) / r(1:ng-1)

       ek(:, i) = hk(:, i) - ck(:, i)
       e(:, i, i1) = h(:, i) - c(:, i, i1)

    end do

  end subroutine oz_solve2

! This routine implements the HNC closure expressed as c = exp(-beta v
! + e) - e - 1 where e = h - c is the indirect correlation function, c
! is the direct correlation function from the Ornstein-Zernicke
! relation, h = g - 1, and g is the pair distribution function.  One
! can show this is equivalent to g = exp(-v + h - c) in Hansen +
! McDonald.  As above, the routine actually works with c' = c + Ulong
! and e' = e - Ulong where Ulong is the long-range part of the
! potential for which the Fourier transform is simple.  This means
! that 'v' in the above expression is the short-range part of the
! potential only.

  subroutine picard_method
    implicit none
    integer :: i1, i0, i

    istep = istep + 1
    i1 = mod(istep-1, nps) + 1
    i0 = i1 - 1; if (i0.eq.0) i0 = nps

    do i = 1, nfnc
       !!  c(:,i,i1) = alpha * ( exp(- ushort(:,i) + e(:,i,i0)) &
       c(:,i,i1) = alpha * ( expnegus(:,i) * exp(e(:,i,i0)) &
            & - e(:,i,i0) - 1.0d0 ) &
            & + (1.0d0 - alpha) * c(:,i,i0)
    end do

  end subroutine picard_method

! The next routine implements the Ng method [K-C Ng, J. Chem. Phys.
! v61, 2680 (1974)] as an accelerated solver for the HNC closure.

  subroutine ng_method
    implicit none
    integer :: i, i1, i0, j, j1, j2, p, nd, icp
    double precision :: dc(ng-1,nfnc,nps-1), de(ng-1,nfnc,nps-1), &
         & a(nps-1,nps-1), x(nps-1), y(nps-1), yy, aux

    integer :: ipiv(nps-1), info  ! DSYSV stuff
    double precision :: work(100) ! DSYSV stuff

    istep = istep + 1
    i1 = mod(istep-1, nps) + 1
    i0 = i1 - 1; if (i0 .eq. 0) i0 = nps

    if (istep .le. nps) then
       nd = istep - 2
    else
       nd = nps - 1
    end if

    do p = 1, nd
       j1 = i0 - p
       if( j1 .le. 0) j1 = nps + j1
       dc(:,:,p) = c(:,:,i0) - c(:,:,j1)
       de(:,:,p) = e(:,:,i0) - e(:,:,j1)
    end do

    a(:,:) = 0.0d0

    x(:) = 0.0d0

    do icp = 1, nfnc
       do j = 1, ng-1
          !!  aux = exp( - ushort(j,icp) + e(j,icp,i0)) - 1.0d0
          aux = expnegus(j,icp) * exp(e(j,icp,i0)) - 1.0d0

          do j1 = 1, nd
             y(j1) = aux * de(j,icp,j1) - dc(j,icp,j1)
          end do

          yy = aux - e(j,icp,i0) - c(j,icp,i0)

          do j1 = 1, nd
             do j2 = j1, nd
                a(j1,j2) = a(j1,j2) + y(j1) * y(j2)
             end do
             x(j1) = x(j1) + y(j1) * yy
          end do

       end do
    end do

    call DSYSV( 'U', nd, 1, a, nps-1, ipiv, x, nps-1, work, &
         & 100, info)

    if (info .gt. 0) then
       print *, 'det=0', (x(i),i=1,nd)
    endif

    do icp = 1, nfnc
       do j = 1, ng-1
          aux = e(j,icp,i0)
          do j1 = 1, nd
             aux = aux - de(j,icp,j1) * x(j1)
          end do
          !!  c(j,icp,i1) = exp( - ushort(j,icp) + aux) - aux - 1.0d0
          c(j,icp,i1) = expnegus(j,icp) * exp(aux) - aux - 1.0d0
       end do
    end do

  end subroutine ng_method

! Calculate the difference between the direct correlation functions
! for the current and previous iteration, used as a convergence test;
! return answer in variable 'error'.

  subroutine conv_test
    implicit none
    integer i1, i0
    i1 = mod(istep - 1, nps) + 1
    i0 = i1 - 1; if (i0 .eq. 0) i0 = nps
    error = sqrt(deltar * sum( (c(:, :, i1) - c(:, :, i0))**2 ))
  end subroutine conv_test

! Basic driver routine for solving HNC: take a number of Picard
! iterations to pump-prime the Ng method.  Stop when error is less
! than tolerance, or when exceed maximum number of iterations.  The
! flag cold_start indicates whether the direct correlation function
! should be re-initialised.  The initial guess to the direct
! correlation function is either zero (start_type = 1), or c = -
! Ushort (start_type = 2), or c = e^(-Ushort)-1 (start_type = 3).  Any
! of these should do in principle, but the initial convergence may be
! different.  Note from above that c is actually defined c' = c +
! Ulong, ie with the long-range part of the potential added.

  subroutine hnc_solve
    implicit none
    integer :: i
    if (cold_start.eq.1) then
       istep = 1
       if (start_type.eq.1) c(:,:,1) = 0.0
       if (start_type.eq.2) c(:,:,1) = - ushort(:,:)
       if (start_type.eq.3) c(:,:,1) = expnegus(:,:) - 1.0
       cold_start = 0
       if (verbose.eq.1) then
          if (start_type.eq.1) print *, "cold start c' = 0"
          if (start_type.eq.2) print *, "cold start c' = -v'"
          if (start_type.eq.3) print *, "cold start c' = e^(-v')-1"
       end if
    else
       if (verbose.eq.1) then
          print *, "warm start c' = previous c'"
       end if
    end if
    call oz_solve
    do i = 1, maxsteps
       if (i .le. npic) then
          call picard_method
       else
          call ng_method
       end if
       call oz_solve
       call conv_test
       if (verbose.eq.1) then
          if (i .le. npic) then
             print *, i, "Picard, error = ", error
          else
             print *, i, "    Ng, error = ", error
          end if
       end if
       if (error .lt. tol) exit
    end do
    if (error .gt. tol) then
       print *, "oz_mod.f90: solve_problem: error > tol"
    else
       if (auto_fns.eq.1) then
          call make_pair_functions
          call make_structure_factors
          call make_thermodynamics
       end if
    end if
  end subroutine hnc_solve

! Given the HNC machinery, the implementation of the RPA is almost
! completely trivial and corresponds to one iteration through the
! Ornstein-Zernike solver given the choice c = - Ushort (HNC
! start_type = 2).

  subroutine rpa_solve
    implicit none
    istep = 1
    c(:,:,1) = - ushort(:,:)
    call oz_solve
    if (auto_fns.eq.1) then
       call make_pair_functions
       call make_structure_factors
       call make_thermodynamics
    end if
    if (verbose.eq.1) then
       print *, "RPA solution, c' = - v'"
    end if
  end subroutine rpa_solve

! The EXP approximation is a development of the RPA approximation, in
! which h --> exp(h)-1.  A full solution requires a follow-up round
! trip through another version of the Ornstein-Zernike relation, to
! obtain the direct and indirect correlation functions.

  subroutine exp_solve
    implicit none
    istep = 1
    c(:,:,1) = - ushort(:,:)
    call oz_solve
    call make_pair_functions
    hr(:,:,:) = exp(hr(:,:,:)) - 1.0d0
    call oz_solve2
    if (auto_fns.eq.1) then
       call make_structure_factors
       call make_thermodynamics
    end if
    if (verbose.eq.1) then
       print *, "EXP solution"
    end if
  end subroutine exp_solve

! Construct the structure factors out of the transform of the total
! correlation function.  Note that ck and ek are available after a
! call to the OZ solver.

  subroutine make_structure_factors
    implicit none
    integer :: i, j, ij
    hk = ck + ek
    do j = 1, ncomp
       do i = 1, j
          ij = i + j*(j-1)/2
          if (i.eq.j) then
             sk(:, i, i) = rho(i) * (1.0 + rho(i) * hk(:, ij))
          else
             sk(:, i, j) = rho(i) * rho(j) * hk(:, ij)
             sk(:, j, i) = sk(:, i, j)
          end if
       end do
    end do
  end subroutine make_structure_factors

! Construct the total correlation functions out of the direct
! correlation functions.  Note that the above routines actually works
! with c' = c + Ulong and e' = e - Ulong where Ulong is the long-range
! part of the potential, but h = g - 1 = e + c = e' + c'.  The pair
! correlation functions are g = 1 + h - the addition of '1' is left
! for the user to implement.

  subroutine make_pair_functions
    implicit none
    integer :: i, j, ij, i1
    i1 = mod(istep-1, nps) + 1
    do j = 1, ncomp
       do i = 1, j
          ij = i + j*(j-1)/2
          hr(:, i, j) = c(:, ij, i1) + e(:, ij, i1)
          hr(:, j, i) = hr(:, i, j)
       end do
    end do
  end subroutine make_pair_functions

! Calculate various thermodynamics properties by spatial integration
! (as contrasted to thermodynamic integration).  We use the trapezium
! rule, taking account where necessary the end-point values (see
! intro), h_0 = 2 h_1 - h_2 at r = 0 (i = 0), and h_n = 0 at r = L (i
! = ng).  Also we have r = i*deltar for i = 1 to ng-1.  Note that the
! above routines actually work with c' = c + Ulong and e' = e - Ulong
! where Ulong is the long-range part of the potential, so we have h =
! g - 1 = e + c = e' + c'.  See also Vrbka et al, J. Chem. Phys. 131,
! 154109 (2009).
!
! The mean-field thermodynamic expressions can often be obtained
! analytically from the potential. In this routine they are calculated
! from species pair contributions, which are themselves calculated in
! the potential routines (which does not have access to the
! densities).

  subroutine make_thermodynamics
    implicit none
    integer :: i, j, ij, i1, irc
    double precision :: rhotot, r1, r2, g1, g2, gc
    double precision :: rhoxx(nfnc), t(nfnc)

    i1 = mod(istep-1, nps) + 1

    ! rhoxx is rho x_i x_j, doubled up for the off-diagonal components

    rhotot = sum(rho)

    do j = 1, ncomp
       do i = 1, j
          ij = i + j*(j-1)/2
          if (i.eq.j) then
             rhoxx(ij) = rho(i)**2 / rhotot
          else
             rhoxx(ij) = 2.0 * rho(i) * rho(j) / rhotot
          end if
       end do
    end do

    ! Calculate the various contributions to the virial-route
    ! pressure.  This is the mean field contribution.

    cf_mf = sum(rhoxx(:) * tp(:))

    ! Evaluate t_ij = - (2pi/3) int_d^inf d(U_ij)/dr h_ij r^3 dr. The
    ! contribution from both end-points r = 0 (i = 0) and r = L (i =
    ! ng) vanishes, hence the sum just consists of the middle part of
    ! the trapezium rule.  If we have set dushort + dulong = 0 within
    ! the hard core, the lower bound is taken care of automatically.

    do i = 1, nfnc
       t(i) =  - twopi * deltar * sum((dushort(:,i) + dulong(:,i)) &
            & * (c(:,i,i1) + e(:,i,i1)) * r(:)**3) / 3.0
    end do

    ! The correlation contribution is sum_ij rho x_i x_j t_ij.

    cf_xc = sum(rhoxx(:) * t(:))

    ! The contact contribution in the case of hard cores.  We
    ! extrapolate the contact value of the pair distribution function
    ! from the two nearest outside points.

    do i = 1, nfnc
       if (diam(i).gt.0.0) then
          irc = nint(diam(i) / deltar)
          r1 = r(irc+1)
          r2 = r(irc+2)
          g1 = 1.0 + c(irc+1,i,i1) + e(irc+1,i,i1)
          g2 = 1.0 + c(irc+2,i,i1) + e(irc+2,i,i1)
          gc = ( (g1 - g2) * diam(i) + g2 * r1 - g1 * r2 ) / (r1 - r2)
          t(i) = twopi * diam(i)**3 * gc / 3.0
       else
          t(i) = 0.0
       end if
    end do

    cf_gc = sum(rhoxx(:) * t(:))

    ! This is the final pressure.

    press = rhotot * (1.0 + cf_gc + cf_mf + cf_xc)

    ! Now we do the compressibility (not to be confused with the above
    ! compressibility factor).  The long range part is accounted for
    ! separately.

    ! Evaluate t_ij = 4 pi int_0^inf c_ij r^2 dr.  Again the
    ! contribution from both endpoints vanishes, hence the sum just
    ! consists of the middle part of the trapezium rule.

    do i = 1, nfnc
       t(i) = fourpi * deltar * sum(c(:,i,i1) * r(:)**2)
    end do

    ! The compressibility is 1 - sum_ij rho x_i x_j t_ij

    comp_xc = - sum(rhoxx(:) * (t(:) - tl(:)))
    comp = 1.0 + comp_xc

    ! Now we do the energy per particle.  First the mean-field
    ! contribution (per particle).

    un_mf = sum(rhoxx(:) * tu(:))

    ! Evaluate t_ij = 2 pi int_0^inf U_ij h_ij r^2 dr.  Note that the
    ! contribution from both end-points again vanishes.  If we have
    ! set ushort + ulong = 0 within the hard core, the lower bound is
    ! taken care of automatically.

    do i = 1, nfnc
       t(i) = twopi * deltar * sum((ushort(:,i) + ulong(:,i)) &
            & * (c(:,i,i1) + e(:,i,i1)) * r(:)**2)
    end do

    ! The extra contribution is sum_ij rho x_i x_j t_ij

    un_xc = sum(rhoxx(:) * t(:))

    ! This is the final energy per particle and density (per unit
    ! volume).

    un = un_mf + un_xc
    uv = rhotot * un

    ! Finally do the chemical potentials (this is valid ONLY for HNC).

    ! Evaluate t_ij = 4 pi int_0^inf (h_ij e_ij / 2 - c_ij) r^2 dr.
    ! Note that the contribution from both end-points again vanishes.
    ! Also we can use c' for c in the second term because charge
    ! neutrality causes the contribution from Ulong to vanish, but
    ! we must use e = e' + Ulong for the first term.

    do i = 1, nfnc
       t(i) = fourpi * deltar * sum((0.5*(c(:,i,i1) + e(:,i,i1)) &
            & * (e(:,i,i1) + ulong(:,i)) - c(:,i,i1)) * r(:)**2)
    end do

    ! The excess chemical potential of the ith component is then sum_j
    ! rho_j t_ij

    muex = 0.0

    do i = 1, ncomp
       do j = 1, ncomp
          if (i.le.j) then
             ij = i + j*(j-1)/2
          else
             ij = j + i*(i-1)/2
          end if
          muex(i) = muex(i) + rho(j) * (t(ij) + tl(ij))
       end do
    end do

    ! Also valid ONLY for HNC is the expression for the free energy
    ! density f = sum_mu rho_mu mu_mu - p (we compute the excess).

    fvex = sum(rho(:) * muex(:)) - rhotot * (cf_mf + cf_xc)
    fnex = sum(rho(:) * muex(:)) / rhotot - (cf_mf + cf_xc)

  end subroutine make_thermodynamics


  subroutine write_thermodynamics
    integer :: i

    if (model_type.eq.3 .or. model_type.eq.4) then
       print *, 'No thermodynamics for this potential type'
    else
       print *, 'Total density = ', sum(rho)
       print *, 'Compressibility factor: mean field contribution = ', cf_mf
       print *, 'Compressibility factor: contact contribution = ', cf_gc
       print *, 'Compressibility factor: correlation contribution = ', cf_xc
       print *, 'Compressibility factor: total = ', 1.0 + cf_mf + cf_gc + cf_xc
       print *, 'Pressure (virial route) = ', press
       print *, 'Excess pressure (virial route) = ', sum(rho) * (cf_mf + cf_xc)
       print *, 'Compressibility: correlation contribution = ', comp_xc
       print *, 'Compressibility: total = ', comp
       print *, 'Internal energy: mean field contribution = ', un_mf
       print *, 'Internal energy: correlation contribution = ', un_xc
       print *, 'Internal energy: un (per particle) = ', un
       print *, 'Internal energy: un / 3 = ', un / 3.0
       print *, 'Internal energy: uv (per unit volume) = ', uv
       do i = 1, ncomp
          print *, 'Chemical potential: species ', i, ' = ', muex(i)
       end do
       print *, 'Excess free energy: fvex (per unit volume) = ', fvex
       print *, 'Excess free energy: fnex (per particle) = ', fnex
    end if

  end subroutine write_thermodynamics

! Routine to reduce A.X = B using Gauss-Jordan elimination, with
! pivoting (see Numerical Recipes for a discussion of this).
!
! The input arrays are A(N, N) and B(N, M).  The output is in B(N, M),
! where X(I, :) = B(PERM(I), :).  An integer return code IRC is
! zero if successful.
!
! Note: this routine was developed independently of the gaussj routine
! in Numerical Recipes.  Differences are that we are somewhat
! profligate with bookkeeping, we don't attempt to overwrite the A
! matrix with anything useful, we provide the user with the pivot
! permutation rather than rearranging B, and we make judicious use of
! FORTRAN 90 language features.
!
! To do the pivoting, we use logical arrays to keep track of which
! rows and columns are valid in the pivot search stage, and an integer
! array PERM(:) to keep track of which column contains the pivot of
! each row.  This integer array then ends up encoding the permutation
! of the rows of B.
!
! One could move the A, B, PERM and auxiliary arrays out into global
! scope but timing tests indicate this is slower than the present
! implementation.

  subroutine axeqb_reduce(a, n, b, m, perm, irc)
    implicit none
    integer :: i, j, ii, jj, p
    integer, intent(in) :: n, m
    integer, intent(out) :: perm(n), irc
    double precision :: alpha, amax, aa
    double precision, parameter :: eps = 1D-10
    double precision, intent(inout) :: a(n, n), b(n, m)
    logical :: row(n), col(n)

    irc = 0

    ! Initially, all rows and all columns are allowed.

    row = .true.
    col = .true.

    do p = 1, n

       ! Search for a suitable pivot in the allowed rows and
       ! columns.  After we have done this p = 1...n times, we
       ! will have run out of pivots and reduced A to a permutation
       ! matrix.

       amax = 0.0
       
       do i = 1, n
          do j = 1, n
             aa = abs(a(i, j))
             if (row(i) .and. col(j) .and. (amax .lt. aa)) then
                amax = aa
                ii = i
                jj = j
             end if
          end do
       end do

       if (amax.lt.eps) then
          irc = 1
          return
       end if

       ! Having found our next pivot, mark the corresponding row and
       ! column as no longer valid in the pivot search, and save the
       ! column that the pivot is in.

       row(ii) = .false.
       col(jj) = .false.
       perm(ii) = jj

       ! Now do the elimination -- first scale the pivot row, then
       ! eliminate the entries that correspond to the pivot column
       ! in all the non-pivot rows.  After each operation the
       ! solution X remains unchanged, but the matrix A is
       ! progressively simplified.

       alpha = 1.0 / a(ii, jj)
       a(ii, :) = alpha * a(ii, :)
       b(ii, :) = alpha * b(ii, :)

       do i = 1, n
          if (i.ne.ii) then
             alpha = a(i, jj)
             a(i, :) = a(i, :) - alpha * a(ii, :)
             b(i, :) = b(i, :) - alpha * b(ii, :)
          end if
       end do

    end do

    ! At this point A_ij = 1 if the ij-th element was chosen as a
    ! pivot, and A_ij = 0 otherwise.  Each row, and each column, of A
    ! therefore contains exactly one unit entry, thus A is a
    ! permutation matrix, also encoded in PERM(:) so that
    ! X(I, :) = B(PERM(I), :).

  end subroutine axeqb_reduce

end module wizard

! End of oz_mod.f90
