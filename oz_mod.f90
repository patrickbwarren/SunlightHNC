! This file is part of SunlightDPD - a home for open source software
! related to the dissipative particle dynamics (DPD) simulation
! method.

! Based on an original code copyright (c) 2007 Lucian Anton.
! Modifications copyright (c) 2008, 2009 Andrey Vlasov.  Additional
! modifications copyright (c) 2009-2018 Unilever UK Central Resources
! Ltd (Registered in England & Wales, Company No 29140; Registered
! Office: Unilever House, Blackfriars, London, EC4P 4BQ, UK).  Later
! modifications copyright (c) 2020-2024 Patrick B Warren
! <patrick.warren@stfc.ac.uk> and STFC.

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

! Full documentation is found in the accompanying LaTeX document.

module wizard

  implicit none

  include "fftw3.f"
  
  integer, parameter :: dp = kind(1.0d0)

  integer, parameter :: & ! Enumeration of current closures
       & NO_CLOSURE  = 0, &
       & HNC_CLOSURE = 1, &
       & RPA_CLOSURE = 2, &
       & MSA_CLOSURE = 3, &
       & EXP_CLOSURE = 4

  integer, parameter :: & ! Enumeration of current potential models
       & NO_MODEL_TYPE = 0, &
       & DPD_GAUSSIAN_CHARGES = 1, &
       & DPD_BESSEL_CHARGES = 2, &
       & DPD_LINEAR_CHARGES = 3, &
       & DPD_SLATER_APPROX_CHARGES = 4, &
       & DPD_SLATER_EXACT_CHARGES = 5, &
       & URPM_WITHOUT_USHORT = 6, &
       & URPM_WITH_USHORT = 7, &
       & RPM_WITHOUT_USHORT = 8, &
       & RPM_WITH_USHORT = 9, &
       & HARD_SPHERES = 10

  integer, parameter :: & ! Enumeration of error codes
       & NO_ERROR = 0, &
       & CONVERGENCE_ERROR = 1, &
       & AXEQB_ERROR = 2, &
       & DSYSV_ERROR = 3, &
       & DGEEV_ERROR = 4, &
       & MISMATCH_ERROR = 5

  character (len=5)  :: version = 'v1.13' ! The current version
  character (len=3)  :: closure_name = '' ! TLA for the last-used closure 
  character (len=32) :: model_name = ''   ! Model name
  character (len=48) :: error_msg = ''    ! Error message

  logical :: suppress_msgs = .false. ! In case of severe numerical instability
  logical :: silent = .false.        ! Prevent printing warning/error messages
  logical :: verbose = .false.       ! Print solver diagnostics
  logical :: cold_start = .true.     ! Force a cold start
  logical :: auto_fns = .true.       ! Calculate stuff at end

  real(kind=dp), parameter :: &
       & pi = 4.0_dp * atan(1.0_dp), &
       & twopi = 2.0_dp * pi, &
       & fourpi = 4.0_dp * pi, &
       & rootpi = sqrt(pi)

  integer :: &
       & start_type = 3,    & ! how to initialise in a cold start
       & model_type = 0,    & ! which potential was (last) chosen
       & closure_type = 0,  & ! last-used closure 
       & istep,             & ! current position in iterative solver
       & ng = 4096,         & ! grid size
       & ncomp = 1,         & ! number of chemical components
       & nfnc = 0,          & ! number of functions, = ncomp (ncomp + 1) / 2
       & nps = 6,           & ! number of previous states used in Ng method
       & npic = 6,          & ! number of Picard steps
       & maxsteps = 100,    & ! max number of steps to take for convergence
       & return_code = 0      ! error code (see above)

  integer*8 :: plan  ! FFTW plan for fast discrete sine transforms

  real(kind=dp) :: &
       & deltar = 0.01_dp,  & ! real space grid spacing
       & deltak,            & ! reciprocal space grid spacing (computed)
       & error,             & ! difference between current and previous solns
       & alpha = 0.2_dp,    & ! Picard method, fraction of new solution
       & tol = 1.0E-12_dp,  & ! Error tolerance for claiming convergence
       & rc = 1.0_dp,       & ! short-range DPD repulsion range
       & lb = 0.0_dp,       & ! long-range Coulomb coupling length
       & sigma = 1.0_dp,    & ! see below (*)
       & sigmap = 1.0_dp,   & ! +- long-range Coulomb smearing length (URPM)
       & kappa = -1.0_dp,   & ! +- long-range Coulomb smoothing parameter (RPM)
       & rgroot = 1.0_dp,   & ! linear charge smearing range (Groot)
       & lbda = 1.0_dp,     & ! Slater charge smearing range (exact)
       & beta = 1.0_dp,     & ! Slater charge smearing range (approx)
       & cf_mf, cf_xc,      & ! the virial route pressure contributions ..
       & cf_gc, pex, press, & ! .. and the virial route pressure
       & comp, comp_xc,     & ! compressibility, and excess
       & aex_rl, aex_rs,    & ! contributions to HNC excess free energy density
       & aex_kl, aex_kd,    & !    --- '' ---
       & aex_ks, aex,       & ! -- '' --, and final excess free energy density
       & deficit,           & ! excess free energy deficit (see docs)
       & uex, uex_mf, uex_xc  ! energy density total and contributions

  ! (*) sigma is used both for the long-range Coulomb smearing length
  ! for the soft potentials, and for hard core diameter for RPM models
  
  real(kind=dp), allocatable :: &
       & rho(:),            & ! density array
       & z(:),              & ! valence array
       & arep(:, :),        & ! repulsion amplitude array
       & dd(:),             & ! hard core diameter array
       & diam(:, :),        & ! ditto, pairwise
       & tp(:),             & ! mean field pressure contribution
       & tu(:),             & ! mean field energy contribution
       & tl(:),             & ! mean field long range potential
       & muex(:),           & ! chemical potential array
       & c(:, :, :),        & ! direct correlation functions (dcfs)
       & e(:, :, :),        & ! indirect correlation functions (icfs)
       & h0(:, :),          & ! reference total correlation functions
       & hr(:, :, :),       & ! current total correlation functions
       & hc(:, :),          & ! current contact values
       & ck(:, :),          & ! transform of dcfs
       & ek(:, :),          & ! transform of icfs
       & sk(:, :, :),       & ! partial structure factors
       & ushort(:, :),      & ! short range potential in real space
       & expnegus(:, :),    & ! exp(-ushort) (includes hard cores)
       & dushort(:, :),     & ! derivative of the same
       & ulong(:, :),       & ! long range potential in real space
       & dulong(:, :),      & ! derivative of the same
       & ulongk(:, :),      & ! long range potential in reciprocal space
       & u0(:),             & ! r = 0 of long range potential (diagonal)
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
    allocate(diam(ncomp, ncomp))
    allocate(dd(nfnc))
    allocate(tp(nfnc))
    allocate(tu(nfnc))
    allocate(tl(nfnc))
    allocate(c(ng-1, nfnc, nps))
    allocate(e(ng-1, nfnc, nps))
    allocate(h0(ng-1, nfnc))
    allocate(hr(ng-1, ncomp, ncomp))
    allocate(hc(ncomp, ncomp))
    allocate(ck(ng-1, nfnc))
    allocate(ek(ng-1, nfnc))
    allocate(sk(ng-1, ncomp, ncomp))
    allocate(ushort(ng-1, nfnc))
    allocate(expnegus(ng-1, nfnc))
    allocate(dushort(ng-1, nfnc))
    allocate(ulong(ng-1, nfnc))
    allocate(dulong(ng-1, nfnc))
    allocate(ulongk(ng-1, nfnc))
    allocate(u0(ncomp))
    allocate(r(ng-1))
    allocate(k(ng-1))
    allocate(fftwx(ng-1))
    allocate(fftwy(ng-1))

    ! Default values

    rho = 0.0_dp
    arep = 0.0_dp
    diam = 0.0_dp
    z = 0.0_dp
    dd = 0.0_dp
    h0 = 0.0_dp

    ! Make grids

    deltak = pi / (real(ng, kind=dp) * deltar)

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
    print *, 'POTENTIAL DETAILS :: model type =', model_type, model_name
    if (model_type.eq.NO_MODEL_TYPE) then
       print *, 'No potential has been selected'
    else if (model_type.le.DPD_SLATER_EXACT_CHARGES) then
       print *, 'DPD potential was selected, matrix A = '
       do i = 1, ncomp
          print *, ' ', arep(i, :)
       end do
       print *, ' valencies, z = ', z
       print *, ' rc = ', rc, ' lb = ', lb
       if (model_type.eq.DPD_GAUSSIAN_CHARGES) then
          print *, ' Gaussian smearing, sigma = ', sigma
       end if
       if (model_type.eq.DPD_BESSEL_CHARGES) then
          print *, ' Bessel smearing, sigma = ', sigma
       end if
       if (model_type.eq.DPD_LINEAR_CHARGES) then
          print *, ' linear (Groot) smearing, R = ', rgroot
          print *, ' equivalent Gaussian sigma = ', &
               & sqrt(2.0_dp/15.0_dp) * rgroot
       end if
       if (model_type.eq.DPD_SLATER_APPROX_CHARGES) then
          print *, ' Slater smearing (approx), beta = ', beta
          print *, ' 1 / beta = ', 1.0_dp / beta, ', &
               & 5 / (8 beta) = ', 0.625_dp / beta
          print *, ' equivalent Gaussian sigma = ', sqrt(0.5_dp) / beta
       end if
       if (model_type.eq.DPD_SLATER_EXACT_CHARGES) then
          print *, ' Slater smearing (exact), lambda = ', lbda
          print *, ' equivalent Gaussian sigma = ', lbda
       end if
    else if (model_type.le.URPM_WITH_USHORT) then
       if (model_type.eq.URPM_WITHOUT_USHORT) then
          print *, 'URPM potential was selected with ushort unused'
       else
          print *, 'URPM potential was selected with ushort used'
       end if
       print *, ' lb = ', lb, ' sigma = ', sigma, ' sigmap = ', sigmap
    else if (model_type.le.RPM_WITH_USHORT) then
       if (model_type.eq.RPM_WITHOUT_USHORT) then
          print *, 'RPM potential was selected with ushort unused'
       else
          print *, 'RPM potential was selected with ushort used'
       end if
       if (kappa.lt.0.0_dp) then
          print *, ' lb = ', lb, ' sigma = ', sigma, ' kappa -> infinity'
       else
          print *, ' lb = ', lb, ' sigma = ', sigma, ' kappa = ', kappa
       end if
       print *, ' valencies, z = ', z
       print *, 'hard core diameter array'
       do i = 1, ncomp
          print *, ' ', diam(i, :)
       end do
       print *, ' hard sphere diameters, dd = ', dd
    else if (model_type.eq.HARD_SPHERES) then
       print *, 'Hard sphere potential with sigma = ', sigma
       print *, 'hard core diameter array'
       do i = 1, ncomp
          print *, ' ', diam(i, :)
       end do
       print *, ' hard sphere diameters, dd = ', dd
    else
       print *, 'Undefined potential'
    end if
    print *, 'SYSTEM DETAILS'
    print *, ' rho = ', rho
    if (sum(rho).gt.0.0_dp) print *, ' x = ', rho(:) / sum(rho)
    print *, ' sum(rho) = ', sum(rho)
    print *, ' sum(rho*z) = ', sum(rho(:)*z(:))
    print *, '====================================================='
  end subroutine write_params

! Build the potential arrays, with parameters rc and arep(:,:) for the
! short-range DPD repulsion, and lb, sigma and z(:) for the long-range
! Gaussian-smeared Coulomb part.  A factor beta = 1/kT is implicit in
! these definitions.  The parameter charge_type is (default) 1 =
! Gaussian, 2 = Bessel, 3 = Groot (linear), 4 = Slater (approx), 5 =
! Slater (exact).  But use the defined integer constants for these!

  subroutine dpd_potential(charge_type)
    implicit none
    integer, intent(in), optional :: charge_type
    integer :: ctype = DPD_GAUSSIAN_CHARGES
    integer :: i, j, ij, irc
    real(kind=dp) :: ulong0, aa(nfnc), zz(nfnc)

    if (present(charge_type) .and. charge_type.gt.0) ctype = charge_type

    model_type = ctype

    ! Sort out some recoded potential parameters.  Also force the
    ! amplitude matrix to be symmetric, set by upper triangular
    ! entries.  This enforces the rule in the documentation and also
    ! simplifies the printing.

    do j = 1, ncomp
       do i = 1, j
          ij = i + j*(j-1)/2
          aa(ij) = arep(i, j)
          zz(ij) = z(i) * z(j)
          if (i.lt.j) arep(j, i) = arep(i, j)
       end do
    end do

    irc = nint(rc/deltar)

    ! Leave out the amplitude, then the function can be re-used
    ! (see below)

    ushort(:,nfnc) = 0.0_dp
    ushort(1:irc,nfnc) = 0.5_dp * (1.0_dp - r(1:irc)/rc)**2

    dushort(:,nfnc) = 0.0_dp
    dushort(1:irc,nfnc) = - (1.0_dp - r(1:irc)/rc) / rc

    ! Gaussian charges

    if (ctype .eq. DPD_GAUSSIAN_CHARGES) then

       ulong(:,nfnc) = (lb / r) * erf(0.5_dp*r/sigma)

       ulong0 = lb / (sigma * rootpi)

       ulongk(:,nfnc) = (fourpi * lb / k**2) * exp(-k**2*sigma**2)

       dulong(:,nfnc) = lb * exp(-0.25_dp*r**2/sigma**2) &
            & / (rootpi * r * sigma) &
            & - lb * erf(0.5_dp*r/sigma) / r**2

       model_name = "DPD_with_Gaussian_charges"
       
    end if

    ! Bessel charges

    if (ctype .eq. DPD_BESSEL_CHARGES) then

       ulong(:,nfnc) = (lb / r) * (1.0_dp - exp(-r/sigma))

       ulong0 = lb / sigma

       ulongk(:,nfnc) = (fourpi * lb / k**2) &
            & * 1.0_dp / (1.0_dp + k**2*sigma**2)

       dulong(:,nfnc) = - (lb / r**2) * (1.0_dp - exp(-r/sigma) &
            & * (1.0_dp + r / sigma))
       
       model_name = "DPD_with_Bessel_charges"

    end if

    ! Linear charge smearing as in Groot [JCP v118, 11265 (2003)].
    ! Note we do not give the real space part here hence the
    ! thermodynamic calculations will be wrong.  This could be fixed
    ! up by doing a numerical FFT of the potential.  Best would be to
    ! separate off the long range 4 pi lb / k^2 first.  TO BE DONE!

    if (ctype .eq. DPD_LINEAR_CHARGES) then

       ulong(:,nfnc) = 0.0_dp

       ulong0 = 0.0_dp

       ulongk(:,nfnc) = (fourpi * lb / k**2) * 144.0_dp &
            & * (2.0_dp - 2.0_dp*cos(k*rgroot) &
            &      - k*rgroot*sin(k*rgroot))**2 &
            &                  / (k**8 * rgroot**8)

       dulong(:,nfnc) = 0.0_dp
   
       model_name = "DPD_with_linear_charges"

    end if

    ! Slater charge smearing as in Gonzales-Melchor et al, [JCP v125,
    ! 224107 (2006)] but here with exact expression for interaction.

    if (ctype .eq. DPD_SLATER_EXACT_CHARGES) then

       ulong(:,nfnc) = (lb / r) * (1.0_dp - exp(-2.0_dp*r/lbda) &
            & * (1.0_dp + 1.375_dp*r/lbda + 0.75_dp*r**2/lbda**2 &
            &   + r**3/(6.0_dp*lbda**3)) )

       ulong0 = 0.625_dp * lb / lbda

       ulongk(:,nfnc) = (fourpi * lb / k**2) * &
            & 1.0_dp / (1.0_dp + k**2*lbda**2/4.0_dp)**4

       dulong(:,nfnc) = - (lb / r**2) * (1.0_dp - exp(-2.0_dp*r/sigma) &
            & * (1.0_dp + 2.0_dp*r/lbda + 2.0_dp*r**2/lbda**2 &
            &     + 7.0_dp*r**3/(6.0_dp*lbda**3) + r**4/(3.0_dp*lbda**4)) )

       model_name = "DPD_with_Slater_(exact)_charges"

    end if

    ! Slater charge smearing as in Gonzales-Melchor et al, [JCP v125,
    ! 224107 (2006)] with original approximate expression.

    if (ctype .eq. DPD_SLATER_APPROX_CHARGES) then

       ulong(:,nfnc) = (lb / r) * (1.0_dp - exp(-2*beta*r) * &
            & (1.0_dp + beta*r) )

       ulong0 = beta * lb

       ulongk(:,nfnc) = (fourpi * lb / k**2) * &
            & 1.0_dp / (1.0_dp + k**2/(4.0_dp*beta**2))**2

       dulong(:,nfnc) = - (lb / r**2) * (1.0_dp - exp(-2.0_dp*beta*r) &
            & * (1.0_dp + 2.0_dp*beta*r*(1.0_dp + beta*r)) )

       model_name = "DPD_with_Slater_(approx)_charges"

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

    tp = pi * rc**3 * aa / 30.0_dp
    tu = tp
    tl = 0.0_dp

    ! The r = 0 diagonal parts of the long range potential
    
    u0(:) = z(:)**2 * ulong0
    
    ! Generate auxiliary function

    expnegus = exp(-ushort)

  end subroutine dpd_potential

! Build the potential arrays for the softened URPM (Gaussian charges),
! with parameters lb, sigma and sigmap.  This expects ncomp = 2, and
! will set z(1) = 1, z(2) = -1.  The logical parameter controls
! whether ushort is used or not.

  subroutine urpm_potential(use_ushort)
    implicit none
    logical, intent(in), optional :: use_ushort
    logical :: uuflag = .false.

    if (present(use_ushort)) uuflag = use_ushort

    if (ncomp.ne.2) then
       return_code = MISMATCH_ERROR
       error_msg = 'mismatch ncomp <> 2 in urpm_potential'
       if (.not.silent) print *, '** error: ', error_msg
       return
    end if

    z(1) = 1.0_dp
    z(2) = -1.0_dp

    ulong(:,1) = lb * erf(0.5_dp*r/sigma) / r

    ulongk(:,1) = fourpi * lb * exp(-k**2*sigma**2) / k**2

    dulong(:,1) = lb * exp(-0.25_dp*r**2/sigma**2) / (rootpi * r * sigma) &
         & - lb * erf(0.5_dp*r/sigma) / r**2

    ulong(:,2) = - lb * erf(0.5_dp*r/sigmap) / r

    ulongk(:,2) = - fourpi * lb * exp(-k**2*sigmap**2) / k**2

    dulong(:,2) = - lb * exp(-0.25_dp*r**2/sigmap**2) / (rootpi * r * sigmap) &
         & + lb * erf(0.5_dp*r/sigmap) / r**2

    ulong(:,3) = ulong(:,1)
    ulongk(:,3) = ulongk(:,1)
    dulong(:,3) = dulong(:,1)

    ushort(:,:) = 0.0_dp
    dushort(:,:) = 0.0_dp

    if (uuflag) then
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

    tp(1) = 0.0_dp
    tp(2) = twopi*lb*(sigmap**2 - sigma**2)
    tp(3) = 0.0_dp

    tu = tp

    ! If not using ushort, we are off the symmetry point condition and
    ! the contribution of the long range part should be incorporated
    ! into the compressibility and chemical potential expressions.

    if (.not.uuflag) then
       tl = 2.0_dp*tp
    else
       tl = 0.0_dp
    end if

    ! The r = 0 diagonal parts of the long range potential
    
    u0(:) = z(:)**2 * lb / (sigma * rootpi)
    
    ! Generate auxiliary function

    expnegus = exp(-ushort)

    if (.not.uuflag) then
       model_type = URPM_WITHOUT_USHORT
       model_name = 'URPM_without_U_short'
    else
       model_type = URPM_WITH_USHORT
       model_name = 'URPM_with_U_short'
    end if

  end subroutine urpm_potential

! Build the potential arrays for the softened RPM (charged hard
! spheres) with parameters lb, sigma and kappa.  This expects ncomp =
! 2 or 3, and will set z(1) = 1, z(2) = -1, and z(3) = 0 if required,
! and hard core diameters to sigma.  The parameter (0 or 1) controls
! whether ushort is used or not.  Using kappa < 0 implies the pure RPM
! case (kappa -> infinity).  The case ncomp = 3 corresponds to the RPM
! in the presence of a neutral hard sphere solvent.

  subroutine rpm_potential(use_ushort)
    implicit none
    integer :: i, j, ij, irc
    logical, intent(in), optional :: use_ushort
    logical :: uuflag = .false.

    if (present(use_ushort)) uuflag = use_ushort

    if (ncomp.eq.1 .or. ncomp.gt.3) then
       return_code = MISMATCH_ERROR
       error_msg = 'mismatch ncomp <> 2 or 3 in rpm_potential'
       if (.not.silent) print *, '** error: ', error_msg
       return
    end if

    z(1) = 1.0_dp
    z(2) = -1.0_dp
    if (ncomp.eq.3) z(3) = 0.0_dp

    ! Set the diameters from an array, or sigma if none are preset.
    ! The remaining logic means sigma, if not set, is set to the
    ! smallest diameter

    if (any(diam.gt.0.0_dp)) then
       do j = 1, ncomp
          do i = 1, j
             ij = i + j*(j-1)/2
             dd(ij) = diam(i, j)
          end do
       end do
       if (sigma.eq.0.0_dp) sigma = minval(dd)
    else
       dd = sigma
    end if

    ! Zero everything to start with, as for hard spheres

    ulong = 0.0_dp
    ulongk = 0.0_dp
    dulong = 0.0_dp
    ushort = 0.0_dp
    dushort = 0.0_dp

    irc = nint(sigma / deltar)

    ulong(:,1) = lb / r
    ulong(1:irc,1) = lb / sigma
    ulongk(:,1) = fourpi * lb * sin(k*sigma)  / (sigma * k**3)
    dulong(:,1) = - lb / r**2
    dulong(1:irc,1) = 0.0_dp

    if (kappa.gt.0.0_dp) then
       ulong(:,2) = - lb * erf(kappa*r) / r
       ulongk(:,2) = - fourpi * lb * exp(-k**2/(4.0_dp*kappa**2)) / k**2
       dulong(:,2) = - 2.0_dp*kappa*lb * exp(-kappa**2*r**2) / (rootpi * r) &
            & + lb * erf(kappa*r) / r**2
    else
       ulong(:,2) = - ulong(:,1)
       ulongk(:,2) = - ulongk(:,1)
       dulong(:,2) = - dulong(:,1)
    end if

    ulong(:,3) = ulong(:,1)
    ulongk(:,3) = ulongk(:,1)
    dulong(:,3) = dulong(:,1)

    ushort(:,:) = 0.0_dp
    dushort(:,:) = 0.0_dp

    if (uuflag) then
       ushort(:,2) = ulong(:,2) + ulong(:,1)
       dushort(:,2) = dulong(:,2) + dulong(:,1)
       ulong(:,2) = - ulong(:,1)
       ulongk(:,2) = - ulongk(:,1)
       dulong(:,2) = - dulong(:,1)
    end if

    ! Generate auxiliary function - this is where HS diam implemented

    expnegus = exp(-ushort)

    do ij = 1, nfnc
       irc = nint(dd(ij) / deltar)
       expnegus(1:irc, ij) = 0.0_dp
    end do

    ! These are the analytic contributions to the thermodynamics.

    tp = 0.0_dp
    tu = 0.0_dp
    tl = 0.0_dp

    if (kappa.gt.0) then
       tp(2) = pi*lb * ( sigma * exp(-kappa**2*sigma**2) / (kappa*rootpi) &
            & + (1/(2.0_dp*kappa**2) - sigma**2/3.0_dp) * erfc(kappa*sigma) )
       tu(2) = pi*lb * ( sigma * exp(-kappa**2*sigma**2) / (kappa*rootpi) &
            & + (1/(2.0_dp*kappa**2) - sigma**2) * erfc(kappa*sigma) )
       if (.not.uuflag) then ! off SYM condition;
          ! missing last term added in v1.11
          tl(2) = pi*lb / kappa**2 - 2.0_dp*pi*lb*sigma**2 / 3.0_dp
       end if
    end if

    ! The r = 0 diagonal parts of the long range potential
    
    u0(:) = z(:)**2 * lb / sigma
    
    if (.not.uuflag) then
       model_type = RPM_WITHOUT_USHORT
       model_name = 'RPM_without_U_short'
    else
       model_type = RPM_WITH_USHORT
       model_name = 'RPM_with_U_short'
    end if

  end subroutine rpm_potential

! Build the potential arrays for hard spheres with diameter sigma.
! This works for arbitrary numbers of components.

  subroutine hs_potential
    implicit none
    integer :: i, j, ij, irc

    ! Set the hard core diameters from an array, or sigma if none are pre-set

    if (any(diam.gt.0.0_dp)) then
       do j = 1, ncomp
          do i = 1, j
             ij = i + j*(j-1)/2
             dd(ij) = diam(i, j)
          end do
       end do
    else
       dd = sigma
    end if

    ulong = 0.0_dp
    ulongk = 0.0_dp
    dulong = 0.0_dp
    ushort = 0.0_dp
    dushort = 0.0_dp

    ! This is where the hard sphere diameters enter

    expnegus = 1.0_dp

    do ij = 1, nfnc
       irc = nint(dd(ij) / deltar)
       expnegus(1:irc, ij) = 0.0_dp
    end do

    ! These are the analytic contributions to the thermodynamics.

    tp = 0.0_dp
    tu = 0.0_dp
    tl = 0.0_dp
    u0 = 0.0_dp
    
    model_type = HARD_SPHERES
    model_name = 'hard_spheres'

  end subroutine hs_potential

! The next routine solves the Ornstein-Zernicke equation to determine
! e = h - c, given c.  We re-partition the long range part of the
! potential so that the routine actually calculates c' = c + Ulong
! and e' = e - Ulong.  This is because the Fourier transform of
! Ulong can be computed in closed form.  Note h = e + c = e' + c'.

  subroutine oz_solve
    implicit none
    integer :: i1, i, j, ij, ik, irc
    integer :: perm(ncomp)
    real(kind=dp) :: &
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
            & / ( 1.0_dp - rho(1) * (ck(:, 1) - ulongk(:, 1)) ) &
            & - ck(:, 1)

    else ! Multicomponent OZ inversion

       ! First set up a unit matrix, and the diagonal R matrix -- see
       ! the documentation for the math here.

       rhomat = 0.0_dp
       unita = 0.0_dp

       do i = 1, ncomp
          rhomat(i,i) = rho(i)
          unita(i,i) = 1.0_dp
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
             return_code = AXEQB_ERROR
             error_msg = 'axeqb encountered singular problem in oz_solve'
             if (.not.suppress_msgs .and. .not.silent) then
                print *, '** error: ', error_msg
                print *, '** further messages of this kind will be suppressed'
                suppress_msgs = .true.
             end if
             return
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
! equation to determine c and e from the reference h0.  In practice as
! always we actually calculate c' = c + Ulong and e' = e - Ulong.
! Note h = e + c = e' + c'.  The result is saved to position 1 in the
! history trajectory.  Note that the offset ulongk in 

  subroutine oz_solve2
    implicit none
    integer :: i, j, ij, ik, irc
    integer :: perm(ncomp)
    real(kind=dp) :: &
         & a(ncomp, ncomp), b(ncomp, ncomp), &
         & hmat(ncomp, ncomp), rhomat(ncomp, ncomp), &
         & unita(ncomp, ncomp), hk(ng-1, nfnc)

    do i=1, nfnc
       fftwx(1:ng-1) = r(1:ng-1) * h0(1:ng-1, i)
       call dfftw_execute(plan)
       hk(1:ng-1, i) =  (twopi * deltar) * fftwy(1:ng-1) / k(1:ng-1)
    end do

    if (ncomp .eq. 1) then

       ! In the one component case the OZ solution is trivial.

       ck(:, 1) = hk(:, 1) / (1.0_dp + rho(1)*hk(:, 1))

    else ! Multicomponent OZ solution

       ! As above set up a unit matrix, and the diagonal R matrix

       rhomat = 0.0_dp
       unita = 0.0_dp

       do i = 1, ncomp
          rhomat(i,i) = rho(i)
          unita(i,i) = 1.0_dp
       end do

       ! Do the matrix calculations for each wavevector k.

       do ik = 1, ng-1

          ! Convert the reciprocal space functions into matrices.

          do j = 1, ncomp
             do i = 1, j
                ij = i + j*(j-1)/2
                hmat(i, j) = hk(ik, ij)
                if (i.lt.j) hmat(j, i) = hmat(i, j)
             end do
          end do

          ! Construct A = I + H . R, and B = H

          a = unita + matmul(hmat, rhomat)
          b = hmat

          ! Solve A.X = B so that X = (I + H.R)^(-1) . H.

          call axeqb_reduce(a, ncomp, b, ncomp, perm, irc)

          if (irc.gt.0) then
             return_code = AXEQB_ERROR
             error_msg = 'axeqb encountered singular problem in oz_solve2'
             if (.not.suppress_msgs .and. .not.silent) then
                print *, '** error: ', error_msg
                print *, '** further messages of this kind will be suppressed'
                suppress_msgs = .true.
             end if
             return
          end if

          ! Now compute C = (I + H.R)^(-1) . H + beta UL
          ! (map back to functions, and unravel the pivoting)
          ! The + beta UL is done afterwards.

          do j = 1, ncomp
             do i = 1, j
                ij = i + j*(j-1)/2
                ck(ik, ij) = b(perm(i), j)
             end do
          end do

       end do

    end if

    ! Do the Fourier back transforms.

    do i = 1, nfnc

       fftwx(1:ng-1) = k(1:ng-1) * ck(1:ng-1, i)
       call dfftw_execute(plan)
       c(1:ng-1, i, 1) =  (deltak / twopi**2) * fftwy(1:ng-1) / r(1:ng-1)

    end do
 
    ! Add + beta UL (since we have the exact real and reciprocal space
    ! functions we don't need to pass this through the DDFT).
    ! Strictly speaking this is a flourish that isn't necessary for
    ! the structural features (pair functions and structure factors)
    ! since it cancels out again.  However several of the
    ! thermodynamic calculations assume this offset is present.

    ck = ck + ulongk
    c(:,:,1) = c(:,:,1) + ulong

    ! Recover the indirect correlation functions.  This means the
    ! results can be used in the structure and thermodynamics routines
    ! as though they had come from the HNC/MSA/RPA solution.

    ek = hk - ck
    e(:, :, 1) = h0(:, :) - c(:, :, 1)

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
! potential only.  As we iterate we move forward in the history
! trajectory, cyclically.

  subroutine hnc_picard
    implicit none
    integer :: i1, i0, i
    istep = istep + 1
    i1 = mod(istep-1, nps) + 1
    i0 = i1 - 1; if (i0.eq.0) i0 = nps
    do i = 1, nfnc
       c(:,i,i1) = alpha * ( expnegus(:,i) * exp(e(:,i,i0)) &
            & - e(:,i,i0) - 1.0_dp ) &
            & + (1.0_dp - alpha) * c(:,i,i0)
    end do
  end subroutine hnc_picard

! The next routine implements the Ng method [K-C Ng, J. Chem. Phys.
! v61, 2680 (1974)] as an accelerated solver for the HNC closure.  As
! we iterate we move forward in the history trajectory, cyclically.

  subroutine hnc_ng
    implicit none
    integer :: i1, i0, j, j1, j2, p, nd, icp
    real(kind=dp) :: dc(ng-1,nfnc,nps-1), de(ng-1,nfnc,nps-1), &
         & a(nps-1,nps-1), x(nps-1), y(nps-1), yy, aux
    integer :: ipiv(nps-1), info
    integer, parameter :: lwork = 100
    real(kind=dp) :: work(lwork)
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
    a(:,:) = 0.0_dp
    x(:) = 0.0_dp
    do icp = 1, nfnc
       do j = 1, ng-1
          aux = expnegus(j,icp) * exp(e(j,icp,i0)) - 1.0_dp
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
    call DSYSV('U', nd, 1, a, nps-1, ipiv, x, nps-1, work, &
         & lwork, info)
    if (info.gt.0) then
       return_code = DSYSV_ERROR
       error_msg = 'DSYSV encountered singular problem in hnc_ng'
       if (.not.suppress_msgs .and. .not.silent) then
          print *, '** error: ', error_msg
          print *, '** further messages of this kind will be suppressed'
          suppress_msgs = .true.
       end if
       return
    end if
    do icp = 1, nfnc
       do j = 1, ng-1
          aux = e(j,icp,i0)
          do j1 = 1, nd
             aux = aux - de(j,icp,j1) * x(j1)
          end do
          c(j,icp,i1) = expnegus(j,icp) * exp(aux) - aux - 1.0_dp
       end do
    end do
  end subroutine hnc_ng

! Basic driver routine for solving HNC: take a number of Picard
! iterations to pump-prime the Ng method.  Stop when error is less
! than tolerance, or when exceed maximum number of iterations.  The
! flag cold_start indicates whether the direct correlation function
! should be re-initialised.  The initial guess to the direct
! correlation function is either zero (start_type = 1), or c = -
! Ushort (start_type = 2), or c = e^(-Ushort)-1 (start_type = 3).  Any
! of these should do in principle, but the initial convergence may be
! different.  Note from above that c is actually defined c' = c +
! Ulong, ie with the long-range part of the potential added.  Note
! that we always start from position 1 in the history trajectory, and
! the history trajectory is pump-primed by Picard steps before
! attempting the Ng accelerator.  At the end, the final solution is
! copied back to position 1 in the history trajectory.

  subroutine hnc_solve
    implicit none
    integer :: i, i1
    return_code = NO_ERROR
    istep = 1
    if (cold_start) then
       if (start_type.eq.1) c(:,:,1) = 0.0_dp
       if (start_type.eq.2) c(:,:,1) = - ushort(:,:)
       if (start_type.eq.3) c(:,:,1) = expnegus(:,:) - 1.0_dp
       cold_start = .false.
       if (verbose) then
          if (start_type.eq.1) print *, "HNC cold start c' = 0"
          if (start_type.eq.2) print *, "HNC cold start c' = -v'"
          if (start_type.eq.3) print *, "HNC cold start c' = e^(-v')-1"
       end if
    else
       if (verbose) print *, "HNC warm start c' = previous c'"
    end if
    call oz_solve
    if (return_code.gt.NO_ERROR) return
    do i = 1, maxsteps
       if (i .le. npic) then
          call hnc_picard
       else
          call hnc_ng
       end if
       call oz_solve
       if (return_code.gt.NO_ERROR) return
       call conv_test
       if (verbose) then
          if (i .le. npic) then
             print *, i, "HNC Picard, error = ", error
          else
             print *, i, "HNC     Ng, error = ", error
          end if
       end if
       if (error .lt. tol) exit
    end do
    if (error .gt. tol) then
       return_code = CONVERGENCE_ERROR
       error_msg = 'error > tol in hnc_solve'
       if (.not.silent) print *, '** warning: ', error_msg
    end if
    i1 = mod(istep-1, nps) + 1
    if (i1.ne.1) then ! copy solution to position 1
       c(:, :, 1) = c(:, :, i1)
       e(:, :, 1) = e(:, :, i1)
    end if
    closure_name = 'HNC'
    closure_type = HNC_CLOSURE
    if (auto_fns) then
       call make_pair_functions
       call make_structure_factors
       call make_thermodynamics
    end if
  end subroutine hnc_solve

! Calculate the difference between the direct correlation functions
! for the current and previous iteration, used as a convergence test;
! return answer in variable 'error'.

  subroutine conv_test
    implicit none
    integer i1, i0
!    real(kind=dp) norm
    i1 = mod(istep-1, nps) + 1
    i0 = i1 - 1; if (i0 .eq. 0) i0 = nps
    error = sqrt(deltar * sum( (c(:, :, i1) - c(:, :, i0))**2 ))
!    norm = sqrt(deltar * sum( c(:, :, i1)**2 ))
!    print *, "conv_test: norm = ", norm
  end subroutine conv_test
  
! This routine implements the MSA closure expressed as c' = - e' - 1
! within the hard core, in similar terms to the HNC closure.  Outwith
! the hard core, c' = - beta v' is left untouched, presuming it is set
! correctly in the MSA initialisation step.  As we iterate we move
! forward in the history trajectory, cyclically.
  
  subroutine msa_picard
    implicit none
    integer :: i1, i0, i, irc
    istep = istep + 1
    i1 = mod(istep-1, nps) + 1
    i0 = i1 - 1; if (i0.eq.0) i0 = nps
    do i = 1, nfnc
       irc = nint(dd(i) / deltar) ! Only work inside the hard core
       c(1:irc,i,i1) = alpha * (  - e(1:irc,i,i0) - 1.0_dp ) &
            & + (1.0_dp - alpha) * c(1:irc,i,i0)
    end do
  end subroutine msa_picard

! The next routine implements the Ng method [K-C Ng, J. Chem. Phys.
! v61, 2680 (1974)] as an accelerated solver for the MSA closure (cf
! HNC above).  As we iterate we move forward in the history
! trajectory, cyclically.
  
  subroutine msa_ng
    implicit none
    integer :: i1, i0, j, j1, j2, p, nd, icp, irc
    real(kind=dp) :: dc(ng-1,nfnc,nps-1), de(ng-1,nfnc,nps-1), &
         & a(nps-1,nps-1), x(nps-1), y(nps-1), yy, aux
    integer :: ipiv(nps-1), info
    integer, parameter :: lwork = 100
    real(kind=dp) :: work(lwork)
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
    a(:,:) = 0.0_dp
    x(:) = 0.0_dp
    do icp = 1, nfnc
       irc = nint(dd(icp) / deltar) ! Only work inside the hard core
       do j = 1, irc
          aux = - 1.0_dp
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
    call DSYSV('U', nd, 1, a, nps-1, ipiv, x, nps-1, work, &
         & lwork, info)
    if (info.gt.0) then
       return_code = DSYSV_ERROR
       error_msg = 'DSYSV encountered singular problem in msa_ng'
       if (.not.suppress_msgs .and. .not.silent) then
          print *, '** error: ', error_msg
          print *, '** further messages of this kind will be suppressed'
          suppress_msgs = .true.
       end if
       return
    end if
    do icp = 1, nfnc
       irc = nint(dd(icp) / deltar) ! Only work inside the hard core
       do j = 1, irc
          aux = e(j,icp,i0)
          do j1 = 1, nd
             aux = aux - de(j,icp,j1) * x(j1)
          end do
          c(j,icp,i1) = - aux - 1.0_dp
       end do
    end do
  end subroutine msa_ng

! Basic driver routine for solving MSA: take a number of Picard
! iterations to pump-prime the Ng method.  Stop when error is less
! than tolerance, or when exceed maximum number of iterations.  The
! flag cold_start indicates whether the direct correlation function
! should be re-initialised.  For a cold start, the initial guess to
! the direct correlation function inside the hard core is c' = -1.
! Irrespective of this, we need to make sure c' is correctly
! initialised outwith the hard core.  Note that we always start from
! position 1 in the history trajectory, and the history trajectory is
! pump-primed by Picard steps before attempting the Ng accelerator.
! At the end, the final solution is copied back to position 1 in the
! history trajectory.

  subroutine msa_solve
    implicit none
    integer :: i, i1, p, irc
    return_code = NO_ERROR
    do i = 1, nfnc
       irc = nint(dd(i) / deltar) ! Reset everywhere outwith hard core
       do p = 1, nps
          c(irc+1:ng-1,i,p) = - ushort(irc+1:ng-1,i)
       end do
    end do
    istep = 1
    if (cold_start) then
       do i = 1, nfnc
          irc = nint(dd(i) / deltar) ! Only initialise inside the hard core
          c(1:irc,i,1) = - 1.0_dp
       end do
       cold_start = .false.
       if (verbose) print *, "MSA cold start c' = -1"
    else
       if (verbose) print *, "MSA warm start c' = previous c'"
    end if
    call oz_solve
    if (return_code.gt.NO_ERROR) return
    do i = 1, maxsteps
       if (i .le. npic) then
          call msa_picard
       else
          call msa_ng
       end if
       call oz_solve
       if (return_code.gt.NO_ERROR) return
       call conv_test
       if (verbose) then
          if (i .le. npic) then
             print *, i, "MSA Picard, error = ", error
          else
             print *, i, "MSA     Ng, error = ", error
          end if
       end if
       if (error .lt. tol) exit
    end do
    if (error .gt. tol) then
       return_code = CONVERGENCE_ERROR
       error_msg = 'error > tol in msa_solve'
       if (.not.silent) print *, '** warning: ', error_msg
    end if
    i1 = mod(istep-1, nps) + 1
    if (i1.ne.1) then ! copy solution to position 1
       c(:, :, 1) = c(:, :, i1)
       e(:, :, 1) = e(:, :, i1)
    end if
    closure_name = 'MSA'
    closure_type = MSA_CLOSURE
    if (auto_fns) then
       call make_pair_functions
       call make_structure_factors
       call make_thermodynamics
    end if
  end subroutine msa_solve

! Given the HNC machinery, the implementation of the RPA is almost
! completely trivial and corresponds to one iteration through the
! Ornstein-Zernike solver given the choice c = - Ushort.  We save the
! result to position 1 in the history trajectory.
  
  subroutine rpa_solve
    implicit none
    return_code = NO_ERROR
    istep = 1
    c(:,:,1) = - ushort(:,:)
    call oz_solve
    if (return_code.gt.NO_ERROR) return
    error = 0.0_dp
    closure_name = 'RPA'
    closure_type = RPA_CLOSURE
    if (auto_fns) then
       call make_pair_functions
       call make_structure_factors
       call make_thermodynamics
    end if
  end subroutine rpa_solve

! Save the reference state, assuming the c and e functions are those
! in the position 1 in the history trajectory.  The corresponding c
! and e functions can be restored by a call to oz_solve2.

  subroutine save_reference
    implicit none
    h0(:, :) = c(:, :, 1) + e(:, :, 1)
  end subroutine save_reference

! The EXP approximation refines the current RPA/MSA solution by using
! the current solution (h = c + e) and a reference solution (h0) to
! send h0 --> (1 + h0) exp(h - h0) - 1.  A round trip through the
! alternate version of the Ornstein-Zernike relation re-calculates the
! direct and indirect correlation functions.  We assume the current
! solution is in position 1 in the history trajectory.  The reference
! state is reset to h0 = 0 after a call to this function, for safety!
  
  subroutine exp_refine
    implicit none
    real(kind=dp) :: hsave(ng-1, nfnc)
    hsave = h0
    h0 = (1.0_dp + h0) * exp(c(:,:,1) + e(:,:,1) - h0) - 1.0_dp
    call oz_solve2
    if (return_code.gt.NO_ERROR) return
    h0 = hsave
    closure_name = 'EXP'
    closure_type = EXP_CLOSURE
    if (auto_fns) then
       call make_pair_functions
       call make_structure_factors
       call make_thermodynamics
    end if
  end subroutine exp_refine

! Construct the structure factors out of the transform of the total
! correlation function.  Note that ck and ek are available after a
! call to the OZ solver.

  subroutine make_structure_factors
    implicit none
    integer :: i, j, ij
    real(kind=dp) :: hk(ng-1, nfnc)
    hk = ck + ek
    do j = 1, ncomp
       do i = 1, j
          ij = i + j*(j-1)/2
          if (i.eq.j) then
             sk(:, i, i) = rho(i) * (1.0_dp + rho(i) * hk(:, ij))
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
! for the user to implement.  We assume the c and e functions are
! those in the position 1 in the history trajectory.

  subroutine make_pair_functions
    implicit none
    integer :: i, j, ij, irc
    real(kind=dp) :: r1, r2, h1, h2

    do j = 1, ncomp
       do i = 1, j
          ij = i + j*(j-1)/2
          hr(:, i, j) = c(:, ij, 1) + e(:, ij, 1)
          hr(:, j, i) = hr(:, i, j)
       end do
    end do

    ! The contact value in the case of hard cores.  We extrapolate
    ! from the two nearest points.  The same extrapolation is done in
    ! make_thermodynamics but we want to keep these subroutines
    ! independent.

    do j = 1, ncomp
       do i = 1, j
          ij = i + j*(j-1)/2
          if (dd(ij).gt.0.0_dp) then
             irc = nint(dd(ij) / deltar)
             r1 = r(irc+1)
             r2 = r(irc+2)
             h1 = c(irc+1, ij, 1) + e(irc+1, ij, 1)
             h2 = c(irc+2, ij, 1) + e(irc+2, ij, 1)
             hc(i, j) = ( (h1 - h2) * dd(ij) + h2 * r1 - h1 * r2 ) / (r1 - r2)
          else
             hc(i, j) = 0.0_dp
          end if
          hc(j, i) = hc(i, j)
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
! 154109 (2009).  We assume the c and e functions are those in the
! position 1 in the history trajectory.
!
! The mean-field thermodynamic expressions can often be obtained
! analytically from the potential. In this routine they are calculated
! from species pair contributions, which are themselves calculated in
! the potential routines (which do not have access to the densities).

  subroutine make_thermodynamics
    implicit none
    integer :: i, j, ij, irc
    real(kind=dp) :: rhotot, r1, r2, h1, h2, hc
    real(kind=dp) :: rhoxx(nfnc), t(nfnc)

    ! rhoxx is rho x_i x_j, doubled up for the off-diagonal components

    rhotot = sum(rho)

    do j = 1, ncomp
       do i = 1, j
          ij = i + j*(j-1)/2
          if (i.eq.j) then
             rhoxx(ij) = rho(i)**2 / rhotot
          else
             rhoxx(ij) = 2.0_dp * rho(i) * rho(j) / rhotot
          end if
       end do
    end do

    ! Calculate the various contributions to the virial-route
    ! pressure.  This is the mean field contribution.

    cf_mf = sum(rhoxx(:) * tp(:))

    ! Evaluate t_ij = - (2pi/3) int_d^inf d(U_ij)/dr h_ij r^3 dr. The
    ! contribution from both end-points r = 0 (i = 0) and r = L (i =
    ! ng) vanishes, hence the sum just consists of the middle part of
    ! the trapezium rule.  With a hard core the first point
    ! contributes half.

    do i = 1, nfnc
       if (dd(i).gt.0.0_dp) then
          irc = nint(dd(i) / deltar)
          t(i) =  - twopi * deltar * sum((dushort(irc+2:,i)+dulong(irc+2:,i)) &
               & * (c(irc+2:,i,1) + e(irc+2:,i,1)) * r(irc+2:)**3) / 3.0_dp &
               & - 0.5_dp * twopi * deltar * ((dushort(irc+1,i) &
               &                                   + dulong(irc+1,i)) &
               & * (c(irc+1,i,1) + e(irc+1,i,1)) * r(irc+1)**3) / 3.0_dp
       else
          t(i) =  - twopi * deltar * sum((dushort(:,i) + dulong(:,i)) &
               & * (c(:,i,1) + e(:,i,1)) * r(:)**3) / 3.0_dp
       end if
    end do

    ! The correlation contribution is sum_ij rho x_i x_j t_ij.

    cf_xc = sum(rhoxx(:) * t(:))

    ! The contact contribution in the case of hard cores.  We
    ! extrapolate the contact value of the pair distribution function
    ! from the two nearest outside points.  The same extrapolation is
    ! done in make_pair_functions but we want to keep these
    ! subroutines independent.

    do i = 1, nfnc
       if (dd(i).gt.0.0_dp) then
          irc = nint(dd(i) / deltar)
          r1 = r(irc+1)
          r2 = r(irc+2)
          h1 = c(irc+1,i,1) + e(irc+1,i,1)
          h2 = c(irc+2,i,1) + e(irc+2,i,1)
          hc = ( (h1 - h2) * dd(i) + h2 * r1 - h1 * r2 ) / (r1 - r2)
          t(i) = twopi * dd(i)**3 * (1.0_dp + hc) / 3.0_dp
       else
          t(i) = 0.0_dp
       end if
    end do

    cf_gc = sum(rhoxx(:) * t(:))

    ! This is the final pressure.

    pex = rhotot * (cf_gc + cf_mf + cf_xc)
    press = rhotot + pex

    ! Now we do the compressibility (not to be confused with the above
    ! compressibility factor).  The long range part is accounted for
    ! separately.

    ! Evaluate t_ij = 4 pi int_0^inf c_ij r^2 dr.  Again the
    ! contribution from both endpoints vanishes, hence the sum just
    ! consists of the middle part of the trapezium rule.

    do i = 1, nfnc
       t(i) = fourpi * deltar * sum(c(:,i,1) * r(:)**2)
    end do

    ! The compressibility is 1 - sum_ij rho x_i x_j t_ij

    comp_xc = - sum(rhoxx(:) * (t(:) - tl(:)))
    comp = 1.0_dp + comp_xc

    ! Now we do the energy density.  First the mean-field
    ! contribution.

    uex_mf = rhotot * sum(rhoxx(:) * tu(:))

    ! Evaluate t_ij = 2 pi int_0^inf U_ij h_ij r^2 dr.  Note that with
    ! a hard core the contribution from the inner end-point is halved
    ! (trapezium rule) and the outer vanishes, otherwise the
    ! contribution from both ends vanishes.

    do i = 1, nfnc
       if (dd(i).gt.0.0_dp) then
          irc = nint(dd(i) / deltar)
          t(i) = twopi * deltar * sum((ushort(irc+2:,i) + ulong(irc+2:,i)) &
               & * (c(irc+2:,i,1) + e(irc+2:,i,1)) * r(irc+2:)**2) & 
               & + 0.5_dp * twopi * deltar * ((ushort(irc+1,i) &
               &                                    + ulong(irc+1,i)) &
               & * (c(irc+1,i,1) + e(irc+1,i,1)) * r(irc+1)**2)
       else
          t(i) = twopi * deltar * sum((ushort(:,i) + ulong(:,i)) &
               & * (c(:,i,1) + e(:,i,1)) * r(:)**2)
       end if
    end do

    ! The extra contribution is sum_ij rho x_i x_j t_ij

    uex_xc = rhotot * sum(rhoxx(:) * t(:))

    ! This is the final energy density.

    uex = uex_mf + uex_xc

    ! Finally do the free energy and chemical potentials for HNC

    if (closure_type.eq.HNC_CLOSURE) then

       ! Evaluate t_ij = 2 pi int_0^inf ( h_ij^2 / 2 - c_ij) r^2 dr where
       ! we use c = c' - Ulong in the second term and pull this out as a
       ! long range contribution.

       do i = 1, nfnc
          t(i) = fourpi * deltar * sum((0.5_dp*(c(:,i,1) + e(:,i,1))**2 &
               & - c(:,i,1)) * r(:)**2)
       end do
       
       aex_rs = 0.5_dp * rhotot * sum(rhoxx(:) * t(:))
       aex_rl = 0.5_dp * rhotot * sum(rhoxx(:) * tl(:))
       
       ! Call out to compute the reciprocal space contributions

       call hnc_kspace_integrals

       ! Final total excess free energy density
    
       aex = aex_rs + aex_rl + aex_kd + aex_ks + aex_kl
       
       ! Evaluate t_ij = 4 pi int_0^inf (h_ij e_ij / 2 - c_ij) r^2 dr.
       ! Note that the contribution from both end-points again vanishes.
       ! We can use c' for c in the second term because the long range
       ! part is accounted for analytically (it is the same integral as
       ! appears in the compressibility), but we must use e = e' + Ulong
       ! for the first term.

       do i = 1, nfnc
          t(i) = fourpi * deltar * sum((0.5_dp*(c(:,i,1) + e(:,i,1)) &
               & * (e(:,i,1) + ulong(:,i)) - c(:,i,1)) * r(:)**2)
       end do

       ! The excess chemical potential of the ith component is then sum_j
       ! rho_j t_ij

       muex = 0.0_dp

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

       ! The free energy density should satisfy
       ! f - sum_mu rho_mu mu_mu + p = 0.
       ! The deficit/excess here is a test of the numerics

       deficit = aex - sum(rho(:) * muex(:)) + pex

    end if

  end subroutine make_thermodynamics

  ! Calculate the reciprocal space contributions to the HNC free
  ! energy.  Note that c = c' - Ulong, and ck is available after a
  ! call to OZ solver.

  subroutine hnc_kspace_integrals
    implicit none
    integer :: i, j, ij, ik, info
    real(kind=dp) :: prefac, &
         & csubumat(ncomp, ncomp), &
         & rhomat(ncomp, ncomp), unita(ncomp, ncomp), &
         & m(ncomp, ncomp), wr(ncomp), wi(ncomp), &
         & lndet(ng-1), tr(ng-1)
    integer, parameter :: lwork = 100
    real(kind=dp) :: work(lwork)
    
    ! Calculate the short range trace fn in k-space
    
    tr = 0.0_dp   
    do i = 1, ncomp
       ij = i + i*(i-1)/2
       tr(:) = tr(:) + rho(i) * ck(:, ij)
    end do
    
    ! Calculate log(det) fn in k-space
    
    if (ncomp .eq. 1) then ! One-component case

       lndet(:) = log(1.0_dp - rho(1)*(ck(:, 1) - ulongk(:, 1)))

    else ! Multicomponent case

       ! First set up a unit matrix, and the diagonal R matrix

       rhomat = 0.0_dp
       unita = 0.0_dp

       do i = 1, ncomp
          rhomat(i, i) = rho(i)
          unita(i, i) = 1.0_dp
       end do

       ! Do the matrix calculations for each wavevector k.

       do ik = 1, ng-1

          ! Unpack the reciprocal space functions into matrices.
          ! Calculate C - beta UL first.

          do j = 1, ncomp
             do i = 1, j
                ij = i + j*(j-1)/2
                csubumat(i, j) = ck(ik, ij) - ulongk(ik, ij)
                if (i.lt.j) then
                   csubumat(j, i) = csubumat(i, j)
                end if
             end do
          end do

          ! Construct M = I - R . (C - beta UL)

          m = unita - matmul(rhomat, csubumat)

          ! Find the eigenvalues of M (which is _not_ symmetric)

          call DGEEV('N', 'N', ncomp, m, ncomp, wr, wi, &
               & m, ncomp, m, ncomp, work, lwork, info)

          if (info.gt.0) then
             return_code = DGEEV_ERROR
             error_msg = 'DGEEV failed to converge in hnc_kspace_integrals'
             if (.not.silent) then
                print *, '** error: ', error_msg
             end if
             return
          end if

          ! Calculate ln det M = 1/2 sum_i ln |lambda_i|^2

          lndet(ik) = 0.5_dp * sum(log(wr(:)**2 + wi(:)**2))

       end do ! loop over k vectors

    end if ! select single component or multicomponent case

    ! Finally the contribution is the integral of ln(det) + trace.
    ! The long range contribution in the trace is taken care of
    ! analytically as it resolves to the r = 0 value of the long range
    ! potential (see documentation).
    
    prefac = 1.0_dp / (fourpi * pi)
    
    aex_kd = prefac * deltak * sum(lndet(:) * k(:)**2)
    aex_ks = prefac * deltak * sum(tr(:) * k(:)**2)
    aex_kl = - 0.5_dp * sum(rho(:) * u0(:))

  end subroutine hnc_kspace_integrals

  subroutine write_thermodynamics
    integer :: i, j, ij

    if (closure_type.eq.NO_CLOSURE) then
       print *, 'No closure = no thermodynamics'
       return
    end if
    
    if (model_type.eq.DPD_LINEAR_CHARGES) then
       print *, 'No thermodynamics for ', model_name
       return
    end if

    print *, 'Thermodynamics for ', model_name
    print *, closure_name, ' closure, convergence = ', error
    print *, 'Compressibility factor: mean field contribution = ', cf_mf
    if (any(dd.gt.0.0_dp)) then
       print *, 'Compressibility factor: contact contribution = ', cf_gc
    end if
    print *, 'Compressibility factor: correlation contribution = ', cf_xc
    print *, 'Compressibility factor: total = ', &
         & 1.0_dp + cf_mf + cf_gc + cf_xc
    print *, 'Excess pressure (virial route) = ', pex
    print *, 'Pressure (virial route) = ', press
    print *, 'Compressibility: correlation contribution = ', comp_xc
    print *, 'Compressibility: total = ', comp
    print *, 'Internal energy density: mean field contribution = ', uex_mf
    print *, 'Internal energy density: correlation contribution = ', uex_xc
    print *, 'Internal energy density = ', uex
    print *, 'Internal energy per particle = ', uex / sum(rho(:))
    print *, 'Internal energy per particle / 3 = ', uex / (3.0_dp*sum(rho(:)))

    do j = 1, ncomp
       do i = 1, j
          ij = i + j*(j-1)/2
          if (dd(ij).gt.0.0_dp) then
             print *, 'Contact g', i, j, ' = 1 + ', hc(i, j), ' = ', 1 + hc(i, j)
          end if
       end do
    end do

    if (closure_type.eq.HNC_CLOSURE) then
       do i = 1, ncomp
          print *, 'HNC: species', i, ' muex = ', muex(i)
       end do
       print *, 'HNC: real space (short),  aex_rs = ', aex_rs
       print *, 'HNC: real space (long),   aex_rl = ', aex_rl
       print *, 'HNC: k-space, log(det),   aex_kd = ', aex_kd
       print *, 'HNC: k-space, tr (short), aex_ks = ', aex_ks
       print *, 'HNC: k-space, tr (long),  aex_kl = ', aex_kl
       print *, 'HNC: total free energy density = ', aex
       print *, 'HNC:   --ditto--  per particle = ', aex / sum(rho(:))
       print *, 'HNC: sum rho muex - pex        = ', sum(rho(:)*muex(:))-pex
       print *, 'HNC: deficit = ', deficit
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
    real(kind=dp) :: alpha, amax, aa
    real(kind=dp), parameter :: eps = 1E-10_dp
    real(kind=dp), intent(inout) :: a(n, n), b(n, m)
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

       amax = 0.0_dp
       
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

       alpha = 1.0_dp / a(ii, jj)
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
