!> Proportionality (linearity) test for the NS-MFP subspace.
!>
!> Fortran counterpart of python/proportionality.py.
!>
!> Why this is load-bearing rather than optional verification
!! ----------------------------------------------------------
!! NS-MFP marches the NONLINEAR solver.  Its output is a Jacobian-vector product only insofar as the O(eps^2) term in
!!
!! dQ'/dt = Q' dF/dQ + Q'^2 d2F/dQ2 + H.O.T.        (Ranjan Eq. 7)
!!
!! is negligible.  Ranjan et al. state this plainly: "linearity is ensured through simple tests based on proportionality between input and output ... Although the value of eps_0 should be evaluated for each problem". Nothing in the algorithm detects a violation: with eps too large the DMD still returns a clean-looking spectrum, but of a contaminated operator.
!!
!! This test is therefore the evidence that the decomposed operator is the linearised one.  It also distinguishes this approach from methods that infer snapshots from nonlinear simulation data, where no such proportionality is available.
!!
!! The test
!! --------
!! Run the same perturbation shape at several amplitudes eps_i.  Under linearity  Q'(t; eps_i) = (eps_i/eps_j) Q'(t; eps_j), so |Q'|/eps must collapse onto one curve and the growth rate must be eps-independent.
!!
!! Failure signatures differ, and so do the fixes: eps too LARGE  -> curves fan out at LATE time (nonlinearity)  -> reduce eps too SMALL  -> curves degrade at EARLY time (round-off)    -> increase
module proportionality_mod
   use table_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private

   public :: dp
   public :: amplitude_run_t, linearity_report_t
   public :: assess_linearity, recommend_amplitude, print_report

   integer, parameter, public :: VERDICT_LINEAR = 0
   integer, parameter, public :: VERDICT_NONLINEAR = 1
   integer, parameter, public :: VERDICT_ROUNDOFF = 2

   !> One member of the amplitude sweep.
   type :: amplitude_run_t
      real(dp) :: eps = 0.0_dp
      real(dp) :: growth_rate = 0.0_dp
      real(dp) :: r_squared = 0.0_dp
      real(dp), allocatable :: times(:)
      real(dp), allocatable :: norms(:)
   end type amplitude_run_t

   type :: linearity_report_t
      integer  :: verdict = VERDICT_LINEAR
      real(dp) :: scaled_spread = 0.0_dp
      real(dp) :: growth_rate_spread = 0.0_dp
      real(dp) :: growth_rate_mean = 0.0_dp
      real(dp) :: early_spread = 0.0_dp
      real(dp) :: late_spread = 0.0_dp
      character(len=512) :: detail = ''
   end type linearity_report_t

contains

   !> Judge whether an amplitude sweep is in the linear regime.
   !>
   !> @param runs        the sweep (at least two amplitudes) @param scaled_tol  allowed relative spread of |Q'|/eps (default 0.02) @param rate_tol    allowed growth-rate spread relative to the mean @param early_frac  fraction of the record treated as "early"
   subroutine assess_linearity(runs, rep, ierr, scaled_tol, rate_tol, &
                               early_frac)
      type(amplitude_run_t), intent(in) :: runs(:)
      type(linearity_report_t), intent(out) :: rep
      integer, intent(out) :: ierr
      real(dp), intent(in), optional :: scaled_tol, rate_tol, early_frac

      integer :: nr, n, i, j, i_early, cnt_e, cnt_l
      real(dp) :: stol, rtol, efrac
      real(dp) :: ref, dev, rel, amax, amin, rmean
      real(dp) :: sum_e, sum_l
      real(dp), allocatable :: scaled(:,:)
      logical :: scaled_ok, rate_ok

      ierr = 0
      nr = size(runs)
      if (nr < 2) then
         ierr = 1; return                   ! need >= 2 amplitudes
      end if

      stol = 0.02_dp;  if (present(scaled_tol)) stol = scaled_tol
      rtol = 0.02_dp;  if (present(rate_tol))   rtol = rate_tol
      efrac = 0.25_dp; if (present(early_frac)) efrac = early_frac

      ! common time window
      n = huge(1)
      do i = 1, nr
         if (.not. allocated(runs(i)%norms)) then
            ierr = 2; return
         end if
         n = min(n, size(runs(i)%norms))
         if (runs(i)%eps <= 0.0_dp) then
            ierr = 3; return                ! amplitude must be positive
         end if
      end do
      if (n < 2) then
         ierr = 4; return
      end if

      allocate(scaled(nr, n))
      do i = 1, nr
         scaled(i, :) = runs(i)%norms(1:n) / runs(i)%eps
      end do

      i_early = max(1, int(efrac * real(n, dp)))

      rep%scaled_spread = 0.0_dp
      sum_e = 0.0_dp; sum_l = 0.0_dp
      cnt_e = 0; cnt_l = 0
      do j = 1, n
         ref = sum(scaled(:, j)) / real(nr, dp)
         do i = 1, nr
            if (ref > 0.0_dp) then
               dev = abs(scaled(i, j) - ref)
               rel = dev / ref
            else
               rel = 0.0_dp
            end if
            rep%scaled_spread = max(rep%scaled_spread, rel)
            ! Attribution uses the window MEAN, not the max: the early and late windows hold different numbers of samples, and the maximum of a noisy sequence grows with sample count, so a max-vs-max comparison would systematically favour the longer window.
            if (j <= i_early) then
               sum_e = sum_e + rel; cnt_e = cnt_e + 1
            else
               sum_l = sum_l + rel; cnt_l = cnt_l + 1
            end if
         end do
      end do
      rep%early_spread = 0.0_dp
      if (cnt_e > 0) rep%early_spread = sum_e / real(cnt_e, dp)
      rep%late_spread = 0.0_dp
      if (cnt_l > 0) rep%late_spread = sum_l / real(cnt_l, dp)

      amax = maxval(runs(:)%growth_rate)
      amin = minval(runs(:)%growth_rate)
      rmean = sum(runs(:)%growth_rate) / real(nr, dp)
      rep%growth_rate_spread = amax - amin
      rep%growth_rate_mean = rmean

      scaled_ok = rep%scaled_spread <= stol
      rate_ok = (abs(rmean) > 0.0_dp) .and. &
                (rep%growth_rate_spread <= rtol * abs(rmean))

      if (scaled_ok .and. rate_ok) then
         rep%verdict = VERDICT_LINEAR
         rep%detail = 'Q'' scales proportionally with eps and the growth '// &
            'rate is amplitude-independent: the snapshots represent the '// &
            'linearised operator.'
      else if (rep%late_spread > rep%early_spread) then
         rep%verdict = VERDICT_NONLINEAR
         rep%detail = 'Collapse degrades at LATE time: the O(eps^2) term '// &
            'is active.  REDUCE the perturbation amplitude and re-run.'
      else
         rep%verdict = VERDICT_ROUNDOFF
         rep%detail = 'Collapse degrades at EARLY time: the perturbation '// &
            'is near machine precision.  INCREASE the amplitude, or '// &
            'tighten solver tolerances (mac.mac_abs_tol, mg.bot_atol).'
      end if
   end subroutine assess_linearity

   !> Amplitude to use for production runs.
   !>
   !> Picks the largest amplitude that is still linear: this maximises the signal-to-round-off margin without entering the nonlinear regime.
   function recommend_amplitude(runs, rep) result(eps)
      type(amplitude_run_t), intent(in) :: runs(:)
      type(linearity_report_t), intent(in) :: rep
      real(dp) :: eps

      select case (rep%verdict)
      case (VERDICT_LINEAR)
         eps = maxval(runs(:)%eps)
      case (VERDICT_NONLINEAR)
         eps = minval(runs(:)%eps) / 10.0_dp
      case default
         eps = maxval(runs(:)%eps) * 10.0_dp
      end select
   end function recommend_amplitude

   subroutine print_report(runs, rep, unit)
      type(amplitude_run_t), intent(in) :: runs(:)
      type(linearity_report_t), intent(in) :: rep
      integer, intent(in), optional :: unit
      integer :: u, i
      character(len=16) :: v
      type(table_t) :: t

      u = 6
      if (present(unit)) u = unit

      select case (rep%verdict)
      case (VERDICT_LINEAR);    v = 'LINEAR'
      case (VERDICT_NONLINEAR); v = 'NONLINEAR'
      case default;             v = 'ROUND-OFF'
      end select

      if (u == 6) then
         call t%init([cel('eps '), cel('growth rate '), cel('r^2 '), &
            cel("|Q'|/eps final")], 'rrrr')
         do i = 1, size(runs)
            call t%row([cel(fmt_e(runs(i)%eps)), cel(fmt_f(runs(i)%growth_rate, 6)), &
               cel(fmt_f(runs(i)%r_squared, 5)), &
               cel(fmt_e(runs(i)%norms(size(runs(i)%norms))/runs(i)%eps))])
         end do
         call t%emit('Proportionality (linearity) test')

         call kv_block('', &
            [cel('scaled-norm spread'), cel('growth-rate spread'), &
               cel('mean growth rate '), cel('verdict ')], &
            [cel(fmt_pct(rep%scaled_spread, 4)), &
             cel(fmt_e(rep%growth_rate_spread)), &
             cel(fmt_e(rep%growth_rate_mean)), cel(trim(v))])
         write(u,'(a)') ''
         write(u,'(a)') trim(rep%detail)
      else
         write(u,'(a)') 'Proportionality (linearity) test'
         write(u,'(a)') repeat('-', 60)
         write(u,'(a14,a16,a11,a18)') 'eps', 'growth rate', 'r^2', &
            '|Q''|/eps final'
         do i = 1, size(runs)
            write(u,'(es14.3,es16.6,f11.5,es18.6)') runs(i)%eps, &
               runs(i)%growth_rate, runs(i)%r_squared, &
               runs(i)%norms(size(runs(i)%norms)) / runs(i)%eps
         end do
         write(u,'(a)') repeat('-', 60)
         write(u,'(a,f9.4,a)') 'scaled-norm spread : ', &
            100.0_dp*rep%scaled_spread, ' %'
         write(u,'(a,es12.4,a,es12.4,a)') 'growth-rate spread : ', &
            rep%growth_rate_spread, '  (mean ', rep%growth_rate_mean, ')'
         write(u,'(a,a)') 'verdict            : ', trim(v)
         write(u,'(a)') trim(rep%detail)
      end if
   end subroutine print_report

end module proportionality_mod
