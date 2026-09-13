!> Tests for proportionality_mod and make_inputs_mod.
!>
!> The proportionality tests construct sweeps in a KNOWN regime - clean linear, nonlinearly contaminated, or round-off dominated - and require the assessment both to flag the failure and to attribute it correctly, since the corrective action (reduce vs increase eps) is opposite in the two failure cases.
!>
!> The deck tests enforce the TWIN INVARIANT: a perturbed deck and its unperturbed twin must differ in exactly one functional setting (ns.radblob).  If anything else differs, the snapshot subtraction that stands in for the NS-MFP body force B_f is no longer a linearised perturbation and the pipeline silently produces a wrong answer.
program test_proportionality
   use testing_mod
   use proportionality_mod
   use make_inputs_mod
   use dispersion_mod, only: growth_rate, PI
   use nsmfp_mod, only: fit_exponential_growth
   implicit none

   type(amplitude_run_t) :: runs(3)
   type(linearity_report_t) :: rep
   integer :: ierr
   real(dp) :: eps_rec

   call section('clean linear sweep')
   call make_sweep(runs, [1.0e-6_dp, 1.0e-5_dp, 1.0e-4_dp], 0.3_dp, &
                   0.0_dp, 0.0_dp)
   call assess_linearity(runs, rep, ierr)
   call check_int(ierr, 0, 'assessment succeeds')
   call check_int(rep%verdict, VERDICT_LINEAR, 'verdict is LINEAR')
   call check(rep%scaled_spread < 1.0e-9_dp, 'scaled curves collapse')
   call check_close_rel(rep%growth_rate_mean, 0.3_dp, 1.0e-6_dp, &
                        'mean growth rate is the planted one')
   eps_rec = recommend_amplitude(runs, rep)
   call check_close(eps_rec, 1.0e-4_dp, 1.0e-15_dp, &
                    'recommends the LARGEST linear amplitude')

   call section('growth rate is amplitude independent when linear')
   call check_close_rel(runs(1)%growth_rate, runs(3)%growth_rate, 1.0e-9_dp, &
                        'smallest and largest eps agree')

   call section('nonlinear contamination')
   ! quadratic term active; contamination grows with amplitude and with time
   call make_sweep(runs, [1.0e-6_dp, 1.0e-3_dp, 1.0e-2_dp], 0.3_dp, &
                   50.0_dp, 0.0_dp)
   call assess_linearity(runs, rep, ierr)
   call check_int(rep%verdict, VERDICT_NONLINEAR, 'verdict is NONLINEAR')
   call check(rep%late_spread > rep%early_spread, &
              'degradation is worse at LATE time')
   call check(recommend_amplitude(runs, rep) < 1.0e-6_dp, &
              'recommends REDUCING the amplitude')

   call section('round-off domination')
   ! additive noise of fixed absolute size; worst for the smallest amplitude and, because the signal grows, worst at EARLY time
   call make_sweep(runs, [1.0e-14_dp, 1.0e-13_dp, 1.0e-12_dp], 0.3_dp, &
                   0.0_dp, 1.0e-14_dp)
   call assess_linearity(runs, rep, ierr)
   call check_int(rep%verdict, VERDICT_ROUNDOFF, 'verdict is ROUND-OFF')
   call check(rep%early_spread >= rep%late_spread, &
              'degradation is worse at EARLY time')
   call check(recommend_amplitude(runs, rep) > 1.0e-12_dp, &
              'recommends INCREASING the amplitude')

   call section('input validation')
   call assess_linearity(runs(1:1), rep, ierr)
   call check_int(ierr, 1, 'rejects a single amplitude')
   runs(1)%eps = -1.0_dp
   call assess_linearity(runs, rep, ierr)
   call check_int(ierr, 3, 'rejects a non-positive amplitude')

   call test_decks()

   call summary('test_proportionality')

contains

   !> Build a synthetic sweep:
   !>   |Q'|(t) = eps e^{rate t} + quad (eps e^{rate t})^2 + noise
   subroutine make_sweep(r, epss, rate, quad, noise)
      type(amplitude_run_t), intent(inout) :: r(:)
      real(dp), intent(in) :: epss(:), rate, quad, noise
      integer, parameter :: N = 120
      integer :: i, j, e
      real(dp) :: t, lin, slope, intercept, r2, ph

      do i = 1, size(r)
         if (allocated(r(i)%times)) deallocate(r(i)%times)
         if (allocated(r(i)%norms)) deallocate(r(i)%norms)
         allocate(r(i)%times(N), r(i)%norms(N))
         r(i)%eps = epss(i)
         do j = 1, N
            t = 5.0_dp * real(j-1,dp) / real(N-1,dp)
            r(i)%times(j) = t
            lin = epss(i) * exp(rate*t)
            ! deterministic pseudo-noise, different per run and sample
            ph = sin(12.9898_dp*real(j,dp) + 78.233_dp*real(i,dp)) * &
               43758.5453_dp
            ph = ph - floor(ph)
            r(i)%norms(j) = abs(lin + quad*lin*lin + noise*(2.0_dp*ph - 1.0_dp))
         end do
         call fit_exponential_growth(r(i)%times, r(i)%norms, r(i)%times(1), &
                                     -1.0_dp, slope, intercept, r2, e)
         r(i)%growth_rate = slope
         r(i)%r_squared = r2
      end do
   end subroutine make_sweep

   subroutine test_decks()
      type(run_plan_t) :: pb, pp
      character(len=:), allocatable :: tb, tp
      character(len=*), parameter :: TMPL = 'inputs/inputs_rayleigh_template.txt'
      integer :: e, ndiff
      logical :: exists

      inquire(file=TMPL, exist=exists)
      if (.not. exists) then
         call section('deck generation (SKIPPED: template not found)')
         write(*,'(a)') '  run the test suite from the project root'
         return
      end if

      call section('deck generation: the twin invariant')
      pb%kr0 = 0.7_dp; pb%eps = 0.0_dp;    pb%role = 'twin'
      pp%kr0 = 0.7_dp; pp%eps = 1.0e-3_dp; pp%role = 'perturbed'
      call plan_run(pb)
      call plan_run(pp)
      call render_deck(pb, TMPL, tb, e)
      call check_int(e, 0, 'twin deck renders')
      call render_deck(pp, TMPL, tp, e)
      call check_int(e, 0, 'perturbed deck renders')
      call check(index(tb, '@') == 0, 'no unsubstituted tokens in twin')
      call check(index(tp, '@') == 0, 'no unsubstituted tokens in perturbed')

      call count_functional_diffs(tb, tp, ndiff)
      call check_int(ndiff, 1, &
                     'exactly ONE functional line differs between twins')
      call check(has_line_starting(tb, 'ns.radblob'), 'twin sets ns.radblob')
      call check(has_line_starting(tp, 'ns.radblob'), &
                 'perturbed sets ns.radblob')

      call section('deck contents')
      call check(has_setting(tp, 'ns.probtype', '41'), 'probtype = 41')
      call check(has_setting(tp, 'ns.axis_dir', '4'), 'axis_dir = 4')
      call check(has_setting(tp, 'geometry.coord_sys', '1'), 'RZ geometry')
      call check(has_setting(tp, 'ns.visual_nddata_format', '1'), &
                 'plotfile output enabled (else nothing to post-process)')
      call check(has_setting(tp, 'amr.max_level', '0'), &
                 'single AMR level, so twins share a grid')
      call check(has_setting(tp, 'amr.LSA_activate', '0'), &
                 'built-in power-iteration LSA disabled')

      call section('deck planning')
      call check(pp%fixed_dt < pp%dt_cap, 'dt below the capillary limit')
      call check_int(mod(pp%max_step, pp%plot_int), 0, &
                     'max_step is a multiple of plot_int')
      call check(pp%n_snapshots >= 2, 'at least two snapshots')
      call check_close_rel(pp%lambda, 2.0_dp*PI/0.7_dp, 1.0e-12_dp, &
                           'box holds exactly one wavelength')
      call check_close_rel(pp%sigma_exact, &
                           growth_rate(0.7_dp,1.0_dp,1.0_dp,1.0_dp), &
                           1.0e-12_dp, 'annotated growth rate is correct')
      call check_int(mod(pp%nr, 8), 0, 'nr is block-friendly')
      call check_int(mod(pp%nz, 8), 0, 'nz is block-friendly')
      call check(pp%dr/pp%dz > 0.5_dp .and. pp%dr/pp%dz < 2.0_dp, &
                 'cells are near-square')

      call section('stable band still yields a finite run')
      pp%kr0 = 1.3_dp
      call plan_run(pp)
      call check(pp%stop_time > 0.0_dp, 'positive stop time')
      call check(pp%max_step > 0, 'positive step count')
      call check_close(pp%sigma_exact, 0.0_dp, 0.0_dp, &
                       'zero growth rate in the stable band')
   end subroutine test_decks

   !> Count differing non-comment, non-blank lines between two decks.
   subroutine count_functional_diffs(a, b, ndiff)
      character(len=*), intent(in) :: a, b
      integer, intent(out) :: ndiff
      character(len=512) :: la, lb
      integer :: pa, pb_, na, nb
      ndiff = 0
      pa = 1; pb_ = 1
      do
         call next_functional(a, pa, la, na)
         call next_functional(b, pb_, lb, nb)
         if (na == 0 .and. nb == 0) exit
         if (na /= nb) then
            ndiff = ndiff + 1
            exit
         end if
         if (trim(la) /= trim(lb)) ndiff = ndiff + 1
      end do
   end subroutine count_functional_diffs

   !> Next line that configures the solver (comments and blanks skipped).
   subroutine next_functional(s, pos, line, found)
      character(len=*), intent(in) :: s
      integer, intent(inout) :: pos
      character(len=*), intent(out) :: line
      integer, intent(out) :: found
      integer :: nl, h
      character(len=512) :: t
      found = 0
      do
         if (pos > len(s)) return
         nl = index(s(pos:), new_line('a'))
         if (nl == 0) then
            t = s(pos:)
            pos = len(s) + 1
         else
            t = s(pos:pos+nl-2)
            pos = pos + nl
         end if
         h = index(t, '#')
         if (h > 0) t = t(1:h-1)
         if (len_trim(adjustl(t)) > 0) then
            line = adjustl(t)
            found = 1
            return
         end if
      end do
   end subroutine next_functional

   logical function has_line_starting(s, key)
      character(len=*), intent(in) :: s, key
      has_line_starting = index(s, new_line('a')//key) > 0 .or. &
                          index(s, key) == 1
   end function has_line_starting

   logical function has_setting(s, key, val)
      character(len=*), intent(in) :: s, key, val
      integer :: p, nl, eq
      character(len=512) :: t
      has_setting = .false.
      p = index(s, new_line('a')//key//' ')
      if (p == 0) p = index(s, new_line('a')//key//'=')
      if (p == 0) return
      p = p + 1
      nl = index(s(p:), new_line('a'))
      if (nl == 0) return
      t = s(p:p+nl-2)
      eq = index(t, '=')
      if (eq == 0) return
      t = adjustl(t(eq+1:))
      p = index(t, '#')
      if (p > 0) t = t(1:p-1)
      has_setting = trim(adjustl(t)) == trim(val)
   end function has_setting

end program test_proportionality
