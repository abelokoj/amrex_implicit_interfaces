!> Generate matched AMReX input decks for the Rayleigh-Plateau NS-MFP study.
!>
!> Fortran counterpart of python/make_inputs.py.
!>
!> Every deck is expanded from inputs/inputs_rayleigh_template.txt so that a perturbed run and its unperturbed twin differ in exactly one token (ns.radblob).  That is the whole point: the twin subtraction stands in for the NS-MFP body force B_f, and it is only valid if nothing else differs. Hand-editing decks is how that invariant gets broken silently, so generate them instead.
module make_inputs_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use dispersion_mod, only: growth_rate, growth_rate_viscous, ohnesorge, PI
   implicit none
   private

   public :: dp
   public :: run_plan_t, plan_run, render_deck, write_deck, describe_plan
   public :: capillary_dt

   !> Everything needed to emit one deck.
   type :: run_plan_t
      character(len=64) :: role = ''
      real(dp) :: kr0 = 0.7_dp, k = 0.0_dp, eps = 0.0_dp, lambda = 0.0_dp
      real(dp) :: r0 = 1.0_dp, sigma = 1.0_dp
      real(dp) :: rho_l = 1.0_dp, rho_g = 0.001225_dp
      real(dp) :: mu_l = 0.02_dp, mu_g = 0.00026_dp
      real(dp) :: rmax = 4.0_dp
      integer  :: cells_per_r0 = 32
      real(dp) :: safety = 0.25_dp, n_periods = 5.0_dp
      integer  :: target_snapshots = 300
      integer  :: nr = 0, nz = 0
      real(dp) :: dr = 0.0_dp, dz = 0.0_dp
      real(dp) :: fixed_dt = 0.0_dp, dt_cap = 0.0_dp
      real(dp) :: stop_time = 0.0_dp
      integer  :: max_step = 0, plot_int = 0, check_int = 0, n_snapshots = 0
      real(dp) :: sigma_exact = 0.0_dp, sigma_visc = 0.0_dp, oh = 0.0_dp
   end type run_plan_t

contains

   !> Brackbill capillary-wave time-step limit, dt < sqrt((rho_l + rho_g) dx^3 / (2 pi sigma)). Surface tension is treated implicitly by the solver, which relaxes this, but we stay under the explicit limit because the perturbation being measured is tiny and temporal error must not compete with it.
   pure function capillary_dt(dx, sigma, rho_l, rho_g) result(dt)
      real(dp), intent(in) :: dx, sigma, rho_l, rho_g
      real(dp) :: dt
      dt = sqrt((rho_l + rho_g) * dx**3 / (2.0_dp * PI * sigma))
   end function capillary_dt

   !> Round a cell count up to a multiple of `block` (at least `block`).
   pure function round_to_block(x, blk) result(n)
      real(dp), intent(in) :: x
      integer, intent(in) :: blk
      integer :: n, m
      m = ceiling(x)
      n = max(blk, ((m + blk - 1) / blk) * blk)
   end function round_to_block

   !> Work out every derived quantity for one deck.
   subroutine plan_run(p)
      type(run_plan_t), intent(inout) :: p
      real(dp) :: t_c, dxmin

      p%k = p%kr0 / p%r0
      p%lambda = 2.0_dp * PI / p%k

      p%dr = p%r0 / real(p%cells_per_r0, dp)
      p%nr = round_to_block(p%rmax / p%dr, 8)
      p%dr = p%rmax / real(p%nr, dp)
      ! aim for near-square cells, then round to a block-friendly count
      p%nz = round_to_block(p%lambda / p%dr, 8)
      p%dz = p%lambda / real(p%nz, dp)

      dxmin = min(p%dr, p%dz)
      p%dt_cap = capillary_dt(dxmin, p%sigma, p%rho_l, p%rho_g)
      p%fixed_dt = p%safety * p%dt_cap

      p%sigma_exact = growth_rate(p%k, p%r0, p%sigma, p%rho_l)
      p%sigma_visc = growth_rate_viscous(p%k, p%r0, p%sigma, p%rho_l, p%mu_l)
      p%oh = ohnesorge(p%mu_l, p%rho_l, p%sigma, p%r0)

      ! March for n_periods e-folding times of the slower (viscous) estimate. In the stable band no e-folding time exists, so fall back on the capillary time scale rather than dividing by ~0.
      t_c = sqrt(p%rho_l * p%r0**3 / p%sigma)
      if (p%sigma_visc > 1.0e-12_dp) then
         p%stop_time = p%n_periods / p%sigma_visc
      else
         p%stop_time = p%n_periods * t_c
      end if

      p%max_step = ceiling(p%stop_time / p%fixed_dt)
      p%plot_int = max(1, p%max_step / p%target_snapshots)
      ! land max_step on an exact multiple of plot_int so both twins produce the same final plotfile
      p%max_step = p%plot_int * (p%max_step / p%plot_int)
      p%check_int = max(p%plot_int * 50, p%plot_int)
      p%n_snapshots = p%max_step / p%plot_int
      p%stop_time = real(p%max_step, dp) * p%fixed_dt
   end subroutine plan_run

   !> Expand the template, substituting every @TOKEN@.
   subroutine render_deck(p, template_path, text, ierr)
      type(run_plan_t), intent(in) :: p
      character(len=*), intent(in) :: template_path
      character(len=:), allocatable, intent(out) :: text
      integer, intent(out) :: ierr

      character(len=1024) :: line
      character(len=:), allocatable :: buf
      integer :: u, ios

      ierr = 0
      open(newunit=u, file=trim(template_path), status='old', &
           action='read', iostat=ios)
      if (ios /= 0) then
         ierr = 1; return
      end if

      buf = ''
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         buf = buf // trim(line) // new_line('a')
      end do
      close(u)

      call subst(buf, '@ROLE@',         trim(p%role))
      call subst(buf, '@KR0@',          fmt_g(p%kr0))
      call subst(buf, '@LAMBDA@',       fmt_g(p%lambda))
      call subst(buf, '@EPS@',          fmt_g(p%eps))
      call subst(buf, '@SIGMA_EXACT@',  fmt_g(p%sigma_exact))
      call subst(buf, '@RMAX@',         fmt_g(p%rmax))
      call subst(buf, '@NR@',           fmt_i(p%nr))
      call subst(buf, '@NZ@',           fmt_i(p%nz))
      call subst(buf, '@DR@',           fmt_g(p%dr))
      call subst(buf, '@DT_CAP@',       fmt_g(p%dt_cap))
      call subst(buf, '@FIXED_DT@',     fmt_g(p%fixed_dt))
      call subst(buf, '@MAX_STEP@',     fmt_i(p%max_step))
      call subst(buf, '@STOP_TIME@',    fmt_g(p%stop_time))
      call subst(buf, '@PLOT_INT@',     fmt_i(p%plot_int))
      call subst(buf, '@CHECK_INT@',    fmt_i(p%check_int))
      call subst(buf, '@SIGMA@',        fmt_g(p%sigma))
      call subst(buf, '@RHO_L@',        fmt_g(p%rho_l))
      call subst(buf, '@RHO_G@',        fmt_g(p%rho_g))
      call subst(buf, '@MU_L@',         fmt_g(p%mu_l))
      call subst(buf, '@MU_G@',         fmt_g(p%mu_g))

      if (index(buf, '@') > 0) then
         ierr = 2                       ! unsubstituted token left behind
      end if
      text = buf
   end subroutine render_deck

   subroutine subst(s, tok, val)
      character(len=:), allocatable, intent(inout) :: s
      character(len=*), intent(in) :: tok, val
      integer :: p
      do
         p = index(s, tok)
         if (p == 0) exit
         s = s(1:p-1) // val // s(p+len(tok):)
      end do
   end subroutine subst

   function fmt_g(x) result(s)
      real(dp), intent(in) :: x
      character(len=:), allocatable :: s
      character(len=64) :: t
      write(t,'(g0.10)') x
      s = trim(adjustl(t))
   end function fmt_g

   function fmt_i(n) result(s)
      integer, intent(in) :: n
      character(len=:), allocatable :: s
      character(len=32) :: t
      write(t,'(i0)') n
      s = trim(adjustl(t))
   end function fmt_i

   subroutine write_deck(p, template_path, out_path, ierr)
      type(run_plan_t), intent(in) :: p
      character(len=*), intent(in) :: template_path, out_path
      integer, intent(out) :: ierr
      character(len=:), allocatable :: text
      integer :: u, ios

      call render_deck(p, template_path, text, ierr)
      if (ierr /= 0) return
      open(newunit=u, file=trim(out_path), status='replace', &
           action='write', iostat=ios)
      if (ios /= 0) then
         ierr = 3; return
      end if
      write(u,'(a)', advance='no') text
      close(u)
   end subroutine write_deck

   subroutine describe_plan(p, unit)
      type(run_plan_t), intent(in) :: p
      integer, intent(in), optional :: unit
      integer :: u
      u = 6
      if (present(unit)) u = unit
      write(u,'(a,f6.3,a,f10.4,a,es10.3)') '  k r0=', p%kr0, &
         '  lambda=', p%lambda, '  eps=', p%eps
      write(u,'(a,i0,a,i0,a,es10.3,a,es10.3)') '    grid ', p%nr, 'x', p%nz, &
         '  dr=', p%dr, ' dz=', p%dz
      write(u,'(a,es10.3,a,es10.3,a,i0,a,i0)') '    dt=', p%fixed_dt, &
         ' (cap limit ', p%dt_cap, ')  steps=', p%max_step, &
         '  snapshots=', p%n_snapshots
      write(u,'(a,f10.6,a,f10.6,a,es10.3)') '    sigma_exact=', &
         p%sigma_exact, '  sigma_visc=', p%sigma_visc, '  Oh=', p%oh
   end subroutine describe_plan

end module make_inputs_mod
