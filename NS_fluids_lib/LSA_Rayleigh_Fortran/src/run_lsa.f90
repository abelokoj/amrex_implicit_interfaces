!> run_lsa -- NS-MFP + DMD stability analysis driver.
!>
!> Fortran counterpart of python/run_lsa.py.
!>
!> Subcommands: analyse    one (perturbed, twin) pair -> growth rate vs theory linearity  amplitude sweep -> proportionality test verdict dispersion wavenumber sweep -> measured vs exact dispersion curve converge   subspace-size convergence of the growth rate
!>
!> Examples: run_lsa analyse --perturbed runs/k0.7_eps0.001 --base runs/k0.7_base \
!>                   --kr0 0.7
!> run_lsa linearity --base runs/k0.7_base --kr0 0.7 \
!>           --case 1e-4=runs/a --case 1e-3=runs/b --case 1e-2=runs/c
!> run_lsa dispersion --root runs --eps 1e-3 --kr0 0.3,0.5,0.7,0.9
program run_lsa
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use dispersion_mod
   use table_mod
   use dmd_mod
   use nsmfp_mod
   use plotfile_mod, only: MAX_NAME
   use proportionality_mod
   implicit none

   character(len=32)  :: cmd
   character(len=512) :: pert_dir, base_dir, root, out_json
   character(len=64)  :: fields(8)
   character(len=64)  :: plt_prefix
   integer :: nfields
   real(dp) :: kr0_list(64), eps_list(64)
   integer  :: n_kr0, n_eps
   character(len=512) :: case_dirs(32)
   real(dp) :: case_eps(32)
   integer :: n_cases
   integer :: skip, stride, maxsnap, rank
   real(dp) :: tol, t_min, t_max
   real(dp) :: r0, sigma, rho_l, mu_l
   integer :: ierr

   call defaults()
   call parse_args(ierr)
   if (ierr /= 0) then
      call usage(); stop 1
   end if

   select case (trim(cmd))
   case ('analyse');    call do_analyse()
   case ('linearity');  call do_linearity()
   case ('dispersion'); call do_dispersion()
   case ('converge');   call do_converge()
   case default;        call usage(); stop 1
   end select

contains

   subroutine defaults()
      cmd = ''
      pert_dir = ''; base_dir = ''; root = 'runs'; out_json = ''
      fields(1) = 'x_velocity'; fields(2) = 'y_velocity'; nfields = 2
      plt_prefix = 'nddataPLT'
      kr0_list(1) = 0.7_dp; n_kr0 = 1
      eps_list(1) = 1.0e-3_dp; n_eps = 1
      n_cases = 0
      skip = 0; stride = 1; maxsnap = 0; rank = 0
      tol = 1.0e-10_dp; t_min = -huge(1.0_dp); t_max = -huge(1.0_dp)
      r0 = 1.0_dp; sigma = 1.0_dp; rho_l = 1.0_dp; mu_l = 0.02_dp
   end subroutine defaults

   !> Build the subspace for one twin pair and extract the growth rate.
   !>
   !> Two independent estimates are produced deliberately: rate_fit - slope of log|Q'(t)|, assumption-light, sees only the
   !>              dominant mode;
   !> rate_dmd - the leading stationary DMD eigenvalue, which also yields
   !>              the rest of the spectrum.
   !> They should agree once the transient has shed; disagreement usually means the fit window still contains the non-modal transient.
   subroutine analyse_pair(pd, bd, rate_dmd, rate_fit, r2, res, ok)
      character(len=*), intent(in) :: pd, bd
      real(dp), intent(out) :: rate_dmd, rate_fit, r2
      type(dmd_result_t), intent(out) :: res
      logical, intent(out) :: ok

      real(dp), allocatable :: X(:,:), times(:), norms(:)
      real(dp) :: dt, intercept, tlo
      integer :: e, n

      ok = .false.
      rate_dmd = 0; rate_fit = 0; r2 = 0

      call build_snapshot_matrix(pd, bd, fields(1:nfields), skip, stride, &
                                 maxsnap, X, times, dt, e, &
                                 verbose=.true., prefix=trim(plt_prefix))
      if (e /= 0) then
         write(*,'(a,i0)') 'ERROR building subspace, ierr=', e
         call explain_nsmfp_error(e)
         return
      end if

      call perturbation_norms(X, norms)
      if (maxval(abs(norms)) <= 0.0_dp) then
         write(*,'(a)') 'ERROR: the perturbation field is identically zero.'
         write(*,'(a)') 'The two runs produced identical output; check that'
         write(*,'(a)') 'the perturbed deck really has ns.radblob /= 0.'
         return
      end if

      n = size(times)
      tlo = t_min
      if (tlo <= -huge(1.0_dp)/2) tlo = times(max(1, n/3))
      call fit_exponential_growth(times, norms, tlo, t_max, &
                                  rate_fit, intercept, r2, e)
      if (e /= 0) then
         write(*,'(a)') 'WARNING: growth-rate fit failed (empty window?)'
      end if

      call dmd_compute(X, dt, rank, tol, res, e)
      if (e /= 0) then
         write(*,'(a,i0)') 'ERROR in DMD, ierr=', e
         return
      end if
      rate_dmd = leading_growth_rate(res)
      ok = .true.
   end subroutine analyse_pair

   subroutine explain_nsmfp_error(e)
      integer, intent(in) :: e
      select case (e)
      case (1); write(*,'(a)') '  no nddataPLT* directories under the perturbed run'
      case (2); write(*,'(a)') '  no nddataPLT* directories under the base run'
      case (3); write(*,'(a)') &
         '  fewer than 2 matching steps; check amr.plot_int'
      case (4); write(*,'(a)') '  too few snapshots after skip/stride'
      case (5); write(*,'(a)') &
         '  incompatible plotfiles (grid/time/domain)'
      case default
         if (e > 2000) then
            write(*,'(a,i0)') '  failed reading a BASE plotfile, code ', e-2000
         else if (e > 1000) then
            write(*,'(a,i0)') &
               '  failed reading a PERTURBED plotfile, code ', e-1000
         end if
      end select
      if (e == 12) write(*,'(a)') '  snapshot spacing is not uniform'
   end subroutine explain_nsmfp_error

   subroutine report(kr0, rate_dmd, rate_fit, r2)
      real(dp), intent(in) :: kr0, rate_dmd, rate_fit, r2
      real(dp) :: inv, vis, spread
      type(table_t) :: t, c
      inv = growth_rate(kr0/r0, r0, sigma, rho_l)
      vis = growth_rate_viscous(kr0/r0, r0, sigma, rho_l, mu_l)

      write(*,'(a)') ''
      call t%init([cel('quantity '), cel('growth rate '), cel('note ')], 'lrl')
      call t%row([cel('DMD: leading stationary mode'), cel(fmt_f(rate_dmd, 6)), &
         cel(' ')])
      call t%row([cel('log-norm slope fit '), cel(fmt_f(rate_fit, 6)), &
         cel('r^2 = '//fmt_f(r2, 5))])
      call t%row([cel('exact inviscid (Rayleigh) '), cel(fmt_f(inv, 6)), cel('theory ')])
      call t%row([cel('viscous estimate '), cel(fmt_f(vis, 6)), cel('theory ')])
      call t%emit('Rayleigh-Plateau growth rate at k r0 = '//fmt_f(kr0, 4))

      call c%init([cel('comparison '), cel('departure '), cel('verdict ')], 'lrl')
      if (inv > 0.0_dp) call c%row([cel('measured vs inviscid '), &
         cel(fmt_pct((rate_dmd-inv)/inv)), cel(' ')])
      if (vis > 0.0_dp) call c%row([cel('measured vs viscous '), &
         cel(fmt_pct((rate_dmd-vis)/vis)), cel(' ')])
      if (abs(rate_dmd) > 0.0_dp) then
         spread = abs(rate_dmd-rate_fit)/abs(rate_dmd)
         if (spread < 0.05_dp) then
            call c%row([cel('DMD vs norm fit '), cel(fmt_pct(spread)), &
               cel('consistent ')])
         else
            call c%row([cel('DMD vs norm fit '), cel(fmt_pct(spread)), &
               cel('DISAGREE: check fit window ')])
         end if
      end if
      call c%emit('')
   end subroutine report

   subroutine do_analyse()
      real(dp) :: rd, rf, r2
      type(dmd_result_t) :: res
      logical :: ok

      call analyse_pair(trim(pert_dir), trim(base_dir), rd, rf, r2, res, ok)
      if (.not. ok) stop 2
      write(*,'(a)') ''
      call dmd_summary(res, 10)
      call report(kr0_list(1), rd, rf, r2)
      if (len_trim(out_json) > 0) call save_json(out_json, kr0_list(1), rd, &
         rf, r2)
      call dmd_free(res)
   end subroutine do_analyse

   subroutine save_json(path, kr0, rd, rf, r2)
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: kr0, rd, rf, r2
      integer :: u
      open(newunit=u, file=trim(path), status='replace', action='write')
      write(u,'(a)') '{'
      write(u,'(a,f12.6,a)') '  "kr0": ', kr0, ','
      write(u,'(a,es20.10,a)') '  "growth_rate_dmd": ', rd, ','
      write(u,'(a,es20.10,a)') '  "growth_rate_fit": ', rf, ','
      write(u,'(a,es20.10,a)') '  "r_squared": ', r2, ','
      write(u,'(a,es20.10,a)') '  "growth_rate_exact_inviscid": ', &
         growth_rate(kr0/r0, r0, sigma, rho_l), ','
      write(u,'(a,es20.10)') '  "growth_rate_viscous_estimate": ', &
         growth_rate_viscous(kr0/r0, r0, sigma, rho_l, mu_l)
      write(u,'(a)') '}'
      close(u)
      write(*,'(a,a)') 'wrote ', trim(path)
   end subroutine save_json

   subroutine do_linearity()
      type(amplitude_run_t), allocatable :: runs(:)
      type(linearity_report_t) :: rep
      real(dp), allocatable :: X(:,:), times(:), norms(:)
      real(dp) :: dt, slope, intercept, r2, tlo
      integer :: i, e, n

      if (n_cases < 2) then
         write(*,'(a)') 'linearity needs at least two --case EPS=DIR options'
         stop 1
      end if
      allocate(runs(n_cases))
      do i = 1, n_cases
         write(*,'(a,es10.3,a,a)') '[eps=', case_eps(i), '] reading ', &
            trim(case_dirs(i))
         call build_snapshot_matrix(trim(case_dirs(i)), trim(base_dir), &
              fields(1:nfields), skip, stride, maxsnap, X, times, dt, e, &
              verbose=.false., prefix=trim(plt_prefix))
         if (e /= 0) then
            write(*,'(a,i0)') 'ERROR building subspace, ierr=', e
            call explain_nsmfp_error(e)
            stop 2
         end if
         call perturbation_norms(X, norms)
         n = size(times)
         tlo = t_min
         if (tlo <= -huge(1.0_dp)/2) tlo = times(1)
         call fit_exponential_growth(times, norms, tlo, t_max, slope, &
                                     intercept, r2, e)
         runs(i)%eps = case_eps(i)
         runs(i)%growth_rate = slope
         runs(i)%r_squared = r2
         allocate(runs(i)%times(n), runs(i)%norms(n))
         runs(i)%times = times
         runs(i)%norms = norms
         deallocate(X, times, norms)
      end do

      call assess_linearity(runs, rep, e)
      if (e /= 0) then
         write(*,'(a,i0)') 'ERROR assessing linearity, ierr=', e
         stop 2
      end if
      write(*,'(a)') ''
      call print_report(runs, rep)
      write(*,'(a)') ''
      write(*,'(a,f10.6)') 'exact inviscid growth rate : ', &
         growth_rate(kr0_list(1)/r0, r0, sigma, rho_l)
      write(*,'(a,f10.6)') 'viscous estimate           : ', &
         growth_rate_viscous(kr0_list(1)/r0, r0, sigma, rho_l, mu_l)
      write(*,'(a,es12.3)') 'recommended amplitude      : ', &
         recommend_amplitude(runs, rep)

      ! Non-zero exit so a driving script notices the sweep is unusable.
      if (rep%verdict /= VERDICT_LINEAR) stop 3
   end subroutine do_linearity

   subroutine do_dispersion()
      real(dp) :: rd, rf, r2, inv, vis, best_k, best_s, kx, sx
      type(dmd_result_t) :: res
      logical :: ok
      integer :: i, nrow
      character(len=512) :: pd, bd
      character(len=64) :: tk, te
      real(dp) :: rows_k(64), rows_m(64), rows_f(64)

      nrow = 0
      do i = 1, n_kr0
         write(tk,'(a)') trim(nice(kr0_list(i)))
         write(te,'(a)') trim(nice(eps_list(1)))
         pd = trim(root)//'/k'//trim(tk)//'_eps'//trim(te)
         bd = trim(root)//'/k'//trim(tk)//'_base'
         write(*,'(a,f6.3,a)') '[k r0=', kr0_list(i), '] reading ...'
         call analyse_pair(trim(pd), trim(bd), rd, rf, r2, res, ok)
         if (.not. ok) then
            write(*,'(a,f6.3,a)') '[k r0=', kr0_list(i), '] FAILED, skipping'
            cycle
         end if
         nrow = nrow + 1
         rows_k(nrow) = kr0_list(i)
         rows_m(nrow) = rd
         rows_f(nrow) = rf
         call dmd_free(res)
      end do

      if (nrow == 0) then
         write(*,'(a)') 'no runs could be analysed'
         stop 2
      end if

      write(*,'(a)') ''
      write(*,'(a8,a14,a14,a14,a14,a12)') 'k r0', 'DMD', 'norm fit', &
         'inviscid', 'viscous', 'err vs vis'
      best_s = -huge(1.0_dp); best_k = 0
      do i = 1, nrow
         inv = growth_rate(rows_k(i)/r0, r0, sigma, rho_l)
         vis = growth_rate_viscous(rows_k(i)/r0, r0, sigma, rho_l, mu_l)
         if (vis > 0.0_dp) then
            write(*,'(f8.3,f14.6,f14.6,f14.6,f14.6,f11.2,a)') rows_k(i), &
               rows_m(i), rows_f(i), inv, vis, &
               100.0_dp*(rows_m(i)-vis)/vis, '%'
         else
            write(*,'(f8.3,f14.6,f14.6,f14.6,f14.6,a12)') rows_k(i), &
               rows_m(i), rows_f(i), inv, vis, '   n/a'
         end if
         if (rows_m(i) > best_s) then
            best_s = rows_m(i); best_k = rows_k(i)
         end if
      end do

      call most_unstable_wavenumber(r0, sigma, rho_l, kx, sx)
      write(*,'(a)') ''
      write(*,'(a,f7.3,a,f10.6,a)') 'most unstable sampled k r0 : ', best_k, &
         '  (sigma = ', best_s, ')'
      write(*,'(a,f7.3,a,f10.6,a)') 'theoretical peak           : ', kx*r0, &
         '  (sigma = ', sx, ')'
   end subroutine do_dispersion

   subroutine do_converge()
      real(dp), allocatable :: X(:,:), times(:)
      real(dp) :: dt, g, vis, lo, hi, mean
      type(dmd_result_t) :: res
      integer :: e, i, n, sizes(6), ns, cnt
      real(dp) :: rates(6)

      call build_snapshot_matrix(trim(pert_dir), trim(base_dir), &
           fields(1:nfields), skip, stride, maxsnap, X, times, dt, e, &
           verbose=.true., prefix=trim(plt_prefix))
      if (e /= 0) then
         write(*,'(a,i0)') 'ERROR building subspace, ierr=', e
         call explain_nsmfp_error(e); stop 2
      end if
      n = size(X, 2)
      sizes = [50, 100, 200, 400, 800, n]
      vis = growth_rate_viscous(kr0_list(1)/r0, r0, sigma, rho_l, mu_l)

      write(*,'(a)') ''
      write(*,'(a10,a8,a16,a14)') 'n_snap', 'rank', 'growth', 'err vs vis'
      cnt = 0
      do i = 1, 6
         ns = sizes(i)
         if (ns < 3 .or. ns > n) cycle
         call dmd_compute(X(:, 1:ns), dt, rank, tol, res, e)
         if (e /= 0) cycle
         g = leading_growth_rate(res)
         cnt = cnt + 1
         rates(cnt) = g
         if (vis > 0.0_dp) then
            write(*,'(i10,i8,f16.6,f13.2,a)') ns, res%rank, g, &
               100.0_dp*(g-vis)/vis, '%'
         else
            write(*,'(i10,i8,f16.6,a14)') ns, res%rank, g, '  n/a'
         end if
         call dmd_free(res)
      end do

      if (cnt > 1) then
         lo = minval(rates(1:cnt)); hi = maxval(rates(1:cnt))
         mean = sum(rates(1:cnt))/real(cnt,dp)
         write(*,'(a)') ''
         write(*,'(a,es12.4,a,f7.2,a)') 'spread over configurations: ', &
            hi-lo, '  (', 100.0_dp*(hi-lo)/max(abs(mean),tiny(1.0_dp)), &
            '% of the mean)'
         write(*,'(a)') 'Converged if this spread is small compared with '// &
            'the difference from theory.'
      end if
   end subroutine do_converge

   function nice(x) result(s)
      real(dp), intent(in) :: x
      character(len=64) :: s, t
      integer :: i, n
      if (x == 0.0_dp) then
         s = '0'; return
      end if
      write(t,'(f24.10)') x
      t = adjustl(t)
      n = len_trim(t)
      do i = n, 1, -1
         if (t(i:i) /= '0') exit
         t(i:i) = ' '
      end do
      n = len_trim(t)
      if (n > 0) then
         if (t(n:n) == '.') t(n:n) = ' '
      end if
      s = trim(t)
   end function nice

   subroutine usage()
      write(*,'(a)') 'usage: run_lsa <analyse|linearity|dispersion|converge> [options]'
      write(*,'(a)') '  --perturbed DIR   perturbed run directory'
      write(*,'(a)') '  --base DIR        unperturbed twin directory'
      write(*,'(a)') '  --kr0 V[,V...]    wavenumber(s) k*r0'
      write(*,'(a)') '  --eps V           perturbation amplitude'
      write(*,'(a)') '  --case EPS=DIR    amplitude case (repeatable)'
      write(*,'(a)') '  --root DIR        root for dispersion sweep'
      write(*,'(a)') '  --fields A,B      plotfile components (default velocities)'
      write(*,'(a)') '  --prefix NAME     plotfile prefix (default nddataPLT)'
      write(*,'(a)') '  --skip N --stride N --max-snapshots N'
      write(*,'(a)') '  --rank N --tol V --t-min V --t-max V'
      write(*,'(a)') '  --r0 V --sigma V --rho-l V --mu-l V'
      write(*,'(a)') '  --save PATH       write results as JSON'
   end subroutine usage

   subroutine parse_args(ierr)
      integer, intent(out) :: ierr
      integer :: na, i, p
      character(len=512) :: a, v

      ierr = 0
      na = command_argument_count()
      if (na < 1) then
         ierr = 1; return
      end if
      call get_command_argument(1, cmd)

      i = 2
      do while (i <= na)
         call get_command_argument(i, a)
         select case (trim(a))
         case ('--perturbed'); call nextv(i,na,v,ierr); pert_dir = v
         case ('--base');      call nextv(i,na,v,ierr); base_dir = v
         case ('--root');      call nextv(i,na,v,ierr); root = v
         case ('--save');      call nextv(i,na,v,ierr); out_json = v
         case ('--kr0')
            call nextv(i,na,v,ierr); call parse_rlist(v, kr0_list, n_kr0)
         case ('--eps')
            call nextv(i,na,v,ierr); call parse_rlist(v, eps_list, n_eps)
         case ('--case')
            call nextv(i,na,v,ierr)
            p = index(v, '=')
            if (p == 0) then
               write(*,'(a)') '--case expects EPS=DIR'; ierr = 2; return
            end if
            n_cases = n_cases + 1
            read(v(1:p-1), *) case_eps(n_cases)
            case_dirs(n_cases) = v(p+1:)
         case ('--prefix')
            call nextv(i,na,v,ierr); plt_prefix = v
         case ('--fields')
            call nextv(i,na,v,ierr); call parse_slist(v, fields, nfields)
         case ('--skip');   call nextv(i,na,v,ierr); read(v,*) skip
         case ('--stride'); call nextv(i,na,v,ierr); read(v,*) stride
         case ('--max-snapshots'); call nextv(i,na,v,ierr); read(v,*) maxsnap
         case ('--rank');   call nextv(i,na,v,ierr); read(v,*) rank
         case ('--tol');    call nextv(i,na,v,ierr); read(v,*) tol
         case ('--t-min');  call nextv(i,na,v,ierr); read(v,*) t_min
         case ('--t-max');  call nextv(i,na,v,ierr); read(v,*) t_max
         case ('--r0');     call nextv(i,na,v,ierr); read(v,*) r0
         case ('--sigma');  call nextv(i,na,v,ierr); read(v,*) sigma
         case ('--rho-l');  call nextv(i,na,v,ierr); read(v,*) rho_l
         case ('--mu-l');   call nextv(i,na,v,ierr); read(v,*) mu_l
         case default
            write(*,'(a,a)') 'unknown option: ', trim(a)
            ierr = 3; return
         end select
         if (ierr /= 0) return
         i = i + 1
      end do
   end subroutine parse_args

   subroutine nextv(i, na, v, ierr)
      integer, intent(inout) :: i
      integer, intent(in) :: na
      character(len=*), intent(out) :: v
      integer, intent(out) :: ierr
      ierr = 0
      i = i + 1
      if (i > na) then
         ierr = 4; return
      end if
      call get_command_argument(i, v)
   end subroutine nextv

   subroutine parse_rlist(s, v, n)
      character(len=*), intent(in) :: s
      real(dp), intent(out) :: v(:)
      integer, intent(out) :: n
      integer :: p
      character(len=512) :: t
      t = s; n = 0
      do
         p = index(t, ',')
         n = n + 1
         if (p == 0) then
            read(t,*) v(n); exit
         end if
         read(t(1:p-1),*) v(n)
         t = t(p+1:)
      end do
   end subroutine parse_rlist

   subroutine parse_slist(s, v, n)
      character(len=*), intent(in) :: s
      character(len=*), intent(out) :: v(:)
      integer, intent(out) :: n
      integer :: p
      character(len=512) :: t
      t = s; n = 0
      do
         p = index(t, ',')
         n = n + 1
         if (p == 0) then
            v(n) = trim(t); exit
         end if
         v(n) = t(1:p-1)
         t = t(p+1:)
      end do
   end subroutine parse_slist

end program run_lsa
