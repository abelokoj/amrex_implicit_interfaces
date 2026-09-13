!> End-to-end test: real AMReX plotfiles -> twin subtraction -> DMD -> sigma.
!>
!> Plotfiles are written in the genuine on-disk AMReX format and read back through the production path (plotfile_mod), so a format or ordering mistake in the hand-written parser shows up as a wrong growth rate rather than passing silently.  The same synthetic files are readable by yt, which is how the format was cross-checked independently of this code.
!>
!> The planted signal mimics the real experiment:
!>
!> Q_unperturbed(t) = drift(t)                       (parasitic currents) Q_perturbed(t)   = drift(t) + eps e^{sigma t} phi(z)
!>
!> so the twin subtraction must remove `drift` exactly and leave a clean exponentially growing capillary mode.  `drift` is made 100x larger than the perturbation - as spurious currents often are in a surface-tension calculation - so the test fails if the subtraction is skipped or mispaired.
program test_end_to_end
   use testing_mod
   use plotfile_mod
   use plotfile_writer_mod
   use nsmfp_mod
   use dmd_mod
   use dispersion_mod, only: growth_rate, PI
   implicit none

   integer, parameter :: NX = 10, NY = 32, NSNAP = 60
   real(dp), parameter :: KR0 = 0.7_dp, R0 = 1.0_dp
   real(dp), parameter :: EPSA = 1.0e-3_dp, DT = 0.05_dp

   real(dp) :: lambda, sexact, lo(2), hi(2)
   character(len=64) :: names(2), want(2)
   character(len=256) :: pdir, bdir
   real(dp), allocatable :: X(:,:), times(:), norms(:)
   real(dp) :: dt_out, slope, intercept, r2
   type(dmd_result_t) :: res
   type(plotfile_t) :: pf
   integer :: ierr, i, j, k, step
   real(dp) :: shape_x(NX,NY), shape_y(NX,NY), d(NX,NY,2)
   real(dp) :: t, growth, expected0, drift_norm
   character(len=MAX_NAME), allocatable :: paths(:)
   integer :: nfound

   lambda = 2.0_dp * PI / (KR0 / R0)
   sexact = growth_rate(KR0/R0, R0, 1.0_dp, 1.0_dp)
   names(1) = 'x_velocity'; names(2) = 'y_velocity'
   want(1) = 'x_velocity';  want(2) = 'y_velocity'
   lo = [0.0_dp, 0.0_dp]
   hi = [4.0_dp*R0, lambda]
   pdir = '/tmp/lsa_e2e/pert'
   bdir = '/tmp/lsa_e2e/base'

   call execute_command_line('rm -rf /tmp/lsa_e2e')
   call execute_command_line('mkdir -p '//trim(pdir)//' '//trim(bdir))

   call build_shapes()
   call write_runs()

   call section('plotfile discovery')
   call find_plotfiles(trim(pdir), paths, nfound)
   call check_int(nfound, NSNAP, 'finds every perturbed plotfile')
   call check(plotfile_step(paths(1)) < plotfile_step(paths(nfound)), &
              'plotfiles are sorted by step number')
   call check_int(plotfile_step('/x/plt00120'), 120, 'parses a step number')
   call check_int(plotfile_step('/x/plt00009'), 9, 'parses a padded step')

   call section('single plotfile read')
   call read_plotfile(trim(paths(4)), pf, ierr, want)
   call check_int(ierr, 0, 'read succeeds')
   call check_int(pf%nx, NX, 'nx correct')
   call check_int(pf%ny, NY, 'ny correct')
   call check_int(pf%ncomp, 2, 'component count correct')
   call check_close(pf%time, 3.0_dp*DT, 1.0e-12_dp, 'time correct')
   call check_int(pf%step, 30, 'step correct')
   call free_plotfile(pf)

   call section('twin subtraction removes the drift')
   call build_snapshot_matrix(trim(pdir), trim(bdir), want, 0, 1, 0, &
                              X, times, dt_out, ierr, verbose=.false.)
   call check_int(ierr, 0, 'subspace assembly succeeds')
   ! Guard: if assembly failed, X is unallocated and every later check would dereference it.  Stop cleanly rather than segfaulting on a bad diagnosis.
   if (ierr /= 0) then
      write(*,'(a)') '  cannot continue without a subspace'
      call summary('test_end_to_end')
   end if
   call check_close(dt_out, DT, 1.0e-12_dp, 'uniform dt detected')
   call check_int(size(X,2), NSNAP, 'all snapshots retained')
   call check_int(size(X,1), NX*NY*2, 'state vector size = nx*ny*ncomp')

   expected0 = EPSA * sqrt(sum(shape_x**2) + sum(shape_y**2))
   call check_close_rel(norm2(X(:,1)), expected0, 1.0e-10_dp, &
                        'first snapshot equals the planted perturbation')
   drift_norm = sqrt(sum(drift(0.0_dp)**2))
   call check(drift_norm > 100.0_dp*norm2(X(:,1)), &
              'drift really was 100x larger than the perturbation')

   call section('growth rate from the norm fit')
   call perturbation_norms(X, norms)
   call fit_exponential_growth(times, norms, times(1), -1.0_dp, &
                               slope, intercept, r2, ierr)
   call check_int(ierr, 0, 'fit succeeds')
   call check_close_rel(slope, sexact, 1.0e-8_dp, &
                        'norm slope equals the exact Rayleigh rate')
   call check(r2 > 0.999999_dp, 'fit is essentially perfect (r^2 > 0.999999)')

   call section('growth rate from DMD (headline check)')
   call dmd_compute(X, dt_out, 0, 1.0e-10_dp, res, ierr)
   call check_int(ierr, 0, 'dmd succeeds')
   call check_close_rel(leading_growth_rate(res), sexact, 1.0e-6_dp, &
                        'DMD recovers the exact Rayleigh growth rate')
   call check(leading_growth_rate(res) > 0.3_dp .and. &
              leading_growth_rate(res) < 0.4_dp, &
              'growth rate has the expected magnitude ~0.343')
   call dmd_free(res)
   deallocate(X, times, norms)

   call section('single-field subspace')
   ! Ranjan et al. note one well-chosen variable can suffice, cutting memory up to 80%.  The pipeline must honour a single-component request.
   call build_snapshot_matrix(trim(pdir), trim(bdir), want(1:1), 0, 1, 0, &
                              X, times, dt_out, ierr, verbose=.false.)
   call check_int(ierr, 0, 'single-field assembly succeeds')
   call check_int(size(X,1), NX*NY, 'state vector holds one component')
   call dmd_compute(X, dt_out, 0, 1.0e-10_dp, res, ierr)
   call check_close_rel(leading_growth_rate(res), sexact, 1.0e-6_dp, &
                        'single field still gives the right growth rate')
   call dmd_free(res); deallocate(X, times)

   call section('skip and stride')
   call build_snapshot_matrix(trim(pdir), trim(bdir), want, 10, 3, 0, &
                              X, times, dt_out, ierr, verbose=.false.)
   call check_int(ierr, 0, 'skip/stride assembly succeeds')
   call check_close(dt_out, 3.0_dp*DT, 1.0e-12_dp, 'stride scales dt')
   call dmd_compute(X, dt_out, 0, 1.0e-10_dp, res, ierr)
   call check_close_rel(leading_growth_rate(res), sexact, 1.0e-6_dp, &
                        'strided subspace gives the same growth rate')
   call dmd_free(res); deallocate(X, times)

   call section('guard rails')
   ! identical runs must give exactly zero: the check the driver uses to detect a perturbed deck that forgot to set radblob
   call build_snapshot_matrix(trim(bdir), trim(bdir), want, 0, 1, 0, &
                              X, times, dt_out, ierr, verbose=.false.)
   call check_int(ierr, 0, 'self-subtraction assembles')
   call check_close(maxval(abs(X)), 0.0_dp, 0.0_dp, &
                    'subtracting a run from itself gives exactly zero')
   deallocate(X, times)

   call test_grid_mismatch()
   call test_time_mismatch()
   call test_missing_dir()
   call test_nonuniform()

   call summary('test_end_to_end')

contains

   subroutine build_shapes()
      integer :: i, j
      real(dp) :: z, r, taper
      do j = 1, NY
         z = (real(j,dp) - 0.5_dp) * lambda / real(NY,dp)
         do i = 1, NX
            r = (real(i,dp) - 0.5_dp) * (4.0_dp*R0) / real(NX,dp)
            taper = exp(-((r - R0)**2) / 0.25_dp)
            shape_x(i,j) = taper * cos(KR0 * z / R0)
            shape_y(i,j) = 0.6_dp * shape_x(i,j)
         end do
      end do
   end subroutine build_shapes

   !> Base-state drift (parasitic currents): large, slow, non-modal.
   function drift(t) result(a)
      real(dp), intent(in) :: t
      real(dp) :: a(NX,NY)
      integer :: i, j
      real(dp) :: z, r
      do j = 1, NY
         z = (real(j,dp) - 0.5_dp) * lambda / real(NY,dp)
         do i = 1, NX
            r = (real(i,dp) - 0.5_dp) * (4.0_dp*R0) / real(NX,dp)
            a(i,j) = 0.5_dp * sin(3.0_dp*r) * cos(2.0_dp*z) * (1.0_dp + &
               0.1_dp*t)
         end do
      end do
   end function drift

   subroutine write_runs()
      integer :: j, step
      real(dp) :: t, growth
      character(len=256) :: p
      do j = 1, NSNAP
         t = real(j-1, dp) * DT
         step = (j-1) * 10
         growth = EPSA * exp(sexact * t)

         d(:,:,1) = drift(t)
         d(:,:,2) = drift(t)
         write(p,'(a,a,i5.5)') trim(bdir), '/nddataPLT', step
         call write_plotfile(trim(p), d, names, t, step, lo, hi)

         d(:,:,1) = drift(t) + growth*shape_x
         d(:,:,2) = drift(t) + growth*shape_y
         write(p,'(a,a,i5.5)') trim(pdir), '/nddataPLT', step
         call write_plotfile(trim(p), d, names, t, step, lo, hi)
      end do
   end subroutine write_runs

   subroutine test_grid_mismatch()
      real(dp) :: a(8,16,1), b(8,32,1), l(2), h(2)
      character(len=64) :: nm(1)
      real(dp), allocatable :: XX(:,:), tt(:)
      real(dp) :: dd
      integer :: e
      nm(1) = 'x_velocity'
      l = [0.0_dp,0.0_dp]; h = [1.0_dp,1.0_dp]
      a = 1.0_dp; b = 1.0_dp
      call execute_command_line('rm -rf /tmp/lsa_gm && mkdir -p /tmp/lsa_gm/a /tmp/lsa_gm/b')
      call write_plotfile('/tmp/lsa_gm/a/nddataPLT00000', a, nm, 0.0_dp, 0, l, h)
      call write_plotfile('/tmp/lsa_gm/a/nddataPLT00010', a, nm, DT, 10, l, h)
      call write_plotfile('/tmp/lsa_gm/b/nddataPLT00000', b, nm, 0.0_dp, 0, l, h)
      call write_plotfile('/tmp/lsa_gm/b/nddataPLT00010', b, nm, DT, 10, l, h)
      call build_snapshot_matrix('/tmp/lsa_gm/a', '/tmp/lsa_gm/b', nm, &
                                 0, 1, 0, XX, tt, dd, e, verbose=.false.)
      call check_int(e, 5, 'grid mismatch is rejected')
   end subroutine test_grid_mismatch

   subroutine test_time_mismatch()
      real(dp) :: a(8,16,1), l(2), h(2)
      character(len=64) :: nm(1)
      real(dp), allocatable :: XX(:,:), tt(:)
      real(dp) :: dd
      integer :: e
      nm(1) = 'x_velocity'
      l = [0.0_dp,0.0_dp]; h = [1.0_dp,1.0_dp]
      a = 1.0_dp
      call execute_command_line('rm -rf /tmp/lsa_tm && mkdir -p /tmp/lsa_tm/a /tmp/lsa_tm/b')
      call write_plotfile('/tmp/lsa_tm/a/nddataPLT00000', a, nm, 0.0_dp, 0, l, h)
      call write_plotfile('/tmp/lsa_tm/a/nddataPLT00010', a, nm, 0.05_dp, 10, l, h)
      call write_plotfile('/tmp/lsa_tm/b/nddataPLT00000', a, nm, 0.0_dp, 0, l, h)
      call write_plotfile('/tmp/lsa_tm/b/nddataPLT00010', a, nm, 0.07_dp, 10, l, h)
      call build_snapshot_matrix('/tmp/lsa_tm/a', '/tmp/lsa_tm/b', nm, &
                                 0, 1, 0, XX, tt, dd, e, verbose=.false.)
      call check_int(e, 5, 'time mismatch is rejected')
   end subroutine test_time_mismatch

   subroutine test_missing_dir()
      character(len=64) :: nm(1)
      real(dp), allocatable :: XX(:,:), tt(:)
      real(dp) :: dd
      integer :: e
      nm(1) = 'x_velocity'
      call build_snapshot_matrix('/tmp/lsa_nope_a', '/tmp/lsa_nope_b', nm, &
                                 0, 1, 0, XX, tt, dd, e, verbose=.false.)
      call check_int(e, 1, 'missing directory is rejected')
   end subroutine test_missing_dir

   subroutine test_nonuniform()
      real(dp) :: tt(4), dd
      integer :: e
      tt = [0.0_dp, 0.1_dp, 0.2_dp, 0.35_dp]
      call uniform_dt(tt, dd, e)
      call check_int(e, 12, 'non-uniform sampling is rejected')
      tt = [0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp]
      call uniform_dt(tt, dd, e)
      call check_int(e, 0, 'uniform sampling is accepted')
      call check_close(dd, 0.1_dp, 1.0e-14_dp, 'dt computed correctly')
      tt = [0.0_dp, 0.1_dp, 0.05_dp, 0.3_dp]
      call uniform_dt(tt, dd, e)
      call check_int(e, 11, 'non-increasing times are rejected')
   end subroutine test_nonuniform

end program test_end_to_end
