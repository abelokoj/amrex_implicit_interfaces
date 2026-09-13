!> NS-MFP snapshot assembly for the Rayleigh-Plateau instability.
!>
!> Fortran counterpart of python/nsmfp.py.  Implements the subspace-generation half of Ranjan, Unnikrishnan & Gaitonde (JCP 403, 2020) on top of the amrex_implicit_interfaces two-phase solver, WITHOUT modifying the solver.
!>
!> The B_f body force, and why the solver requires no modification
!! -----------------------------------------------------
!! Ranjan et al. add a restraining body force
!!
!! B_f = dQbar/dt - F(Qbar)                              (their Eq. 3)
!!
!! whose sole role, visible in their Eq. (6), is to annihilate the bracketed base-state drift term.  We achieve the same cancellation by running an unperturbed TWIN alongside the perturbed case from an identical restart, with identical discretisation, and subtracting snapshot by snapshot:
!!
!! Q'(t_n) = Q_perturbed(t_n) - Q_unperturbed(t_n).
!!
!! Both runs suffer the same drift, so the difference removes it to machine precision and Q' obeys the linearised dynamics to O(eps^2), exactly as in their Eq. (7)-(8).  This is the "subtract the base flow from the total flow" step of their Fig. 1; the drift that B_f cancels by stored forcing is here cancelled by a computed twin.
!!
!! For this problem the substitution is more than a convenience.  A quiescent liquid column is an EXACT equilibrium, so the only base-state drift is discretisation error - principally the parasitic (spurious) currents of the discrete curvature/surface-tension operator.  Those are precisely what the twin measures and removes, and unlike a frozen B_f the twin tracks their evolution.
!!
!! Validity conditions, enforced by plotfile_mod::compatible:
!! - identical grid, domain, field set and times (ns.fixed_dt, max_level=0)
!! - same MPI rank count for both runs
module nsmfp_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use plotfile_mod
   implicit none
   private

   public :: dp
   public :: build_snapshot_matrix, perturbation_norms
   public :: fit_exponential_growth, uniform_dt

contains

   !> Assemble the NS-MFP subspace Q'(t_n) = Q_pert(t_n) - Q_unpert(t_n).
   !>
   !> Snapshots are paired strictly by STEP NUMBER, so a plotfile missing on one side cannot silently shift the pairing.
   !>
   !> @param pert_dir, base_dir  directories holding the twin runs' plt* output @param fields   component names forming the state vector @param skip     discard this many leading snapshots (shed the transient) @param stride   use every stride-th snapshot @param maxsnap  cap on retained snapshots (<= 0 for no cap) @param X        (n_space, n_snap) subspace, allocated here @param times    snapshot times @param dt       uniform sampling interval @param ierr     0 on success
   subroutine build_snapshot_matrix(pert_dir, base_dir, fields, skip, stride, &
                                    maxsnap, X, times, dt, ierr, verbose, &
                                    prefix)
      character(len=*), intent(in) :: pert_dir, base_dir
      character(len=*), intent(in) :: fields(:)
      integer, intent(in) :: skip, stride, maxsnap
      real(dp), allocatable, intent(out) :: X(:,:), times(:)
      real(dp), intent(out) :: dt
      integer, intent(out) :: ierr
      logical, intent(in), optional :: verbose
      !> Plotfile directory prefix.  The solver does NOT honour amr.plot_file: with ns.visual_nddata_format=1 it writes nddataPLT00000023 (and MOF_PLT*), not plt00023.  Confirmed by running the solver.  Default therefore is 'nddataPLT'.
      character(len=*), intent(in), optional :: prefix

      character(len=MAX_NAME), allocatable :: ppaths(:), bpaths(:)
      integer :: np, nb, i, j, k, n_pair, n_use, nsp, idx
      integer, allocatable :: pidx(:), bidx(:)
      type(plotfile_t) :: pa, pb
      character(len=256) :: msg
      logical :: verb
      real(dp) :: nrm
      character(len=32) :: pfx

      verb = .false.
      if (present(verbose)) verb = verbose
      pfx = 'nddataPLT'
      if (present(prefix)) pfx = prefix
      ierr = 0
      dt = 0.0_dp

      call find_plotfiles(pert_dir, ppaths, np, trim(pfx))
      call find_plotfiles(base_dir, bpaths, nb, trim(pfx))
      if (np == 0) then
         ierr = 1; return                    ! no nddataPLT* under perturbed dir
      end if
      if (nb == 0) then
         ierr = 2; return                    ! no nddataPLT* under base dir
      end if

      ! ---- pair by step number ------------------------------------------
      allocate(pidx(np), bidx(np))
      n_pair = 0
      do i = 1, np
         do j = 1, nb
            if (plotfile_step(ppaths(i)) == plotfile_step(bpaths(j))) then
               n_pair = n_pair + 1
               pidx(n_pair) = i
               bidx(n_pair) = j
               exit
            end if
         end do
      end do
      if (n_pair < 2) then
         ierr = 3; return                    ! need >= 2 matching steps
      end if

      ! ---- apply skip / stride / cap ------------------------------------
      n_use = 0
      do i = skip + 1, n_pair, max(stride, 1)
         n_use = n_use + 1
         if (maxsnap > 0 .and. n_use >= maxsnap) exit
      end do
      if (n_use < 2) then
         ierr = 4; return                    ! too few after skip/stride
      end if

      ! ---- read and subtract ---------------------------------------------
      allocate(times(n_use))
      k = 0
      do i = skip + 1, n_pair, max(stride, 1)
         k = k + 1
         if (k > n_use) exit

         call read_plotfile(trim(ppaths(pidx(i))), pa, ierr, fields)
         if (ierr /= 0) then
            ierr = 1000 + ierr; return
         end if
         call read_plotfile(trim(bpaths(bidx(i))), pb, ierr, fields)
         if (ierr /= 0) then
            ierr = 2000 + ierr; return
         end if
         if (.not. compatible(pa, pb, 1.0e-9_dp, msg)) then
            if (verb) write(*,'(a)') 'incompatible plotfiles: '//trim(msg)
            ierr = 5; return
         end if

         nsp = pa%nx * pa%ny * pa%ncomp
         if (.not. allocated(X)) allocate(X(nsp, n_use))

         idx = 0
         do j = 1, pa%ncomp
            X(idx+1:idx+pa%nx*pa%ny, k) = &
               reshape(pa%data(:,:,j) - pb%data(:,:,j), [pa%nx*pa%ny])
            idx = idx + pa%nx * pa%ny
         end do
         times(k) = pa%time

         if (verb .and. (mod(k-1, 25) == 0 .or. k == n_use)) then
            nrm = sqrt(sum(X(:,k)**2))
            write(*,'(a,i0,a,i0,a,es12.5,a,es12.5)') '  [', k, '/', n_use, &
               '] t=', times(k), '  |Q''|=', nrm
         end if

         call free_plotfile(pa)
         call free_plotfile(pb)
      end do

      call uniform_dt(times, dt, ierr)
      if (ierr /= 0) return

      if (verb) then
         write(*,'(a,i0,a,i0,a,es12.5)') 'subspace: ', size(X,1), &
            ' dof x ', size(X,2), ' snapshots, dt=', dt
      end if
   end subroutine build_snapshot_matrix

   !> Sampling interval, insisting the snapshots are evenly spaced.
   !>
   !> DMD assumes a constant dt; non-uniform sampling biases every eigenvalue, so this is checked rather than silently averaged.
   subroutine uniform_dt(times, dt, ierr)
      real(dp), intent(in) :: times(:)
      real(dp), intent(out) :: dt
      integer, intent(out) :: ierr
      integer :: i, n
      real(dp) :: d, spread

      ierr = 0
      n = size(times)
      if (n < 2) then
         ierr = 10; dt = 0.0_dp; return
      end if
      do i = 2, n
         if (times(i) <= times(i-1)) then
            ierr = 11; return                ! not strictly increasing
         end if
      end do
      dt = (times(n) - times(1)) / real(n - 1, dp)
      spread = 0.0_dp
      do i = 2, n
         d = times(i) - times(i-1)
         spread = max(spread, abs(d - dt))
      end do
      if (spread > 1.0e-6_dp * max(dt, tiny(1.0_dp))) then
         ierr = 12                           ! non-uniform sampling
      end if
   end subroutine uniform_dt

   !> Column-wise L2 norms of the subspace.
   !>
   !> On a linear, exponentially growing mode this is a straight line on a log plot with slope equal to the growth rate, giving an independent check on the DMD result.
   subroutine perturbation_norms(X, norms)
      real(dp), intent(in) :: X(:,:)
      real(dp), allocatable, intent(out) :: norms(:)
      integer :: j
      allocate(norms(size(X,2)))
      do j = 1, size(X,2)
         norms(j) = sqrt(sum(X(:,j)**2))
      end do
   end subroutine perturbation_norms

   !> Least-squares slope of log|Q'| versus t over an optional window.
   !>
   !> A poor r^2 usually means the window still contains the initial transient, or that nonlinearity has set in at the late end.
   subroutine fit_exponential_growth(times, norms, t_min, t_max, &
                                     slope, intercept, r2, ierr)
      real(dp), intent(in) :: times(:), norms(:)
      real(dp), intent(in) :: t_min, t_max
      real(dp), intent(out) :: slope, intercept, r2
      integer, intent(out) :: ierr

      integer :: i, n, m
      real(dp) :: st, sy, stt, sty, t, y, ybar, ssres, sstot, pred

      ierr = 0
      n = size(times)
      m = 0
      st = 0; sy = 0; stt = 0; sty = 0
      do i = 1, n
         if (norms(i) <= 0.0_dp) cycle
         if (times(i) < t_min) cycle
         if (t_max > t_min .and. times(i) > t_max) cycle
         t = times(i); y = log(norms(i))
         m = m + 1
         st = st + t; sy = sy + y
         stt = stt + t*t; sty = sty + t*y
      end do
      if (m < 2) then
         ierr = 20; slope = 0; intercept = 0; r2 = 0; return
      end if

      slope = (real(m,dp)*sty - st*sy) / (real(m,dp)*stt - st*st)
      intercept = (sy - slope*st) / real(m,dp)

      ybar = sy / real(m,dp)
      ssres = 0; sstot = 0
      do i = 1, n
         if (norms(i) <= 0.0_dp) cycle
         if (times(i) < t_min) cycle
         if (t_max > t_min .and. times(i) > t_max) cycle
         y = log(norms(i))
         pred = slope*times(i) + intercept
         ssres = ssres + (y - pred)**2
         sstot = sstot + (y - ybar)**2
      end do
      if (sstot > 0.0_dp) then
         r2 = 1.0_dp - ssres/sstot
      else
         r2 = 0.0_dp
      end if
   end subroutine fit_exponential_growth

end module nsmfp_mod
