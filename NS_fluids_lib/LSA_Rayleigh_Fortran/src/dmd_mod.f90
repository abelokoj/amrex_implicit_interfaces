!> SVD-based Dynamic Mode Decomposition (Schmid 2010).
!>
!> Fortran counterpart of python/dmd.py.  Implements the variant Ranjan, Unnikrishnan & Gaitonde (JCP 403, 2020) use in section 2.2 to extract stability modes from an NS-MFP snapshot subspace, with the amplitude ranking of Jovanovic, Schmid & Nichols (2014).
!>
!> Algorithm, for snapshots X (n_space x n_snap) sampled every dt: A = X(:,1:m-1),  B = X(:,2:m) A = U S V^H                                   (LAPACK dgesvd) S_tilde = U_r^H B V_r S_r^{-1}                (r = truncation rank) S_tilde w = mu w                              (LAPACK zgeev) Phi = B V_r S_r^{-1} w                        (exact DMD modes) omega = log(mu)/dt Re(omega) is the growth rate, Im(omega) the circular frequency.
!>
!> For Rayleigh-Plateau the mode of interest is STATIONARY (Im(omega) = 0): the capillary instability grows monotonically rather than oscillating, hence leading_growth_rate's stationary_only default.
!> ORIENTATION
!>
!> Dynamic Mode Decomposition: the step that turns a sequence of flow fields into growth rates.
!>
!> The idea in one sentence: assume each snapshot follows from the previous one by multiplication by a fixed matrix A, then obtain the eigenvalues of A without ever forming it. If an eigenvalue lambda satisfies |lambda| > 1 the corresponding mode grows from snapshot to snapshot, and sigma = log(lambda)/dt converts that per-snapshot factor into a rate per unit time.
!>
!> The obstacle is size. Each snapshot holds thousands of values, so A would hold millions, and forming it is neither possible nor necessary. The projection below reduces the problem to the few directions the snapshots actually span, found by a singular value decomposition, and takes the eigenvalues of the resulting small matrix. Those are the eigenvalues sought; the large matrix never exists. This is the sense in which the method is Jacobian-free.
module dmd_mod
   use table_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private

   public :: dp
   public :: dmd_result_t, dmd_compute, dmd_free
   public :: leading_growth_rate, dmd_summary

   !> Outcome of a DMD, sorted by descending amplitude.
   type :: dmd_result_t
      integer :: n_modes = 0
      integer :: rank = 0
      real(dp) :: dt = 0.0_dp
      complex(dp), allocatable :: ritz(:)       !< discrete eigenvalues mu
      complex(dp), allocatable :: omega(:)      !< log(mu)/dt
      complex(dp), allocatable :: modes(:,:)    !< (n_space, n_modes)
      real(dp),    allocatable :: amplitude(:)  !< |b|
      real(dp),    allocatable :: svals(:)      !< singular values of A
   end type dmd_result_t

contains

   subroutine dmd_free(res)
      type(dmd_result_t), intent(inout) :: res
      if (allocated(res%ritz))      deallocate(res%ritz)
      if (allocated(res%omega))     deallocate(res%omega)
      if (allocated(res%modes))     deallocate(res%modes)
      if (allocated(res%amplitude)) deallocate(res%amplitude)
      if (allocated(res%svals))     deallocate(res%svals)
      res%n_modes = 0
      res%rank = 0
   end subroutine dmd_free

   !> Compute the DMD of a snapshot matrix.
   !>
   !> @param X      (n_space, n_snap) snapshots; column j is time (j-1)*dt @param dt     sampling interval @param rank   truncation rank; <= 0 selects by tolerance @param tol    relative singular-value cutoff when rank <= 0 @param res    result @param ierr   0 on success, non-zero on failure
   subroutine dmd_compute(X, dt, rank, tol, res, ierr)
      real(dp), intent(in) :: X(:,:)
      real(dp), intent(in) :: dt, tol
      integer,  intent(in) :: rank
      type(dmd_result_t), intent(out) :: res
      integer, intent(out) :: ierr

      integer :: n, m, mm, r, i, j, info, lwork
      real(dp), allocatable :: A(:,:), B(:,:), U(:,:), VT(:,:), S(:)
      real(dp), allocatable :: rwork(:), work(:)
      complex(dp), allocatable :: St(:,:), w(:,:), vl(:,:), cwork(:)
      complex(dp), allocatable :: BVS(:,:), rhs(:), lhs(:,:)
      complex(dp), allocatable :: mu(:)
      real(dp), allocatable :: BVSr(:,:)
      integer, allocatable :: order(:)

      ierr = 0
      n = size(X, 1)
      m = size(X, 2)

      if (m < 2) then
         ierr = 1; return            ! need at least two snapshots
      end if
      if (dt <= 0.0_dp) then
         ierr = 2; return
      end if
      if (n < 1) then
         ierr = 3; return
      end if
      if (any(X /= X)) then          ! NaN check
         ierr = 4; return
      end if

      mm = m - 1
      allocate(A(n, mm), B(n, mm))
      A = X(:, 1:mm)
      B = X(:, 2:m)

      ! ---- SVD of A:  A = U S V^T -------------------------------------
      allocate(U(n, min(n, mm)), VT(min(n, mm), mm), S(min(n, mm)))
      lwork = max(1, 10 * max(n, mm))
      allocate(work(lwork))
      call dgesvd('S', 'S', n, mm, A, n, S, U, n, VT, min(n, mm), &
                  work, lwork, info)
      if (info /= 0) then
         ierr = 10 + info; return
      end if
      deallocate(work)

      ! ---- truncation rank --------------------------------------------
      if (rank > 0) then
         r = min(rank, size(S))
      else
         r = 0
         do i = 1, size(S)
            if (S(i) > tol * S(1)) r = r + 1
         end do
         r = max(r, 1)
      end if

      ! ---- BVS = B V_r S_r^{-1}  (real) -------------------------------- V = VT^T, so V(:,i) = VT(i,:).
      allocate(BVSr(n, r))
      BVSr = 0.0_dp
      do i = 1, r
         do j = 1, mm
            BVSr(:, i) = BVSr(:, i) + B(:, j) * VT(i, j)
         end do
         BVSr(:, i) = BVSr(:, i) / S(i)
      end do

      ! ---- S_tilde = U_r^T * BVS  (r x r) ------------------------------
      allocate(St(r, r))
      do i = 1, r
         do j = 1, r
            St(i, j) = cmplx(dot_product(U(:, i), BVSr(:, j)), 0.0_dp, dp)
         end do
      end do

      ! ---- eigen-decomposition of S_tilde ------------------------------
      allocate(mu(r), w(r, r), vl(1, 1))
      lwork = max(1, 4 * r)
      allocate(cwork(lwork), rwork(max(1, 2 * r)))
      call zgeev('N', 'V', r, St, r, mu, vl, 1, w, r, cwork, lwork, rwork, info)
      if (info /= 0) then
         ierr = 20 + info; return
      end if
      deallocate(cwork, rwork)

      ! ---- exact DMD modes:  Phi = BVS * w ------------------------------
      allocate(BVS(n, r))
      BVS = cmplx(BVSr, 0.0_dp, dp)
      allocate(res%modes(n, r))
      res%modes = matmul(BVS, w)

      ! ---- amplitudes: least-squares fit of mode 1 to snapshot 1 -------
      allocate(lhs(n, r), rhs(max(n, r)))
      lhs = res%modes
      rhs = cmplx(0.0_dp, 0.0_dp, dp)
      rhs(1:n) = cmplx(X(:, 1), 0.0_dp, dp)
      allocate(res%amplitude(r))
      call solve_lstsq(lhs, rhs, n, r, res%amplitude, info)
      if (info /= 0) then
         ierr = 30 + info; return
      end if

      ! ---- continuous-time eigenvalues ---------------------------------
      allocate(res%ritz(r), res%omega(r))
      res%ritz = mu
      do i = 1, r
         if (abs(mu(i)) > 0.0_dp) then
            res%omega(i) = log(mu(i)) / dt
         else
            ! a mode annihilated in one step: -inf growth
            res%omega(i) = cmplx(-huge(1.0_dp), 0.0_dp, dp)
         end if
      end do

      res%n_modes = r
      res%rank = r
      res%dt = dt
      allocate(res%svals(size(S)))
      res%svals = S

      ! ---- sort by descending amplitude --------------------------------
      allocate(order(r))
      call argsort_desc(res%amplitude, order)
      res%ritz      = res%ritz(order)
      res%omega     = res%omega(order)
      res%amplitude = res%amplitude(order)
      res%modes     = res%modes(:, order)
   end subroutine dmd_compute

   !> Least-squares solve of lhs * b = rhs, returning |b|.
   subroutine solve_lstsq(lhs, rhs, n, r, absb, info)
      complex(dp), intent(inout) :: lhs(:,:), rhs(:)
      integer, intent(in) :: n, r
      real(dp), intent(out) :: absb(:)
      integer, intent(out) :: info

      complex(dp), allocatable :: cwork(:)
      real(dp), allocatable :: rwork(:), sv(:)
      integer :: lwork, rnk, i

      allocate(sv(min(n, r)), rwork(max(1, 5 * min(n, r))))
      lwork = max(1, 2 * min(n, r) + max(n, r) + 64 * max(n, r))
      allocate(cwork(lwork))
      call zgelss(n, r, 1, lhs, size(lhs, 1), rhs, size(rhs), sv, &
                  -1.0_dp, rnk, cwork, lwork, rwork, info)
      if (info == 0) then
         do i = 1, r
            absb(i) = abs(rhs(i))
         end do
      end if
   end subroutine solve_lstsq

   !> Indices that sort v in descending order (insertion sort; r is small).
   pure subroutine argsort_desc(v, idx)
      real(dp), intent(in) :: v(:)
      integer, intent(out) :: idx(:)
      integer :: i, j, key
      do i = 1, size(v)
         idx(i) = i
      end do
      do i = 2, size(v)
         key = idx(i)
         j = i - 1
         do while (j >= 1)
            if (v(idx(j)) >= v(key)) exit
            idx(j + 1) = idx(j)
            j = j - 1
         end do
         idx(j + 1) = key
      end do
   end subroutine argsort_desc

   !> Largest growth rate among physically retained modes.
   !>
   !> For Rayleigh-Plateau use stationary_only = .true. (the default sense): the capillary mode is non-oscillatory, and restricting to it avoids picking up weak spurious oscillatory pairs.  Returns -huge if nothing survives the filter.
   function leading_growth_rate(res, stationary_only, freq_tol, amp_frac) &
         result(s)
      type(dmd_result_t), intent(in) :: res
      logical,  intent(in), optional :: stationary_only
      real(dp), intent(in), optional :: freq_tol, amp_frac
      real(dp) :: s
      logical  :: stat_only
      real(dp) :: ftol, afrac, amax, g, f
      integer  :: i

      stat_only = .true.
      if (present(stationary_only)) stat_only = stationary_only
      ftol = 1.0e-6_dp
      if (present(freq_tol)) ftol = freq_tol
      afrac = 1.0e-3_dp
      if (present(amp_frac)) afrac = amp_frac

      s = -huge(1.0_dp)
      if (res%n_modes < 1) return
      amax = maxval(res%amplitude)

      do i = 1, res%n_modes
         if (res%amplitude(i) < afrac * amax) cycle
         f = aimag(res%omega(i))
         if (stat_only .and. abs(f) > ftol) cycle
         g = real(res%omega(i), dp)
         if (g > s) s = g
      end do
   end function leading_growth_rate

   !> Print a table of the leading modes.
   subroutine dmd_summary(res, nshow, unit)
      type(dmd_result_t), intent(in) :: res
      integer, intent(in), optional  :: nshow, unit
      integer :: n, u, i
      type(table_t) :: t

      n = min(10, res%n_modes)
      if (present(nshow)) n = min(nshow, res%n_modes)
      u = 6
      if (present(unit)) u = unit

      if (u == 6) then
         call t%init([cel('mode '), cel('growth sigma'), cel('freq omega '), &
            cel('|lambda| '), cel('amplitude ')], 'rrrrr')
         do i = 1, n
            call t%row([cel(fmt_i(i-1)), cel(fmt_f(real(res%omega(i), dp), 6)), &
               cel(fmt_f(aimag(res%omega(i))/(2.0_dp*3.141592653589793_dp), 6)), &
               cel(fmt_f(abs(res%ritz(i)), 6)), cel(fmt_e(res%amplitude(i)))])
         end do
         call t%emit('Decomposition: rank '//fmt_i(res%rank)//', dt '// &
                     fmt_e(res%dt)//', '//fmt_i(res%n_modes)//' modes')
      else
         write(u, '(a,i0,a,es12.5,a,i0)') 'rank=', res%rank, '  dt=', res%dt, &
            '  n_modes=', res%n_modes
         write(u, '(a4,a16,a16,a12,a14)') '#', 'growth sigma', 'freq omega', &
            '|lambda|', 'amplitude'
         do i = 1, n
            write(u, '(i4,es16.6,es16.6,f12.6,es14.4)') i - 1, &
               real(res%omega(i), dp), aimag(res%omega(i)) / &
               (2.0_dp * 3.141592653589793_dp), &
               abs(res%ritz(i)), res%amplitude(i)
         end do
      end if
   end subroutine dmd_summary

end module dmd_mod
