!> Interface-mode amplitude: the correct observable for an interfacial instability.  Fortran counterpart of python/interface_mode.py.
!>
!> WHY THIS MODULE EXISTS (learned from running the solver)
!! --------------------------------------------------------
!! The obvious NS-MFP observable, |Q'(t)| = |Q_pert - Q_base|, fails for Rayleigh-Plateau in two different ways depending on the field chosen, and both failures only became visible with real solver output:
!!
!! 1. VELOCITY FIELDS give Q'(0) = 0 exactly.  The perturbation lives on the interface, so both twins start from rest.  That breaks the DMD amplitude ranking, which projects each mode onto the FIRST snapshot: all amplitudes come out zero and the ordering is meaningless. Velocities are still a valid state vector, but only after the startup transient is skipped.
!!
!! 2. THE RAW LEVEL SET dilutes the signal.  L0101 is a signed distance function, so Q' ~ eps*sin(k z) across the WHOLE domain while only the near-interface part participates in the instability.  In a real run the raw level-set norm stayed flat to 0.6% while the interface amplitude genuinely grew by 37%.  A flat norm does NOT prove a stable interface - it can mean the observable is wrong.
!!
!! The fix is to measure what the dispersion relation predicts: the amplitude of the k-Fourier component of the interface radius r(z,t).  For r = r0 + a(t) cos(k z), linear theory gives a(t) = a(0) exp(sigma t), so a log-linear fit of a(t) returns sigma directly, immune to both the zero initial condition and the static-offset dilution.
!!
!! Use alongside DMD, not instead of it: DMD gives the spectrum and mode shapes, this gives a robust scalar growth rate to check it against.
!> ORIENTATION
!>
!> This module measures the shape of the interface, which is the quantity the growth rate is fitted to.
!>
!> The solver stores the interface implicitly, as the zero contour of a level-set field that is negative inside the liquid and positive outside. To find the interface radius at one axial station, walk outwards along that row until the sign changes, then interpolate linearly between the two straddling cell centres. Repeating at every station gives the profile r(z), and the amplitude of its first Fourier harmonic is the observable the analysis uses.
!>
!> That observable was chosen deliberately. A velocity-based measure has no signal at the first instant, because both twin calculations start from rest and their difference is identically zero. The raw level-set difference spreads the discrepancy across the whole domain and dilutes it. The interface amplitude is the quantity the theory predicts.
module interface_mode_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use plotfile_mod
   use dispersion_mod, only: PI
   implicit none
   private

   public :: dp
   public :: interface_radius, mode_amplitude, amplitude_series
   public :: fit_growth_rate_amp, local_growth_rate, has_converged

contains

   !> Interface radius r(z) from the zero contour of the level set.
   !>
   !> Linear interpolation of the first sign change in each axial column. Columns with no sign change (interface outside the domain) are set to a negative sentinel rather than a silently wrong value.
   subroutine interface_radius(pltdir, ls_field, time, z, radius, nz, ierr)
      character(len=*), intent(in) :: pltdir, ls_field
      real(dp), intent(out) :: time
      real(dp), allocatable, intent(out) :: z(:), radius(:)
      integer, intent(out) :: nz, ierr

      type(plotfile_t) :: pf
      character(len=64) :: want(1)
      integer :: i, j, nr, ihit
      real(dp) :: rmax, zmax, denom, f
      real(dp), allocatable :: r(:)

      want(1) = ls_field
      call read_plotfile(pltdir, pf, ierr, want)
      if (ierr /= 0) then
         nz = 0
         return
      end if

      nr = pf%nx
      nz = pf%ny
      rmax = pf%prob_hi(1)
      zmax = pf%prob_hi(2)
      time = pf%time

      allocate(r(nr), z(nz), radius(nz))
      do i = 1, nr
         r(i) = (real(i,dp) - 0.5_dp) * rmax / real(nr,dp)
      end do
      do j = 1, nz
         z(j) = (real(j,dp) - 0.5_dp) * zmax / real(nz,dp)
      end do

      radius = -1.0_dp                       ! sentinel: not found
      do j = 1, nz
         ihit = 0
         do i = 1, nr - 1
            if (pf%data(i,j,1) * pf%data(i+1,j,1) < 0.0_dp) then
               ihit = i
               exit
            end if
         end do
         if (ihit > 0) then
            denom = pf%data(ihit,j,1) - pf%data(ihit+1,j,1)
            if (denom /= 0.0_dp) then
               f = pf%data(ihit,j,1) / denom
            else
               f = 0.0_dp
            end if
            radius(j) = r(ihit) + f * (r(ihit+1) - r(ihit))
         end if
      end do

      call free_plotfile(pf)
   end subroutine interface_radius

   !> Amplitude of the `harmonic`-th Fourier component of r(z).
   !>
   !> a = 2 |<(r - <r>) exp(-i 2 pi n z / lambda)>|
   !>
   !> The factor 2 makes `a` the amplitude of a cos(k z) perturbation, matching the dispersion-relation convention.  The mean is removed first so a drifting mean radius (a mass-conservation error, say) cannot leak into the mode amplitude.
   function mode_amplitude(z, radius, wavelength, harmonic) result(a)
      real(dp), intent(in) :: z(:), radius(:), wavelength
      integer, intent(in) :: harmonic
      real(dp) :: a
      integer :: j, n
      real(dp) :: rbar, ang, sr, si

      n = 0
      rbar = 0.0_dp
      do j = 1, size(radius)
         if (radius(j) >= 0.0_dp) then
            n = n + 1
            rbar = rbar + radius(j)
         end if
      end do
      if (n < 4) then
         a = -1.0_dp                          ! not enough valid columns
         return
      end if
      rbar = rbar / real(n,dp)

      sr = 0.0_dp; si = 0.0_dp
      do j = 1, size(radius)
         if (radius(j) >= 0.0_dp) then
            ang = 2.0_dp * PI * real(harmonic,dp) * z(j) / wavelength
            sr = sr + (radius(j) - rbar) * cos(ang)
            si = si - (radius(j) - rbar) * sin(ang)
         end if
      end do
      a = 2.0_dp * sqrt(sr*sr + si*si) / real(n,dp)
   end function mode_amplitude

   !> Interface mode amplitude a(t) over every plotfile in a run.
   !>
   !> `meanr` is returned because a drifting mean radius is the usual signal of a mass-conservation problem, which would invalidate the growth rate.
   subroutine amplitude_series(run_dir, wavelength, ls_field, harmonic, &
                               times, amps, meanr, n, ierr, prefix)
      character(len=*), intent(in) :: run_dir, ls_field
      real(dp), intent(in) :: wavelength
      integer, intent(in) :: harmonic
      real(dp), allocatable, intent(out) :: times(:), amps(:), meanr(:)
      integer, intent(out) :: n, ierr
      character(len=*), intent(in), optional :: prefix

      character(len=MAX_NAME), allocatable :: paths(:)
      real(dp), allocatable :: z(:), radius(:)
      character(len=32) :: pfx
      integer :: i, nz, e, j, cnt
      real(dp) :: t, s

      pfx = 'nddataPLT'
      if (present(prefix)) pfx = prefix

      call find_plotfiles(run_dir, paths, n, trim(pfx))
      if (n == 0) then
         ierr = 1                             ! no plotfiles found
         return
      end if

      allocate(times(n), amps(n), meanr(n))
      do i = 1, n
         call interface_radius(trim(paths(i)), ls_field, t, z, radius, nz, e)
         if (e /= 0) then
            ierr = 100 + e
            return
         end if
         times(i) = t
         amps(i) = mode_amplitude(z, radius, wavelength, harmonic)
         s = 0.0_dp; cnt = 0
         do j = 1, nz
            if (radius(j) >= 0.0_dp) then
               s = s + radius(j); cnt = cnt + 1
            end if
         end do
         meanr(i) = merge(s/real(max(cnt,1),dp), -1.0_dp, cnt > 0)
         deallocate(z, radius)
      end do

      ! Sort by TIME, not by step.  For twin pairing, step number is the right key (matching steps must be subtracted).  For a time series it is not: checkpoint restarts can emit a plotfile whose step number does not follow the time ordering of the surrounding files, so a step-sorted series can be non-monotonic in t.  Fits that select by time value are unaffected, but anything using "the last element" would be wrong.
      call sort_by_time(times, amps, meanr, n)
      ierr = 0
   end subroutine amplitude_series

   !> Insertion-sort three parallel arrays by ascending time.
   pure subroutine sort_by_time(times, amps, meanr, n)
      real(dp), intent(inout) :: times(:), amps(:), meanr(:)
      integer, intent(in) :: n
      integer :: i, j
      real(dp) :: kt, ka, km
      do i = 2, n
         kt = times(i); ka = amps(i); km = meanr(i)
         j = i - 1
         do while (j >= 1)
            if (times(j) <= kt) exit
            times(j+1) = times(j); amps(j+1) = amps(j); meanr(j+1) = meanr(j)
            j = j - 1
         end do
         times(j+1) = kt; amps(j+1) = ka; meanr(j+1) = km
      end do
   end subroutine sort_by_time

   !> Log-linear fit of a(t) over an optional window.
   subroutine fit_growth_rate_amp(times, amps, t_min, t_max, sigma, r2, ierr)
      real(dp), intent(in) :: times(:), amps(:), t_min, t_max
      real(dp), intent(out) :: sigma, r2
      integer, intent(out) :: ierr

      integer :: i, m
      real(dp) :: st, sy, stt, sty, t, y, ybar, ssr, sst, pred, b

      m = 0; st = 0; sy = 0; stt = 0; sty = 0
      do i = 1, size(times)
         if (amps(i) <= 0.0_dp) cycle
         if (times(i) < t_min) cycle
         if (t_max > t_min .and. times(i) > t_max) cycle
         t = times(i); y = log(amps(i))
         m = m + 1; st = st + t; sy = sy + y
         stt = stt + t*t; sty = sty + t*y
      end do
      if (m < 3) then
         ierr = 1; sigma = 0; r2 = 0; return
      end if
      sigma = (real(m,dp)*sty - st*sy) / (real(m,dp)*stt - st*st)
      b = (sy - sigma*st) / real(m,dp)
      ybar = sy / real(m,dp)
      ssr = 0; sst = 0
      do i = 1, size(times)
         if (amps(i) <= 0.0_dp) cycle
         if (times(i) < t_min) cycle
         if (t_max > t_min .and. times(i) > t_max) cycle
         y = log(amps(i)); pred = sigma*times(i) + b
         ssr = ssr + (y-pred)**2
         sst = sst + (y-ybar)**2
      end do
      r2 = merge(1.0_dp - ssr/sst, 0.0_dp, sst > 0.0_dp)
      ierr = 0
   end subroutine fit_growth_rate_amp

   !> Instantaneous growth rate d(ln a)/dt over sliding windows.
   !>
   !> This is the diagnostic that says whether the run has reached asymptotic modal growth.  If sigma_local is still CLIMBING at the end of the record, the flow is still in the startup transient and any single fitted growth rate underestimates the true one - which is exactly what the first real runs of this problem showed.
   subroutine local_growth_rate(times, amps, width, centres, sigmas, nout)
      real(dp), intent(in) :: times(:), amps(:), width
      real(dp), allocatable, intent(out) :: centres(:), sigmas(:)
      integer, intent(out) :: nout

      integer, parameter :: NW = 8
      integer :: i, e
      real(dp) :: lo, hi, c, s, r2

      allocate(centres(NW), sigmas(NW))
      nout = 0
      lo = minval(times) + 0.5_dp*width
      hi = maxval(times) - 0.5_dp*width
      if (hi <= lo) return

      do i = 1, NW
         c = lo + (hi-lo)*real(i-1,dp)/real(NW-1,dp)
         call fit_growth_rate_amp(times, amps, c-0.5_dp*width, &
                                  c+0.5_dp*width, s, r2, e)
         if (e == 0) then
            nout = nout + 1
            centres(nout) = c
            sigmas(nout) = s
         end if
      end do
   end subroutine local_growth_rate

   !> Heuristic: has the local growth rate stopped climbing?
   !>
   !> Compares the last two sliding-window estimates.  True is NECESSARY but not sufficient - also check grid convergence.
   function has_converged(times, amps, width, rel_tol, sigma_last, rel) &
         result(ok)
      real(dp), intent(in) :: times(:), amps(:), width, rel_tol
      real(dp), intent(out) :: sigma_last, rel
      logical :: ok
      real(dp), allocatable :: c(:), s(:)
      integer :: n

      call local_growth_rate(times, amps, width, c, s, n)
      if (n < 2) then
         ok = .false.; sigma_last = -1.0_dp; rel = -1.0_dp
         return
      end if
      sigma_last = s(n)
      rel = abs(s(n) - s(n-1)) / max(abs(s(n)), 1.0e-30_dp)
      ok = (rel <= rel_tol)
   end function has_converged

end module interface_mode_mod
