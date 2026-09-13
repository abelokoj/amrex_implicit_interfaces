!> Minimal assertion harness for the Fortran test suite.
!>
!> Deliberately tiny and dependency-free: adding a unit-test framework would complicate the build for no benefit at this size.  Each test program calls `check*` repeatedly and ends with `summary`, which sets a non-zero exit status when anything failed so `make test` reports it.
module testing_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private

   public :: dp
   public :: check, check_close, check_close_rel, check_int, check_str
   public :: check_close_mixed
   public :: summary, section

   integer, save :: n_pass = 0
   integer, save :: n_fail = 0

contains

   subroutine section(name)
      character(len=*), intent(in) :: name
      write(*,'(a)') ''
      write(*,'(a,a)') '== ', trim(name)
   end subroutine section

   subroutine check(cond, name)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: name
      if (cond) then
         n_pass = n_pass + 1
         write(*,'(a,a)') '  ok   ', trim(name)
      else
         n_fail = n_fail + 1
         write(*,'(a,a)') '  FAIL ', trim(name)
      end if
   end subroutine check

   !> Absolute tolerance comparison.
   subroutine check_close(got, want, atol, name)
      real(dp), intent(in) :: got, want, atol
      character(len=*), intent(in) :: name
      if (abs(got - want) <= atol) then
         n_pass = n_pass + 1
         write(*,'(a,a)') '  ok   ', trim(name)
      else
         n_fail = n_fail + 1
         write(*,'(a,a)') '  FAIL ', trim(name)
         write(*,'(a,es22.14,a,es22.14,a,es10.2)') '       got ', got, &
            '  want ', want, '  atol ', atol
      end if
   end subroutine check_close

   !> Relative tolerance comparison.
   subroutine check_close_rel(got, want, rtol, name)
      real(dp), intent(in) :: got, want, rtol
      character(len=*), intent(in) :: name
      real(dp) :: denom
      denom = max(abs(want), tiny(1.0_dp))
      if (abs(got - want) / denom <= rtol) then
         n_pass = n_pass + 1
         write(*,'(a,a)') '  ok   ', trim(name)
      else
         n_fail = n_fail + 1
         write(*,'(a,a)') '  FAIL ', trim(name)
         write(*,'(a,es22.14,a,es22.14,a,es10.2)') '       got ', got, &
            '  want ', want, '  rtol ', rtol
      end if
   end subroutine check_close_rel

   !> Passes if EITHER the absolute or the relative tolerance is met.
   !>
   !> Needed whenever the expected value may be exactly zero: a pure relative test against 0 can never pass, however small the actual error, because any round-off in the reference is infinitely large in relative terms.
   subroutine check_close_mixed(got, want, atol, rtol, name)
      real(dp), intent(in) :: got, want, atol, rtol
      character(len=*), intent(in) :: name
      real(dp) :: adiff
      adiff = abs(got - want)
      if (adiff <= atol .or. adiff <= rtol * abs(want)) then
         n_pass = n_pass + 1
         write(*,'(a,a)') '  ok   ', trim(name)
      else
         n_fail = n_fail + 1
         write(*,'(a,a)') '  FAIL ', trim(name)
         write(*,'(a,es22.14,a,es22.14)') '       got ', got, '  want ', want
      end if
   end subroutine check_close_mixed

   subroutine check_int(got, want, name)
      integer, intent(in) :: got, want
      character(len=*), intent(in) :: name
      if (got == want) then
         n_pass = n_pass + 1
         write(*,'(a,a)') '  ok   ', trim(name)
      else
         n_fail = n_fail + 1
         write(*,'(a,a)') '  FAIL ', trim(name)
         write(*,'(a,i0,a,i0)') '       got ', got, '  want ', want
      end if
   end subroutine check_int

   subroutine check_str(got, want, name)
      character(len=*), intent(in) :: got, want, name
      if (trim(got) == trim(want)) then
         n_pass = n_pass + 1
         write(*,'(a,a)') '  ok   ', trim(name)
      else
         n_fail = n_fail + 1
         write(*,'(a,a)') '  FAIL ', trim(name)
         write(*,'(a,a,a,a)') '       got [', trim(got), '] want [', trim(want)
      end if
   end subroutine check_str

   !> Print totals and stop with a non-zero status if anything failed.
   subroutine summary(title)
      character(len=*), intent(in) :: title
      write(*,'(a)') ''
      write(*,'(a,a,a,i0,a,i0,a)') '--- ', trim(title), ': ', n_pass, &
         ' passed, ', n_fail, ' failed'
      if (n_fail > 0) then
         stop 1
      end if
   end subroutine summary

end module testing_mod
