!> Tests for dispersion_mod.
!>
!> The module reproduces the truncated power series used by the solver's own run2d/rayleigh_capillary_growth.F90, so testing it against that same series would be circular.  Instead I_n(x) is cross-checked against an INDEPENDENT evaluation (numerical integration of the standard integral representation
!>
!> I_n(x) = (1/pi) * integral_0^pi exp(x cos t) cos(n t) dt
!>
!> by the trapezoid rule, which converges spectrally for this smooth periodic integrand).  Agreement to ~1e-14 pins the series independently.
!>
!> The headline value, sigma = 0.34332680120934783 at the single-precision k = 0.7 the original program actually uses, is checked against the number printed by that program when compiled and run.
program test_dispersion
   use testing_mod
   use dispersion_mod
   implicit none

   real(dp) :: k, k_max, s_max, lam, s1, s2, s4
   integer :: i

   call section('modified Bessel I_n against an independent integral')
   do i = 0, 10
      k = real(i, dp) / 10.0_dp
      ! Mixed tolerance: I_1(0) is exactly zero, so a purely relative test against the quadrature's round-off can never pass there.
      call check_close_mixed(mod_bessel_i(k, 0), bessel_i_quad(k, 0), &
                             1.0e-14_dp, 1.0e-13_dp, 'I_0 matches quadrature')
      call check_close_mixed(mod_bessel_i(k, 1), bessel_i_quad(k, 1), &
                             1.0e-14_dp, 1.0e-13_dp, 'I_1 matches quadrature')
   end do

   call section('growth rate against the original Fortran utility')
   ! rayleigh_capillary_growth.F90 stores k in SINGLE precision, so the value it actually evaluates is 0.69999998807907104, and it prints 0.34332680120934783.  Feeding the same value must reproduce it exactly.
   call check_close(growth_rate(0.69999998807907104_dp, 1.0_dp, 1.0_dp, &
      1.0_dp), &
                    0.34332680120934783_dp, 1.0e-15_dp, &
                    'reproduces the shipped utility to 1e-15')

   call section('stability boundary')
   call check(growth_rate_squared(0.999_dp, 1.0_dp, 1.0_dp, 1.0_dp) > 0.0_dp, &
              'unstable just below k r0 = 1')
   call check_close(growth_rate_squared(1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp), &
                    0.0_dp, 1.0e-12_dp, 'marginal at k r0 = 1')
   call check(growth_rate_squared(1.001_dp, 1.0_dp, 1.0_dp, 1.0_dp) < 0.0_dp, &
              'stable just above k r0 = 1')
   call check_close(growth_rate(1.5_dp, 1.0_dp, 1.0_dp, 1.0_dp), 0.0_dp, &
                    0.0_dp, 'growth rate is zero (not NaN) in the stable band')

   call section('most unstable mode')
   call most_unstable_wavenumber(1.0_dp, 1.0_dp, 1.0_dp, k_max, s_max)
   ! The .F90 header states "anecdotedly, k r0 = .7 is the critical value".
   call check_close(k_max, 0.697_dp, 5.0e-3_dp, 'peak at k r0 ~= 0.697')
   lam = most_unstable_wavelength(1.0_dp, 1.0_dp, 1.0_dp)
   ! inputs.growthrate.LSA uses yblob = 9.0 for "wave length = 2 pi r0/.7"
   call check_close(lam, 9.0_dp, 0.1_dp, 'lambda ~= 9.0 r0')
   call check_close_rel(lam, 2.0_dp*PI/0.7_dp, 0.01_dp, &
                        'lambda ~= 2 pi r0 / 0.7')
   call check(s_max >= growth_rate(0.5_dp,1.0_dp,1.0_dp,1.0_dp), &
              'peak growth exceeds k r0 = 0.5')
   call check(s_max >= growth_rate(0.9_dp,1.0_dp,1.0_dp,1.0_dp), &
              'peak growth exceeds k r0 = 0.9')

   call section('scaling with physical parameters')
   s1 = growth_rate(0.7_dp, 1.0_dp, 1.0_dp, 1.0_dp)
   s2 = growth_rate(0.7_dp, 1.0_dp, 4.0_dp, 1.0_dp)
   call check_close_rel(s2, 2.0_dp*s1, 1.0e-12_dp, 'sigma ~ sqrt(surface tension)')
   s2 = growth_rate(0.7_dp, 1.0_dp, 1.0_dp, 4.0_dp)
   call check_close_rel(s2, 0.5_dp*s1, 1.0e-12_dp, 'sigma ~ 1/sqrt(density)')
   ! at fixed k*r0, doubling r0 scales sigma by r0^{-3/2}
   s4 = growth_rate(0.35_dp, 2.0_dp, 1.0_dp, 1.0_dp)
   call check_close_rel(s4, s1*2.0_dp**(-1.5_dp), 1.0e-12_dp, &
                        'sigma ~ r0^(-3/2) at fixed k r0')

   call section('viscous correction')
   call check_close(growth_rate_viscous(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.0_dp), &
                    s1, 1.0e-14_dp, 'reduces to inviscid as mu -> 0')
   call check(growth_rate_viscous(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp) < s1, &
              'viscosity reduces the growth rate')
   call check(growth_rate_viscous(0.7_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp) > &
      0.0_dp, &
              'still unstable at Oh = 0.02')
   ! stable band: must be real and non-positive, never NaN
   s2 = growth_rate_viscous(1.5_dp,1.0_dp,1.0_dp,1.0_dp,0.02_dp)
   call check(s2 == s2, 'viscous rate is not NaN in the stable band')
   call check(s2 <= 0.0_dp, 'viscous rate is non-positive in the stable band')

   call section('Ohnesorge number and box length')
   call check_close(ohnesorge(0.02_dp,1.0_dp,1.0_dp,1.0_dp), 0.02_dp, &
                    1.0e-14_dp, 'Oh = mu for unit rho, sigma, r0')
   call check_close(box_length_for_mode(0.7_dp, 1), 2.0_dp*PI/0.7_dp, &
                    1.0e-12_dp, 'box holds one wavelength')
   call check_close(box_length_for_mode(0.7_dp, 8), 8.0_dp*2.0_dp*PI/0.7_dp, &
                    1.0e-12_dp, 'box holds eight wavelengths')

   call summary('test_dispersion')

contains

   !> I_n(x) by trapezoid quadrature of the integral representation. Independent of the power series under test.
   function bessel_i_quad(x, n) result(val)
      real(dp), intent(in) :: x
      integer, intent(in) :: n
      real(dp) :: val, t, h
      integer :: i, m
      m = 20000
      h = PI / real(m, dp)
      val = 0.0_dp
      do i = 0, m
         t = real(i, dp) * h
         if (i == 0 .or. i == m) then
            val = val + 0.5_dp * exp(x*cos(t)) * cos(real(n,dp)*t)
         else
            val = val + exp(x*cos(t)) * cos(real(n,dp)*t)
         end if
      end do
      val = val * h / PI
   end function bessel_i_quad

end program test_dispersion
