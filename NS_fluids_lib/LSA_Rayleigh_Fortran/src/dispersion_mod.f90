!> Exact linear dispersion relation for the Rayleigh-Plateau instability.
!>
!> Fortran counterpart of python/dispersion.py, and a generalisation of the solver's own utility NS_fluids_lib/run2d/rayleigh_capillary_growth.F90, which implements the inviscid Rayleigh (1878) result quoted in Popinet (2009) section 6.5.
!>
!> For an axisymmetric liquid column of radius r0, density rho and surface tension sigma, perturbed as r(z,t) = r0 + eps exp(s t) cos(k z),
!>
!> s^2 = (sigma/(rho r0^3)) * x * I1(x)/I0(x) * (1 - x^2),   x = k r0
!>
!> unstable for x < 1, peaking at x ~= 0.697 (lambda ~= 9.02 r0).
!>
!> The modified Bessel functions use the same truncated power series as the original .F90 so this module reproduces it bit for bit; test_dispersion checks that claim against an independent continued-fraction evaluation.
module dispersion_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   private

   public :: dp
   public :: mod_bessel_i, growth_rate, growth_rate_squared
   public :: growth_rate_viscous, ohnesorge
   public :: most_unstable_wavenumber, most_unstable_wavelength
   public :: box_length_for_mode

   real(dp), parameter, public :: PI = 3.141592653589793238462643383279502884_dp

contains

   !> Modified Bessel function of the first kind, I_n(x), n >= 0.
   !>
   !> Truncated power series, matching mod_bessel_first_kind in rayleigh_capillary_growth.F90: I_n(x) = sum_{m=0}^{M} (x/2)^{n+2m} / (m! (m+n)!) 100 terms is far more than needed on 0 <= x <= 1, the only range that matters for the unstable band.
   pure function mod_bessel_i(x, n) result(val)
      real(dp), intent(in) :: x
      integer,  intent(in) :: n
      real(dp) :: val
      real(dp) :: m_fact, mpn_fact, xpower
      integer  :: i

      m_fact = 1.0_dp
      mpn_fact = 1.0_dp
      do i = 1, n
         mpn_fact = mpn_fact * real(i, dp)
      end do

      xpower = (0.5_dp * x) ** n
      val = 0.0_dp
      do i = 0, 100
         if (i > 0) then
            m_fact = m_fact * real(i, dp)
            mpn_fact = mpn_fact * real(i + n, dp)
            xpower = xpower * x * x / 4.0_dp
         end if
         val = val + xpower / (m_fact * mpn_fact)
      end do
   end function mod_bessel_i

   !> s^2(k); negative in the stable band x > 1 (oscillatory, not growing).
   pure function growth_rate_squared(k, r0, sigma, rho) result(s2)
      real(dp), intent(in) :: k, r0, sigma, rho
      real(dp) :: s2, x, i0, i1

      x = k * r0
      if (x <= 0.0_dp) then
         s2 = 0.0_dp
         return
      end if
      i0 = mod_bessel_i(x, 0)
      i1 = mod_bessel_i(x, 1)
      s2 = sigma * i1 * x * (1.0_dp - x * x) / (rho * r0**3 * i0)
   end function growth_rate_squared

   !> Inviscid growth rate; zero rather than imaginary in the stable band.
   pure function growth_rate(k, r0, sigma, rho) result(s)
      real(dp), intent(in) :: k, r0, sigma, rho
      real(dp) :: s, s2

      s2 = growth_rate_squared(k, r0, sigma, rho)
      s = merge(sqrt(s2), 0.0_dp, s2 > 0.0_dp)
   end function growth_rate

   !> Ohnesorge number Oh = mu / sqrt(rho sigma r0).
   pure function ohnesorge(mu, rho, sigma, r0) result(oh)
      real(dp), intent(in) :: mu, rho, sigma, r0
      real(dp) :: oh
      oh = mu / sqrt(rho * sigma * r0)
   end function ohnesorge

   !> Growth rate with the leading viscous correction.
   !>
   !> Solves  s^2 + 3 nu k^2 s - s_inviscid^2 = 0  for the positive root (Chandrasekhar, Hydrodynamic and Hydromagnetic Stability, ch. XII). Reduces to the inviscid result as mu -> 0.  In the stable band the discriminant can go negative (a decaying oscillation); we report only the real growth part, so it is clamped at zero.
   pure function growth_rate_viscous(k, r0, sigma, rho, mu) result(s)
      real(dp), intent(in) :: k, r0, sigma, rho, mu
      real(dp) :: s, nu, a, s2, disc

      nu = mu / rho
      s2 = growth_rate_squared(k, r0, sigma, rho)
      a = 1.5_dp * nu * k * k
      disc = max(a * a + s2, 0.0_dp)
      s = -a + sqrt(disc)
   end function growth_rate_viscous

   !> Fastest-growing wavenumber, by dense sampling plus parabolic refinement.
   subroutine most_unstable_wavenumber(r0, sigma, rho, k_max, s_max)
      real(dp), intent(in)  :: r0, sigma, rho
      real(dp), intent(out) :: k_max, s_max
      integer,  parameter   :: n = 200001
      real(dp) :: x, k, s, dk, y0, y1, y2, denom, delta
      real(dp) :: kbest, sbest, kprev, sprev, kprev2, sprev2
      integer  :: i

      kbest = 0.0_dp; sbest = -1.0_dp
      kprev = 0.0_dp; sprev = 0.0_dp
      kprev2 = 0.0_dp; sprev2 = 0.0_dp
      y0 = 0.0_dp; y1 = 0.0_dp; y2 = 0.0_dp
      dk = (1.0_dp / r0) / real(n - 1, dp)

      do i = 1, n
         x = real(i - 1, dp) / real(n - 1, dp)
         if (x <= 0.0_dp) x = 1.0e-9_dp
         k = x / r0
         s = growth_rate(k, r0, sigma, rho)
         if (s > sbest) then
            sbest = s
            kbest = k
            ! remember the bracketing samples for the parabolic refinement
            y0 = sprev2; y1 = sprev; y2 = s
         end if
         kprev2 = kprev; sprev2 = sprev
         kprev = k;      sprev = s
      end do

      ! Parabolic vertex through the three samples around the discrete max. y1 here is the sample just before the max, so shift accordingly.
      y0 = growth_rate(max(kbest - dk, 1.0e-12_dp), r0, sigma, rho)
      y1 = sbest
      y2 = growth_rate(min(kbest + dk, 1.0_dp / r0), r0, sigma, rho)
      denom = y0 - 2.0_dp * y1 + y2
      if (abs(denom) > 0.0_dp) then
         delta = 0.5_dp * (y0 - y2) / denom
      else
         delta = 0.0_dp
      end if

      k_max = kbest + delta * dk
      s_max = growth_rate(k_max, r0, sigma, rho)
   end subroutine most_unstable_wavenumber

   !> Wavelength of the fastest-growing mode.
   function most_unstable_wavelength(r0, sigma, rho) result(lam)
      real(dp), intent(in) :: r0, sigma, rho
      real(dp) :: lam, k_max, s_max
      call most_unstable_wavenumber(r0, sigma, rho, k_max, s_max)
      lam = 2.0_dp * PI / k_max
   end function most_unstable_wavelength

   !> Axial box length holding exactly n wavelengths of wavenumber k.
   pure function box_length_for_mode(k, n) result(L)
      real(dp), intent(in) :: k
      integer,  intent(in) :: n
      real(dp) :: L
      L = 2.0_dp * PI * real(n, dp) / k
   end function box_length_for_mode

end module dispersion_mod
