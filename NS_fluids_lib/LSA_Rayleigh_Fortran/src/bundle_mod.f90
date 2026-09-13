!> Portable figure data: interface profiles and field slices as plain text.
!>
!> An amplitude record is enough to draw the growth curve and the local growth rate, but two of the four figures reported for a case need more than that. The interface figure needs the radius profile r(z) at several instants, and the field figure needs a plotfile variable on the (r,z) plane together with the interface contour laid over it. Both were previously obtainable only by reading the plotfiles, which means keeping the plotfiles; a sweep of sixteen cases leaves tens of gigabytes on the machine that ran it, and transferring that to a workstation merely to redraw a legend is not practical.
!>
!> This module writes a small text format carrying exactly what the figures need. A profile bundle for a typical case is a few hundred kilobytes, so the figure data of a whole sweep travels in a directory small enough to copy over a login session. The format is identical to the one `python/bundle.py` writes in the Python tree, and the plotting scripts of both trees read it through that one reader, so a bundle written here is plotted by the same code that plots a bundle written there.
!>
!> The format is deliberately crude: a keyword line, then free-format reals. It is written here with ordinary formatted output and parsed in Python by splitting on whitespace. Nothing in it depends on a library that might be absent on a cluster.
!>
!> A profile bundle, `profiles.dat`:
!>
!>     # LSA profile bundle v1
!>     WAVELENGTH   8.9760000000000000E+00
!>     NZ           72
!>     NBLOCKS      120
!>     Z
!>       <NZ reals>
!>     BLOCK  1   <time>
!>       <NZ radius values>
!>     ...
!>     END
!>
!> A field bundle, `field_<name>.dat`, carries NR and NZ, then the R and Z axes, then one block of NR*NZ reals per instant with the radial index varying fastest.
!>
!> A radius of -1 marks an axial station at which no interface was found, which happens when the interface leaves the radial extent of the domain. The Python reader turns those into NaN so that a plot leaves a gap rather than drawing a spurious excursion to the axis. That sentinel is exactly what `interface_radius` in interface_mode_mod already returns, so nothing is translated on the way out.
module bundle_mod
   use plotfile_mod, only: dp
   implicit none
   private

   public :: dp
   public :: write_profile_bundle, write_field_bundle

   !> Reals per line. Five keeps the file readable in a pager, and list-directed input ignores line breaks entirely, so the choice is cosmetic.
   integer, parameter :: PER_LINE = 5

   !> Output edit descriptor. The three-digit exponent is not cosmetic. With the default two-digit form, `es24.16`, a value whose exponent needs three digits is written with the exponent letter dropped, so 4.66e-310 emerges as `4.6598606362008991-310`. Fortran reads that back correctly, since list-directed input accepts the form, but no other language does, and the Python reader that draws every figure would reject the file. Forcing `e3` makes the letter unconditional. The width is 25 rather than 24 because a negative mantissa with a three-digit exponent needs one column more than the default form.
   character(len=*), parameter :: REAL_FMT = '(5es25.16e3)'

contains

   !> Write the numbers of a block, five to a line.
   subroutine write_reals(unit_out, values, n)
      integer, intent(in) :: unit_out, n
      real(dp), intent(in) :: values(n)
      integer :: i

      do i = 1, n, PER_LINE
         write(unit_out,REAL_FMT) values(i:min(i + PER_LINE - 1, n))
      end do
   end subroutine write_reals

   !> Write a profile bundle.
   !>
   !> `radius(j,it)` is the interface radius at axial station j and snapshot it, in the layout `interface_radius` produces column by column. `ierr` is non-zero only if the file cannot be opened, since everything else here is a formatted write to a unit already known good.
   subroutine write_profile_bundle(path, wavelength, z, nz, times, radius, &
                                   nt, ierr)
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: wavelength
      integer, intent(in) :: nz, nt
      real(dp), intent(in) :: z(nz), times(nt), radius(nz, nt)
      integer, intent(out) :: ierr

      integer :: unit_out, it

      open(newunit=unit_out, file=trim(path), status='replace', &
           action='write', iostat=ierr)
      if (ierr /= 0) return

      write(unit_out,'(a)') '# LSA profile bundle v1'
      write(unit_out,'(a,es25.16e3)') 'WAVELENGTH ', wavelength
      write(unit_out,'(a,i0)') 'NZ ', nz
      write(unit_out,'(a,i0)') 'NBLOCKS ', nt
      write(unit_out,'(a)') 'Z'
      call write_reals(unit_out, z, nz)
      do it = 1, nt
         write(unit_out,'(a,i0,a,es25.16e3)') 'BLOCK ', it, ' ', times(it)
         call write_reals(unit_out, radius(:, it), nz)
      end do
      write(unit_out,'(a)') 'END'
      close(unit_out)
   end subroutine write_profile_bundle

   !> Write a field bundle.
   !>
   !> `plane(i,j,it)` is the variable at radial index i, axial index j and snapshot it. The radial index varies fastest in the file, which is the same order the Python writer emits, so the two files agree term by term.
   subroutine write_field_bundle(path, field, r, nr, z, nz, times, plane, &
                                 nt, ierr)
      character(len=*), intent(in) :: path, field
      integer, intent(in) :: nr, nz, nt
      real(dp), intent(in) :: r(nr), z(nz), times(nt), plane(nr, nz, nt)
      integer, intent(out) :: ierr

      integer :: unit_out, it, j

      open(newunit=unit_out, file=trim(path), status='replace', &
           action='write', iostat=ierr)
      if (ierr /= 0) return

      write(unit_out,'(a)') '# LSA field bundle v1'
      write(unit_out,'(a,a)') 'FIELD ', trim(field)
      write(unit_out,'(a,i0)') 'NR ', nr
      write(unit_out,'(a,i0)') 'NZ ', nz
      write(unit_out,'(a,i0)') 'NBLOCKS ', nt
      write(unit_out,'(a)') 'R'
      call write_reals(unit_out, r, nr)
      write(unit_out,'(a)') 'Z'
      call write_reals(unit_out, z, nz)
      do it = 1, nt
         write(unit_out,'(a,i0,a,es25.16e3)') 'BLOCK ', it, ' ', times(it)
         ! Written column by column so that the radial index varies fastest across the whole block, matching the Fortran-order ravel the Python writer performs on the same array.
         do j = 1, nz
            call write_reals(unit_out, plane(:, j, it), nr)
         end do
      end do
      write(unit_out,'(a)') 'END'
      close(unit_out)
   end subroutine write_field_bundle

end module bundle_mod
