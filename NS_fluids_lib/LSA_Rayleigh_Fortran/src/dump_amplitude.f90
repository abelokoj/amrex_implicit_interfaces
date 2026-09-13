!> Write the interface-mode amplitude record of a calculation as plain text.
!>
!> The record is the quantity the dispersion relation predicts, namely the amplitude of the k-Fourier component of the interface radius; see interface_mode_mod for why a field norm is unsuitable for an interfacial instability. Output is three columns: time, mode amplitude, and mean interface radius. The last of these is included because a drifting mean signals a mass-conservation error, which would invalidate the growth rate.
!>
!> Writing plain text rather than plotting directly keeps this implementation free of any Python dependency; the figure scripts in plotting/ read these files, and they may equally be plotted elsewhere.
!>
!> ./bin/dump_amplitude <rundir> <wavelength> <outfile>
program dump_amplitude
   use interface_mode_mod
   implicit none

   character(len=512) :: rundir, lambda_text, outfile
   real(dp), allocatable :: times(:), amps(:), mean_radius(:)
   real(dp) :: wavelength
   integer :: n, ierr, i, unit_out

   if (command_argument_count() < 3) then
      write(*,'(a)') 'usage: dump_amplitude <rundir> <wavelength> <outfile>'
      stop 1
   end if
   call get_command_argument(1, rundir)
   call get_command_argument(2, lambda_text)
   call get_command_argument(3, outfile)
   read(lambda_text, *) wavelength

   ! The level set of the first material, L0101, carries the interface; the first harmonic is the mode imposed by the deck.
   call amplitude_series(trim(rundir), wavelength, 'L0101', 1, &
                         times, amps, mean_radius, n, ierr)
   if (ierr /= 0) then
      write(*,'(a,i0)') 'ERROR reading the run, ierr=', ierr
      if (ierr == 1) write(*,'(a)') &
         '  no plotfiles found; the solver writes nddataPLT* directories'
      stop 1
   end if

   open(newunit=unit_out, file=trim(outfile), status='replace', &
        action='write')
   write(unit_out,'(a)') '# time  amplitude  mean_radius'
   do i = 1, n
      write(unit_out,'(3es24.16)') times(i), amps(i), mean_radius(i)
   end do
   close(unit_out)

   write(*,'(a,i0,a,a)') 'wrote ', n, ' rows to ', trim(outfile)
end program dump_amplitude
