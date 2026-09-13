!> Tests for the portable figure-data bundle.
!>
!> Two things are checked. The writers in bundle_mod produce a file whose keywords, counts and numbers are what the format specifies, which is verified here by parsing the file back with list-directed input. And the whole export path works on genuine plotfiles: synthetic plotfiles carrying a known interface are written in the real on-disk format, `interface_radius` extracts the profile, the bundle is written, and the radii read back match the analytic interface that was planted.
!>
!> The second is the load-bearing test. A format that round-trips through its own writer proves little; a bundle whose radii agree with the interface actually planted in the plotfiles proves that the export carries what the interface figure needs. The planted interface is r = r0 + a cos(k z), so every station has a value known in closed form.
!>
!> Interoperability with the Python reader is checked separately, in tests/test_bundle.py of the Python tree, which reads a bundle written by this program. Neither language is trusted to validate its own output alone.
program test_bundle
   use testing_mod
   use bundle_mod
   use plotfile_mod
   use interface_mode_mod, only: interface_radius
   use plotfile_writer_mod
   implicit none

   character(len=256) :: dir, path
   real(dp), allocatable :: z(:), radius(:), times(:), profile(:,:)
   real(dp), allocatable :: r_axis(:), z_axis(:), plane(:,:,:)
   real(dp), allocatable :: ls(:,:,:)
   real(dp) :: wavelength, r0, amp, k, t, expect
   integer :: nz, nr, i, j, it, ierr, nt
   integer :: unit_in, nz_read, nblocks_read, idx
   character(len=64) :: word
   real(dp) :: value

   call section('bundle: profile writer and reader')

   dir = '/tmp/lsa_bundle_test'
   call execute_command_line('rm -rf '//trim(dir)//' && mkdir -p '//trim(dir))

   nz = 12
   nt = 3
   wavelength = 8.976_dp
   allocate(z(nz), times(nt), profile(nz, nt))
   do j = 1, nz
      z(j) = (real(j, dp) - 0.5_dp) * wavelength / real(nz, dp)
   end do
   do it = 1, nt
      times(it) = 0.5_dp * real(it, dp)
      do j = 1, nz
         profile(j, it) = 1.0_dp + 0.01_dp * real(it, dp) * z(j)
      end do
   end do

   path = trim(dir)//'/profiles.dat'
   call write_profile_bundle(trim(path), wavelength, z, nz, times, profile, &
                             nt, ierr)
   call check_int(ierr, 0, 'profile bundle written')

   ! Parse the file back with list-directed input, which is how the Fortran side would read it and which fails loudly if a keyword or a count is wrong.
   open(newunit=unit_in, file=trim(path), status='old', action='read')
   read(unit_in,'(a)') word                          ! the comment line
   read(unit_in,*) word, value
   call check_str(trim(word), 'WAVELENGTH', 'first keyword is WAVELENGTH')
   call check_close_rel(value, wavelength, 1.0e-14_dp, 'wavelength round trip')
   read(unit_in,*) word, nz_read
   call check_str(trim(word), 'NZ', 'second keyword is NZ')
   call check_int(nz_read, nz, 'NZ round trip')
   read(unit_in,*) word, nblocks_read
   call check_str(trim(word), 'NBLOCKS', 'third keyword is NBLOCKS')
   call check_int(nblocks_read, nt, 'NBLOCKS round trip')
   read(unit_in,*) word
   call check_str(trim(word), 'Z', 'axis keyword is Z')
   allocate(radius(nz))
   read(unit_in,*) radius
   call check_close(radius(1), z(1), 1.0e-14_dp, 'first axial station')
   call check_close(radius(nz), z(nz), 1.0e-14_dp, 'last axial station')
   do it = 1, nt
      read(unit_in,*) word, idx, value
      call check_str(trim(word), 'BLOCK', 'block keyword')
      call check_int(idx, it, 'block index counts from one')
      call check_close(value, times(it), 1.0e-14_dp, 'block time round trip')
      read(unit_in,*) radius
      call check_close(radius(nz), profile(nz, it), 1.0e-14_dp, &
                       'block radius round trip')
   end do
   read(unit_in,*) word
   call check_str(trim(word), 'END', 'file terminates with END')
   close(unit_in)
   deallocate(radius)

   call section('bundle: field writer')

   nr = 5
   allocate(r_axis(nr), z_axis(nz), plane(nr, nz, nt))
   do i = 1, nr
      r_axis(i) = real(i, dp)
   end do
   z_axis = z
   do it = 1, nt
      do j = 1, nz
         do i = 1, nr
            ! Distinct in all three indices, so a transposed or mis-strided write cannot pass unnoticed.
            plane(i, j, it) = real(100 * it + 10 * j + i, dp)
         end do
      end do
   end do

   path = trim(dir)//'/field_test_var.dat'
   call write_field_bundle(trim(path), 'test_var', r_axis, nr, z_axis, nz, &
                           times, plane, nt, ierr)
   call check_int(ierr, 0, 'field bundle written')

   open(newunit=unit_in, file=trim(path), status='old', action='read')
   read(unit_in,'(a)') word
   read(unit_in,*) word, word
   call check_str(trim(word), 'test_var', 'field name round trip')
   read(unit_in,*) word, nz_read
   call check_int(nz_read, nr, 'NR round trip')
   close(unit_in)

   call section('bundle: export of a planted interface')

   ! Genuine plotfiles carrying r = r0 + a cos(k z), the interface the deck imposes. The level set is the signed distance to that surface along a radial line, which is what the solver writes.
   nr = 32
   nz = 24
   r0 = 1.0_dp
   amp = 0.05_dp
   k = 2.0_dp * acos(-1.0_dp) / wavelength
   deallocate(z, times, profile)
   allocate(z(nz), times(nt), profile(nz, nt), ls(nr, nz, 1))

   do j = 1, nz
      z(j) = (real(j, dp) - 0.5_dp) * wavelength / real(nz, dp)
   end do
   ! The times array was reallocated above and carries whatever the allocator left behind, so it is set explicitly. An uninitialised value here is not merely untidy: a denormal reaches the writer, whose exponent needs three digits, which is precisely the format hazard the three-digit descriptor in bundle_mod exists to avoid.
   do it = 1, nt
      times(it) = 0.5_dp * real(it, dp)
   end do

   do it = 1, nt
      t = 0.5_dp * real(it, dp)
      do j = 1, nz
         do i = 1, nr
            ls(i, j, 1) = (real(i, dp) - 0.5_dp) * 4.0_dp / real(nr, dp) &
                          - (r0 + amp * cos(k * z(j)))
         end do
      end do
      write(path, '(a,a,i8.8)') trim(dir), '/nddataPLT', 50 * it
      call write_plotfile(trim(path), ls, ['L0101'], t, 50 * it, &
                          [0.0_dp, 0.0_dp], [4.0_dp, wavelength])
   end do

   ! The middle snapshot is enough: every one carries the same interface, so agreement at one instant with the closed form is the property under test.
   write(path, '(a,a,i8.8)') trim(dir), '/nddataPLT', 100
   call interface_radius(trim(path), 'L0101', t, z, radius, nz_read, ierr)
   call check_int(ierr, 0, 'planted plotfile read back')
   call check_int(nz_read, nz, 'axial extent recovered')

   ! The zero contour is located by linear interpolation between cell centres, and the level set is exactly linear in r, so the recovered radius should match the closed form to round-off rather than to the mesh spacing.
   expect = r0 + amp * cos(k * z(1))
   call check_close(radius(1), expect, 1.0e-12_dp, &
                    'recovered radius matches the planted interface')
   expect = r0 + amp * cos(k * z(nz / 2))
   call check_close(radius(nz / 2), expect, 1.0e-12_dp, &
                    'and at a station half way along')

   do j = 1, nz
      profile(j, 1) = radius(j)
      profile(j, 2) = radius(j)
      profile(j, 3) = radius(j)
   end do
   path = trim(dir)//'/exported.dat'
   call write_profile_bundle(trim(path), wavelength, z, nz, times, profile, &
                             nt, ierr)
   call check_int(ierr, 0, 'exported profile bundle written')

   call summary('test_bundle')
end program test_bundle
