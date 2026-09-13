!> Export the figure data of a calculation as plain text, once.
!>
!> `dump_amplitude` writes the amplitude record, which carries the growth curve and the local growth rate. Two further figures need more: the interface figure needs the radius profile r(z) at several instants, and the field figure needs a plotfile variable on the (r,z) plane with the interface laid over it. This program extracts both and writes them in the portable format defined in bundle_mod, so that every figure of a case can afterwards be drawn from a few hundred kilobytes of text rather than from the plotfiles.
!>
!> The separation matters because the plotfiles usually stay where the calculation ran. A sweep leaves tens of gigabytes on a cluster that carries no matplotlib, so the figures are drawn on a workstation; moving the plotfiles there is impractical, while moving a bundle is a single copy. Running this program once, on the machine holding the plotfiles, is what makes that possible, and `scripts/run_all.sh` invokes it for every case.
!>
!> Profiles are exported for every snapshot, since the file remains small. Field slices are exported only at the requested number of instants, spread evenly across the record, because a slice is far larger than a profile and only a few panels are ever drawn.
!>
!>     ./bin/dump_profiles <rundir> <wavelength> <outdir> [field] [npanel]
!>
!> The bundle is read by `plotting/plot_results.py --bundle <outdir>`, which is the same reader the Python tree uses, so a bundle from either implementation plots identically.
program dump_profiles
   use bundle_mod
   use plotfile_mod
   use interface_mode_mod, only: interface_radius
   implicit none

   character(len=512) :: rundir, lambda_text, outdir, field, npanel_text
   character(len=MAX_NAME), allocatable :: paths(:)
   character(len=64) :: want(1)
   real(dp), allocatable :: z(:), radius(:), times(:), profile(:,:)
   real(dp), allocatable :: r_axis(:), z_axis(:), plane(:,:,:)
   real(dp), allocatable :: panel_times(:)
   integer, allocatable :: order(:), chosen(:)
   type(plotfile_t) :: pf
   real(dp) :: wavelength, t, target_t
   integer :: n, nz, nr, i, j, it, ierr, npanel, best, nchosen
   logical :: have_field

   if (command_argument_count() < 3) then
      write(*,'(a)') 'usage: dump_profiles <rundir> <wavelength> <outdir> ' // &
                     '[field] [npanel]'
      stop 1
   end if
   call get_command_argument(1, rundir)
   call get_command_argument(2, lambda_text)
   call get_command_argument(3, outdir)
   read(lambda_text, *) wavelength

   have_field = command_argument_count() >= 4
   field = ''
   if (have_field) call get_command_argument(4, field)

   npanel = 3
   if (command_argument_count() >= 5) then
      call get_command_argument(5, npanel_text)
      read(npanel_text, *) npanel
   end if
   if (npanel < 1) npanel = 1

   call find_plotfiles(trim(rundir), paths, n, 'nddataPLT')
   if (n == 0) then
      write(*,'(a)') 'ERROR: no plotfiles found; the solver writes ' // &
                     'nddataPLT* directories'
      stop 1
   end if

   ! ---- interface profiles, every snapshot ---------------------------
   ! The first plotfile establishes the axial grid, which every snapshot of a run shares because the mesh is fixed.
   call interface_radius(trim(paths(1)), 'L0101', t, z, radius, nz, ierr)
   if (ierr /= 0) then
      write(*,'(a,i0)') 'ERROR reading the first plotfile, ierr=', ierr
      stop 1
   end if

   allocate(times(n), profile(nz, n))
   times(1) = t
   profile(:, 1) = radius
   do i = 2, n
      call interface_radius(trim(paths(i)), 'L0101', t, z, radius, nz, ierr)
      if (ierr /= 0) then
         write(*,'(a,a)') 'ERROR reading ', trim(paths(i))
         stop 1
      end if
      times(i) = t
      profile(:, i) = radius
   end do

   ! A restart can emit a plotfile whose step number does not follow the time ordering of its neighbours, so the record is sorted by time and not by name.
   call sort_index(times, n, order)
   times = times(order)
   profile = profile(:, order)
   paths = paths(order)

   call write_profile_bundle(trim(outdir)//'/profiles.dat', wavelength, &
                             z, nz, times, profile, n, ierr)
   if (ierr /= 0) then
      write(*,'(a)') 'ERROR: cannot write '//trim(outdir)//'/profiles.dat'
      stop 1
   end if
   write(*,'(a,i0,a,a)') 'wrote ', n, ' interface profiles to ', &
                         trim(outdir)//'/profiles.dat'

   ! ---- field slices, a few instants ---------------------------------
   if (.not. have_field) stop 0

   nchosen = min(npanel, n)
   allocate(chosen(nchosen), panel_times(nchosen))
   do i = 1, nchosen
      if (nchosen == 1) then
         target_t = times(n)
      else
         target_t = times(1) + (times(n) - times(1)) * &
                    real(i - 1, dp) / real(nchosen - 1, dp)
      end if
      best = 1
      do j = 2, n
         if (abs(times(j) - target_t) < abs(times(best) - target_t)) best = j
      end do
      chosen(i) = best
      panel_times(i) = times(best)
   end do

   want(1) = trim(field)
   call read_plotfile(trim(paths(chosen(1))), pf, ierr, want)
   if (ierr /= 0) then
      write(*,'(a,a)') 'ERROR: cannot read field ', trim(field)
      stop 1
   end if
   nr = pf%nx
   allocate(r_axis(nr), z_axis(pf%ny), plane(nr, pf%ny, nchosen))
   do i = 1, nr
      r_axis(i) = (real(i, dp) - 0.5_dp) * pf%prob_hi(1) / real(nr, dp)
   end do
   do j = 1, pf%ny
      z_axis(j) = (real(j, dp) - 0.5_dp) * pf%prob_hi(2) / real(pf%ny, dp)
   end do
   plane(:, :, 1) = pf%data(:, :, 1)
   call free_plotfile(pf)

   do it = 2, nchosen
      call read_plotfile(trim(paths(chosen(it))), pf, ierr, want)
      if (ierr /= 0) then
         write(*,'(a,a)') 'ERROR: cannot read ', trim(paths(chosen(it)))
         stop 1
      end if
      plane(:, :, it) = pf%data(:, :, 1)
      call free_plotfile(pf)
   end do

   call write_field_bundle(trim(outdir)//'/field_'//trim(field)//'.dat', &
                           trim(field), r_axis, nr, z_axis, size(z_axis), &
                           panel_times, plane, nchosen, ierr)
   if (ierr /= 0) then
      write(*,'(a)') 'ERROR: cannot write the field bundle'
      stop 1
   end if
   write(*,'(a,i0,a,a,a,a)') 'wrote ', nchosen, ' slices of ', trim(field), &
                             ' to ', trim(outdir)//'/field_'// &
                             trim(field)//'.dat'

contains

   !> Indices that sort `values` into increasing order, by insertion.
   !>
   !> A record holds at most a few hundred snapshots, so the quadratic cost is irrelevant and the absence of a dependency is worth more than the asymptotics.
   subroutine sort_index(values, m, idx)
      real(dp), intent(in) :: values(:)
      integer, intent(in) :: m
      integer, allocatable, intent(out) :: idx(:)
      integer :: a, b, key

      allocate(idx(m))
      do a = 1, m
         idx(a) = a
      end do
      do a = 2, m
         key = idx(a)
         b = a - 1
         do while (b >= 1)
            if (values(idx(b)) <= values(key)) exit
            idx(b + 1) = idx(b)
            b = b - 1
         end do
         idx(b + 1) = key
      end do
   end subroutine sort_index

end program dump_profiles
