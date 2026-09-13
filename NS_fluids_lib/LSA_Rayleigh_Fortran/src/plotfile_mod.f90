!> Native AMReX plotfile reader.
!>
!> The Python implementation delegates plotfile parsing to yt.  There is no equivalent Fortran library, so this module parses the format directly. That makes it the highest-risk component here, and it is consequently the most heavily tested: tests/test_plotfile.f90 writes plotfiles in the real on-disk format and reads them back, and the same synthetic files are cross-checked against yt in the Python tree.
!>
!> Format (AMReX WriteMultiLevelPlotfile, as produced when the deck sets ns.visual_nddata_format = 1):
!>
!> pltNNNNN/Header          global metadata pltNNNNN/Level_0/Cell_H  MultiFab header: boxes and FAB offsets pltNNNNN/Level_0/Cell_D_XXXXX   FABs: ASCII header line + raw doubles
!>
!> Global Header layout: line 1     version string, "HyperCLaw-V1.1" line 2     ncomp next ncomp lines   variable names then       dimensionality, time, finest_level,
!>              prob_lo, prob_hi, ref_ratios, index space,
!>              level_steps, cell sizes, coord_sys, bwidth
!> then per level: "lev ngrids time", level_steps,
!>              per grid 2 reals per dimension, then the level directory
!>
!> Cell_H layout: version, how, ncomp, nghost, "(N 0", N box specs, ")", N, N lines "FabOnDisk: <file> <offset>", then min/max blocks
!>
!> Each FAB begins with an ASCII line such as FAB ((8, (64 11 52 0 1 12 0 1023)),(8, (8 7 6 5 4 3 2 1)))((0,0) (7,11) (0,0)) 2 The second byte-order block gives endianness: descending "8 7 ... 1" is little endian, ascending "1 2 ... 8" is big endian.  Components follow the header in Fortran (column-major) order, one after another, no ghost cells.
!>
!> Restrictions, matching the recommended LSA configuration (amr.max_level = 0): single level, 2-D, double precision.  Multiple grids per level ARE supported, so the reader works for any MPI decomposition.
module plotfile_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64, int64
   implicit none
   private

   public :: dp
   public :: plotfile_t, read_plotfile, free_plotfile
   public :: find_plotfiles, plotfile_step, compatible

   integer, parameter, public :: MAX_NAME = 512

   !> One plotfile, flattened onto a uniform level-0 grid.
   type :: plotfile_t
      character(len=MAX_NAME) :: path = ''
      real(dp) :: time = 0.0_dp
      integer  :: step = -1
      integer  :: ncomp = 0
      integer  :: nx = 0, ny = 0
      real(dp) :: prob_lo(2) = 0.0_dp, prob_hi(2) = 0.0_dp
      character(len=64), allocatable :: names(:)
      !> data(i,j,c), i = 1..nx, j = 1..ny, c = 1..ncomp
      real(dp), allocatable :: data(:,:,:)
   end type plotfile_t

contains

   subroutine free_plotfile(pf)
      type(plotfile_t), intent(inout) :: pf
      if (allocated(pf%names)) deallocate(pf%names)
      if (allocated(pf%data))  deallocate(pf%data)
      pf%ncomp = 0; pf%nx = 0; pf%ny = 0
   end subroutine free_plotfile

   !> Step number encoded in a plotfile name, e.g. "plt00120" -> 120. Returns -1 if there are no trailing digits.
   function plotfile_step(path) result(step)
      character(len=*), intent(in) :: path
      integer :: step
      integer :: i, j, n, base
      character(len=MAX_NAME) :: nm

      ! strip any directory component
      nm = path
      base = 0
      do i = len_trim(path), 1, -1
         if (path(i:i) == '/') then
            base = i
            exit
         end if
      end do
      if (base > 0) nm = path(base+1:)

      n = len_trim(nm)
      j = n
      do while (j >= 1)
         if (nm(j:j) < '0' .or. nm(j:j) > '9') exit
         j = j - 1
      end do
      if (j == n) then
         step = -1
      else
         read(nm(j+1:n), *) step
      end if
   end function plotfile_step

   !> List plotfile directories under `dir`, sorted by step number.
   !>
   !> Uses a directory listing captured to a temporary file, because standard Fortran has no directory-traversal intrinsic.
   subroutine find_plotfiles(dir, paths, nfound, prefix)
      character(len=*), intent(in) :: dir
      character(len=MAX_NAME), allocatable, intent(out) :: paths(:)
      integer, intent(out) :: nfound
      character(len=*), intent(in), optional :: prefix

      character(len=16) :: pfx
      character(len=MAX_NAME) :: tmpfile, line
      character(len=MAX_NAME), allocatable :: tmp(:)
      integer, allocatable :: steps(:), order(:)
      integer :: u, ios, n, i, pid
      character(len=32) :: pidstr

      ! Default confirmed by running the solver: with ns.visual_nddata_format=1 it writes nddataPLT00000023, NOT plt00023, and it does not honour amr.plot_file.
      pfx = 'nddataPLT'
      if (present(prefix)) pfx = prefix

      call get_pid(pid)
      write(pidstr, '(i0)') pid
      tmpfile = '/tmp/.lsa_plt_list_' // trim(pidstr)

      ! -d lists directories only; 2>/dev/null suppresses "no match"
      call execute_command_line( &
         'ls -d ' // trim(dir) // '/' // trim(pfx) // '[0-9]* ' // &
         '2>/dev/null > ' // trim(tmpfile))

      n = 0
      open(newunit=u, file=trim(tmpfile), status='old', action='read', &
           iostat=ios)
      if (ios /= 0) then
         nfound = 0
         allocate(paths(0))
         return
      end if
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) > 0) n = n + 1
      end do
      rewind(u)

      allocate(tmp(max(n,1)), steps(max(n,1)))
      i = 0
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         i = i + 1
         tmp(i) = line
         steps(i) = plotfile_step(line)
      end do
      close(u, status='delete')

      allocate(order(max(n,1)))
      call argsort_int_asc(steps(1:n), order(1:n))
      allocate(paths(max(n,1)))
      do i = 1, n
         paths(i) = tmp(order(i))
      end do
      nfound = n
   end subroutine find_plotfiles

   subroutine get_pid(pid)
      integer, intent(out) :: pid
      integer :: u, ios
      real :: r
      ! getpid() is a GNU extension; fall back on a random tag so that concurrent runs do not collide on the temporary listing file.
      call random_number(r)
      pid = int(r * 1.0e6) + 1
      u = 0; ios = 0
   end subroutine get_pid

   pure subroutine argsort_int_asc(v, idx)
      integer, intent(in)  :: v(:)
      integer, intent(out) :: idx(:)
      integer :: i, j, key
      do i = 1, size(v)
         idx(i) = i
      end do
      do i = 2, size(v)
         key = idx(i)
         j = i - 1
         do while (j >= 1)
            if (v(idx(j)) <= v(key)) exit
            idx(j+1) = idx(j)
            j = j - 1
         end do
         idx(j+1) = key
      end do
   end subroutine argsort_int_asc

   !> Read a plotfile.  ierr is 0 on success.
   subroutine read_plotfile(path, pf, ierr, want_names)
      character(len=*), intent(in) :: path
      type(plotfile_t), intent(out) :: pf
      integer, intent(out) :: ierr
      character(len=*), intent(in), optional :: want_names(:)

      character(len=MAX_NAME) :: line, leveldir, fabfile
      integer :: u, ios, i, j, ncomp, ndim, finest, ngrids
      integer :: lo(2), hi(2), nx, ny
      real(dp) :: t
      integer :: ig
      integer, allocatable :: glo(:,:), ghi(:,:)
      integer(int64), allocatable :: goff(:)
      character(len=MAX_NAME), allocatable :: gfile(:)
      integer :: sel_ncomp
      integer, allocatable :: sel(:)

      ierr = 0
      pf%path = path

      open(newunit=u, file=trim(path)//'/Header', status='old', &
           action='read', iostat=ios)
      if (ios /= 0) then
         ierr = 1; return                     ! no Header: not a plotfile
      end if

      read(u,'(a)',iostat=ios) line           ! version string
      if (ios /= 0) then; ierr = 2; return; end if
      read(u,*,iostat=ios) ncomp
      if (ios /= 0 .or. ncomp < 1) then; ierr = 3; return; end if

      allocate(pf%names(ncomp))
      do i = 1, ncomp
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) then; ierr = 4; return; end if
         pf%names(i) = trim(adjustl(line))
      end do

      read(u,*,iostat=ios) ndim
      if (ios /= 0) then; ierr = 5; return; end if
      if (ndim /= 2) then
         ierr = 100 + ndim                    ! only 2-D supported
         return
      end if
      read(u,*,iostat=ios) t
      if (ios /= 0) then; ierr = 6; return; end if
      read(u,*,iostat=ios) finest
      if (ios /= 0) then; ierr = 7; return; end if
      if (finest /= 0) then
         ierr = 200 + finest                  ! only single-level supported
         return
      end if

      read(u,*,iostat=ios) pf%prob_lo(1), pf%prob_lo(2)
      if (ios /= 0) then; ierr = 8; return; end if
      read(u,*,iostat=ios) pf%prob_hi(1), pf%prob_hi(2)
      if (ios /= 0) then; ierr = 9; return; end if

      read(u,'(a)',iostat=ios) line           ! ref ratios (blank at level 0)
      read(u,'(a)',iostat=ios) line           ! index space
      if (ios /= 0) then; ierr = 10; return; end if
      call parse_box(line, lo, hi, ierr)
      if (ierr /= 0) return
      nx = hi(1) - lo(1) + 1
      ny = hi(2) - lo(2) + 1

      read(u,'(a)',iostat=ios) line           ! level steps
      read(u,'(a)',iostat=ios) line           ! cell sizes
      read(u,'(a)',iostat=ios) line           ! coord_sys
      read(u,'(a)',iostat=ios) line           ! bwidth
      read(u,*,iostat=ios) i, ngrids          ! "0 ngrids time"
      if (ios /= 0 .or. ngrids < 1) then; ierr = 11; return; end if
      read(u,'(a)',iostat=ios) line           ! level steps
      pf%step = plotfile_step(path)
      do ig = 1, ngrids                        ! grid extents, 2 lines each
         read(u,'(a)',iostat=ios) line
         read(u,'(a)',iostat=ios) line
      end do
      read(u,'(a)',iostat=ios) line           ! e.g. "Level_0/Cell"
      if (ios /= 0) then; ierr = 12; return; end if
      leveldir = trim(path) // '/' // trim(adjustl(line))
      close(u)

      pf%time = t
      pf%nx = nx
      pf%ny = ny

      ! ---- which components to keep ------------------------------------
      if (present(want_names)) then
         sel_ncomp = size(want_names)
         allocate(sel(sel_ncomp))
         do i = 1, sel_ncomp
            sel(i) = -1
            do j = 1, ncomp
               if (trim(pf%names(j)) == trim(want_names(i))) then
                  sel(i) = j; exit
               end if
            end do
            if (sel(i) < 0) then
               ierr = 300 + i                 ! requested field not present
               return
            end if
         end do
      else
         sel_ncomp = ncomp
         allocate(sel(sel_ncomp))
         do i = 1, sel_ncomp
            sel(i) = i
         end do
      end if
      pf%ncomp = sel_ncomp

      ! ---- Cell_H: boxes and FAB offsets --------------------------------
      call read_cell_h(trim(leveldir)//'_H', glo, ghi, goff, gfile, &
                       ngrids, ncomp, ierr)
      if (ierr /= 0) return

      allocate(pf%data(nx, ny, sel_ncomp))
      pf%data = 0.0_dp

      do ig = 1, ngrids
         ! strip the trailing "/Cell" to get the containing directory
         call dirname(leveldir, fabfile)
         fabfile = trim(fabfile) // '/' // trim(gfile(ig))
         call read_fab_into(fabfile, goff(ig), glo(:,ig), ghi(:,ig), &
                            lo, ncomp, sel, pf%data, ierr)
         if (ierr /= 0) return
      end do
   end subroutine read_plotfile

   subroutine dirname(path, dir)
      character(len=*), intent(in) :: path
      character(len=*), intent(out) :: dir
      integer :: i
      dir = '.'
      do i = len_trim(path), 1, -1
         if (path(i:i) == '/') then
            dir = path(1:i-1)
            return
         end if
      end do
   end subroutine dirname

   !> Parse "((lo1,lo2) (hi1,hi2) (t1,t2))" into lo and hi.
   subroutine parse_box(line, lo, hi, ierr)
      character(len=*), intent(in) :: line
      integer, intent(out) :: lo(2), hi(2), ierr
      character(len=MAX_NAME) :: s
      integer :: i, n
      ierr = 0
      s = line
      ! replace all bracket and comma characters with blanks, then read the first four integers: lo1 lo2 hi1 hi2
      n = len_trim(s)
      do i = 1, n
         if (s(i:i) == '(' .or. s(i:i) == ')' .or. s(i:i) == ',') s(i:i) = ' '
      end do
      read(s, *, iostat=ierr) lo(1), lo(2), hi(1), hi(2)
      if (ierr /= 0) ierr = 400
   end subroutine parse_box

   !> Read the MultiFab header Cell_H.
   subroutine read_cell_h(fname, glo, ghi, goff, gfile, ngrids, ncomp, ierr)
      character(len=*), intent(in) :: fname
      integer, allocatable, intent(out) :: glo(:,:), ghi(:,:)
      integer(int64), allocatable, intent(out) :: goff(:)
      character(len=MAX_NAME), allocatable, intent(out) :: gfile(:)
      integer, intent(inout) :: ngrids
      integer, intent(in) :: ncomp
      integer, intent(out) :: ierr

      integer :: u, ios, i, nb, p
      character(len=MAX_NAME) :: line, tok
      integer :: lo(2), hi(2)
      integer :: dummy_ncomp

      ierr = 0
      open(newunit=u, file=trim(fname), status='old', action='read', &
           iostat=ios)
      if (ios /= 0) then; ierr = 500; return; end if

      read(u,'(a)',iostat=ios) line          ! version
      read(u,'(a)',iostat=ios) line          ! how
      read(u,*,iostat=ios) dummy_ncomp
      if (ios /= 0) then; ierr = 501; return; end if
      if (dummy_ncomp /= ncomp) then
         ierr = 502; return                  ! Header and Cell_H disagree
      end if
      read(u,'(a)',iostat=ios) line          ! nghost

      ! "(N 0"  -> N is the number of boxes
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) then; ierr = 503; return; end if
      tok = adjustl(line)
      if (tok(1:1) == '(') tok = tok(2:)
      read(tok, *, iostat=ios) nb
      if (ios /= 0 .or. nb < 1) then; ierr = 504; return; end if
      ngrids = nb

      allocate(glo(2,nb), ghi(2,nb), goff(nb), gfile(nb))
      do i = 1, nb
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) then; ierr = 505; return; end if
         call parse_box(line, lo, hi, ierr)
         if (ierr /= 0) return
         glo(:,i) = lo
         ghi(:,i) = hi
      end do

      read(u,'(a)',iostat=ios) line          ! ")"
      read(u,'(a)',iostat=ios) line          ! number of grids again
      do i = 1, nb
         read(u,'(a)',iostat=ios) line       ! "FabOnDisk: <file> <offset>"
         if (ios /= 0) then; ierr = 506; return; end if
         line = adjustl(line)
         p = index(line, ':')
         if (p > 0) line = adjustl(line(p+1:))
         read(line, *, iostat=ios) gfile(i), goff(i)
         if (ios /= 0) then; ierr = 507; return; end if
      end do
      close(u)
   end subroutine read_cell_h

   !> Read one FAB and scatter its components into `out`.
   subroutine read_fab_into(fname, offset, blo, bhi, dlo, ncomp, sel, out, ierr)
      character(len=*), intent(in) :: fname
      integer(int64), intent(in) :: offset
      integer, intent(in) :: blo(2), bhi(2), dlo(2), ncomp, sel(:)
      real(dp), intent(inout) :: out(:,:,:)
      integer, intent(out) :: ierr

      integer :: u, ios, i, j, c, k, bnx, bny, ncell
      integer(int64) :: pos
      logical :: little
      real(dp), allocatable :: buf(:)
      integer :: ii, jj

      ierr = 0
      bnx = bhi(1) - blo(1) + 1
      bny = bhi(2) - blo(2) + 1
      ncell = bnx * bny

      open(newunit=u, file=trim(fname), status='old', access='stream', &
           form='unformatted', action='read', iostat=ios)
      if (ios /= 0) then; ierr = 600; return; end if

      ! Skip the ASCII FAB header: read from `offset` to the first newline. The stored offset points at the START of the FAB, i.e. at the header.
      pos = offset + 1_int64          ! stream positions are 1-based
      call scan_fab_header(u, pos, little, ierr)
      if (ierr /= 0) then
         close(u); return
      end if

      allocate(buf(ncell))
      do c = 1, ncomp
         ! Only read components we actually want; skip the rest by advancing.
         k = 0
         do i = 1, size(sel)
            if (sel(i) == c) k = i
         end do
         if (k == 0) then
            pos = pos + int(ncell, int64) * 8_int64
            cycle
         end if
         read(u, pos=pos, iostat=ios) buf
         if (ios /= 0) then
            close(u); ierr = 601; return
         end if
         if (.not. little) call byteswap8(buf)
         do j = 1, bny
            jj = blo(2) - dlo(2) + j
            do i = 1, bnx
               ii = blo(1) - dlo(1) + i
               out(ii, jj, k) = buf(i + (j-1)*bnx)
            end do
         end do
         pos = pos + int(ncell, int64) * 8_int64
      end do
      close(u)
   end subroutine read_fab_into

   !> Advance `pos` past the FAB ASCII header and determine endianness.
   subroutine scan_fab_header(u, pos, little, ierr)
      integer, intent(in) :: u
      integer(int64), intent(inout) :: pos
      logical, intent(out) :: little
      integer, intent(out) :: ierr

      character(len=1) :: ch
      character(len=1024) :: hdr
      integer :: n, ios, p1, p2

      ierr = 0
      little = .true.
      hdr = ''
      n = 0
      do
         read(u, pos=pos, iostat=ios) ch
         if (ios /= 0) then; ierr = 610; return; end if
         pos = pos + 1_int64
         if (ch == char(10)) exit
         n = n + 1
         if (n <= len(hdr)) hdr(n:n) = ch
         if (n > 1000) then; ierr = 611; return; end if
      end do

      ! Endianness lives in the SECOND "(8, ( ... ))" block of the real descriptor: descending "8 7 6 5 4 3 2 1" is little endian, ascending "1 2 3 4 5 6 7 8" is big endian.
      p1 = index(hdr(1:n), '),(')
      if (p1 > 0) then
         p2 = index(hdr(p1:n), '(8, (')
         if (p2 > 0) then
            p2 = p1 + p2 + 4
            little = (hdr(p2:p2) == '8')
         end if
      end if
   end subroutine scan_fab_header

   !> Reverse the byte order of an array of doubles.
   subroutine byteswap8(a)
      real(dp), intent(inout) :: a(:)
      integer(int64) :: t
      integer :: i
      integer(int64) :: b(8), v
      do i = 1, size(a)
         t = transfer(a(i), t)
         v = 0_int64
         v = ior(ishft(iand(t, int(z'00000000000000FF', int64)), 56), v)
         v = ior(ishft(iand(t, int(z'000000000000FF00', int64)), 40), v)
         v = ior(ishft(iand(t, int(z'0000000000FF0000', int64)), 24), v)
         v = ior(ishft(iand(t, int(z'00000000FF000000', int64)),  8), v)
         v = ior(iand(ishft(t, -8),  int(z'00000000FF000000', int64)), v)
         v = ior(iand(ishft(t, -24), int(z'0000000000FF0000', int64)), v)
         v = ior(iand(ishft(t, -40), int(z'000000000000FF00', int64)), v)
         v = ior(iand(ishft(t, -56), int(z'00000000000000FF', int64)), v)
         a(i) = transfer(v, a(i))
         b = 0_int64
      end do
   end subroutine byteswap8

   !> Check that two plotfiles may legitimately be subtracted.
   !>
   !> Guards the silent-garbage failure mode where the perturbed and twin runs drifted onto different grids, times or field sets, so that the difference is meaningless rather than a linearised perturbation.
   function compatible(a, b, time_tol, msg) result(ok)
      type(plotfile_t), intent(in) :: a, b
      real(dp), intent(in) :: time_tol
      character(len=*), intent(out) :: msg
      logical :: ok
      real(dp) :: scale

      ok = .false.
      msg = ''
      if (a%ncomp /= b%ncomp) then
         msg = 'component count mismatch'; return
      end if
      if (a%nx /= b%nx .or. a%ny /= b%ny) then
         msg = 'grid mismatch: the two runs used different grids; '// &
               'set amr.max_level=0 and identical n_cell'
         return
      end if
      if (any(abs(a%prob_lo - b%prob_lo) > 0.0_dp) .or. &
          any(abs(a%prob_hi - b%prob_hi) > 0.0_dp)) then
         msg = 'domain extent mismatch between the two runs'; return
      end if
      scale = max(1.0_dp, abs(a%time), abs(b%time))
      if (abs(a%time - b%time) > time_tol * scale) then
         msg = 'time mismatch: both runs must use the same ns.fixed_dt '// &
               'and plot cadence'
         return
      end if
      ok = .true.
   end function compatible

end module plotfile_mod
