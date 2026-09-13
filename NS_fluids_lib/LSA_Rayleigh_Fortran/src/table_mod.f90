!> Formatted table output.
!>
!> Every report this program produces is a small table, and the reader of a thesis appendix or a batch log should be able to take in a growth rate, its residual and its departure from theory without counting columns. The routines here render aligned tables with a rule above and below the header, so that the alignment survives being pasted into a document or a log file.
!>
!> The renderer is deliberately plain ASCII rather than box-drawing characters. Output is read through terminals, log files, `less`, and text editors on machines whose locale settings are not known in advance, and a table that degrades into replacement characters is worse than one that was never drawn. Column widths are measured from the contents rather than fixed in advance, so a column of six-digit growth rates and a column of run labels of unequal length both come out aligned.
!>
!> Two shapes cover everything reported here. `table_t` is the column table, used where several cases are compared. `kv_block` is the aligned key-value block, used where a single case is summarised and there is one value per line; it is the same idea with the header suppressed and the first column left-aligned.
module table_mod
   implicit none
   private

   integer, parameter :: dp = kind(1.0d0)

   !> Maximum width of a single cell, and the maximum number of columns. Both are generous for the reports here and keep the implementation free of nested allocatable strings, which older compilers handle inconsistently.
   integer, parameter, public :: TCELL = 48
   integer, parameter, public :: TMAXCOL = 12

   public :: table_t, kv_block, rule, fmt_f, fmt_e, fmt_i, fmt_pct, cel

   type :: table_t
      private
      integer :: ncol = 0
      integer :: nrow = 0
      character(len=TCELL) :: head(TMAXCOL) = ''
      !> Alignment of each column: 'l' or 'r'. Numbers are right-aligned so that the decimal points line up; labels are left-aligned.
      character(len=1) :: align(TMAXCOL) = 'r'
      character(len=TCELL), allocatable :: cell(:,:)
   contains

      procedure :: init  => table_init
      procedure :: row   => table_row
      procedure :: emit  => table_emit
   end type table_t

contains

   !> Coerce a string to the fixed cell width. Fortran requires every element of an array constructor to have the same length, so cells and headers assembled inline are passed through this rather than padded by hand, which would otherwise have to be redone whenever a label changed.
   pure function cel(s) result(out)
      character(len=*), intent(in) :: s
      character(len=TCELL) :: out
      out = s
   end function cel

   !> Begin a table. `headers` names the columns and fixes their number; `aligns` is an optional string of 'l' and 'r', one character per column, defaulting to left for the first column and right for the rest, which is the layout every report here wants.
   subroutine table_init(this, headers, aligns)
      class(table_t), intent(inout) :: this
      character(len=*), intent(in) :: headers(:)
      character(len=*), intent(in), optional :: aligns
      integer :: i

      this%ncol = min(size(headers), TMAXCOL)
      this%nrow = 0
      do i = 1, this%ncol
         this%head(i) = headers(i)
         if (present(aligns)) then
            if (i <= len(aligns)) then
               this%align(i) = aligns(i:i)
            else
               this%align(i) = 'r'
            end if
         else
            if (i == 1) then
               this%align(i) = 'l'
            else
               this%align(i) = 'r'
            end if
         end if
      end do
      if (allocated(this%cell)) deallocate(this%cell)
      allocate(this%cell(64, this%ncol))
      this%cell = ''
   end subroutine table_init

   !> Append one row. Cells are supplied already formatted, which keeps the numeric formatting decisions with the code that knows what the number means rather than in the renderer.
   subroutine table_row(this, cells)
      class(table_t), intent(inout) :: this
      character(len=*), intent(in) :: cells(:)
      character(len=TCELL), allocatable :: bigger(:,:)
      integer :: i, n

      if (this%nrow >= size(this%cell, 1)) then
         n = 2*size(this%cell, 1)
         allocate(bigger(n, this%ncol))
         bigger = ''
         bigger(1:this%nrow, :) = this%cell(1:this%nrow, :)
         call move_alloc(bigger, this%cell)
      end if

      this%nrow = this%nrow + 1
      do i = 1, min(size(cells), this%ncol)
         this%cell(this%nrow, i) = cells(i)
      end do
   end subroutine table_row

   !> Render the table. `title`, when present and non-empty, is written above it.
   subroutine table_emit(this, title)
      class(table_t), intent(in) :: this
      character(len=*), intent(in), optional :: title
      integer :: w(TMAXCOL), i, j, total
      character(len=:), allocatable :: line

      if (this%ncol <= 0) return

      ! Width of each column is the widest of its header and its cells.
      do j = 1, this%ncol
         w(j) = len_trim(this%head(j))
         do i = 1, this%nrow
            w(j) = max(w(j), len_trim(this%cell(i, j)))
         end do
      end do

      ! Two spaces of gutter between columns.
      total = 0
      do j = 1, this%ncol
         total = total + w(j)
      end do
      total = total + 2*(this%ncol - 1)

      if (present(title)) then
         if (len_trim(title) > 0) then
            write(*,'(a)') ''
            write(*,'(a)') trim(title)
         end if
      end if

      line = ''
      do j = 1, this%ncol
         line = line//pad(this%head(j), w(j), this%align(j))
         if (j < this%ncol) line = line//'  '
      end do
      write(*,'(a)') line
      write(*,'(a)') repeat('-', total)

      do i = 1, this%nrow
         line = ''
         do j = 1, this%ncol
            line = line//pad(this%cell(i, j), w(j), this%align(j))
            if (j < this%ncol) line = line//'  '
         end do
         write(*,'(a)') line
      end do
      write(*,'(a)') repeat('-', total)
   end subroutine table_emit

   !> An aligned key-value block, for summarising a single case. Labels are padded to a common width so that the values form a column, and `unit`, when supplied, follows the value.
   subroutine kv_block(title, labels, values, width)
      character(len=*), intent(in) :: title
      character(len=*), intent(in) :: labels(:)
      character(len=*), intent(in) :: values(:)
      integer, intent(in), optional :: width
      integer :: i, wl, wv, n

      n = min(size(labels), size(values))
      if (n <= 0) return

      wl = 0
      wv = 0
      do i = 1, n
         wl = max(wl, len_trim(labels(i)))
         wv = max(wv, len_trim(values(i)))
      end do
      if (present(width)) wl = max(wl, width)

      if (len_trim(title) > 0) then
         write(*,'(a)') ''
         write(*,'(a)') trim(title)
         write(*,'(a)') repeat('-', wl + 2 + wv)
      end if

      do i = 1, n
         write(*,'(a)') pad(labels(i), wl, 'l')//'  '//pad(values(i), wv, 'r')
      end do
   end subroutine kv_block

   !> A horizontal rule of the given width, used to close a block opened by kv_block.
   subroutine rule(n)
      integer, intent(in) :: n
      write(*,'(a)') repeat('-', max(n, 1))
   end subroutine rule

   !> Pad a string to a width, left- or right-aligned. Content wider than the field is returned untruncated, since a mangled number is worse than a misaligned one.
   pure function pad(s, w, how) result(out)
      character(len=*), intent(in) :: s
      integer, intent(in) :: w
      character(len=1), intent(in) :: how
      character(len=:), allocatable :: out
      integer :: n

      n = len_trim(s)
      if (n >= w) then
         out = trim(s)
      else if (how == 'l') then
         out = trim(s)//repeat(' ', w - n)
      else
         out = repeat(' ', w - n)//trim(s)
      end if
   end function pad

   !> Fixed-point number as a trimmed string, with `d` decimal places.
   pure function fmt_f(v, d) result(out)
      real(dp), intent(in) :: v
      integer, intent(in) :: d
      character(len=:), allocatable :: out
      character(len=TCELL) :: buf
      character(len=16) :: f

      write(f,'(a,i0,a,i0,a)') '(f', TCELL, '.', d, ')'
      write(buf, f) v
      out = trim(adjustl(buf))
   end function fmt_f

   !> Number in exponent form, for amplitudes and residuals spanning many decades.
   pure function fmt_e(v) result(out)
      real(dp), intent(in) :: v
      character(len=:), allocatable :: out
      character(len=TCELL) :: buf

      write(buf,'(es12.3)') v
      out = trim(adjustl(buf))
   end function fmt_e

   !> Integer as a trimmed string.
   pure function fmt_i(n) result(out)
      integer, intent(in) :: n
      character(len=:), allocatable :: out
      character(len=TCELL) :: buf

      write(buf,'(i0)') n
      out = trim(adjustl(buf))
   end function fmt_i

   !> A percentage, given the fraction. Written with a trailing sign character so that a column of them reads as percentages without a unit in the header.
   pure function fmt_pct(frac, d) result(out)
      real(dp), intent(in) :: frac
      integer, intent(in), optional :: d
      character(len=:), allocatable :: out
      character(len=TCELL) :: buf
      character(len=16) :: f
      integer :: nd

      nd = 2
      if (present(d)) nd = d
      write(f,'(a,i0,a,i0,a)') '(f', TCELL, '.', nd, ')'
      write(buf, f) 100.0_dp*frac
      out = trim(adjustl(buf))//' %'
   end function fmt_pct

end module table_mod
