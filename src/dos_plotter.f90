module dos_plotter
    !! DOS visualization: interactive ASCII terminal plot + CSV export
    !! Interactive controls: ↑↓ scroll, +/- zoom, ← → pan, R reset, Q quit
    use castep_config, only: dp
    use term_utils, only: C_RED, C_GREEN, C_YELLOW, C_CYAN, C_BOLD, C_DIM, C_RESET, C_AXIS, &
        get_term_size, draw_line
    implicit none
    private

    integer, parameter, public :: &
        DOS_MODE_ASCII  = 1, &
        DOS_MODE_EXPORT = 3

    public :: plot_dos_ascii, write_dos_csv
    public :: plot_pdos_ascii, write_pdos_csv

contains

    subroutine plot_dos_ascii(energy_grid, dos_data, nspin, e_fermi, smearing, &
            term_w_in, term_h_in, y_center_in, y_half_in, e_center_in, half_range_in, &
            xlabel, xunit, hide_fermi, title)
        real(dp), intent(in) :: energy_grid(:)
        real(dp), intent(in) :: dos_data(:,:)
        integer,  intent(in) :: nspin
        real(dp), intent(in) :: e_fermi, smearing
        integer,  intent(in) :: term_w_in, term_h_in
        real(dp), intent(in), optional :: y_center_in, y_half_in, e_center_in, half_range_in
        character(len=*), intent(in), optional :: xlabel, xunit
        logical, intent(in), optional :: hide_fermi
        character(len=*), intent(in), optional :: title

        integer  :: ne, nw, nh, ix, iy, ie, is, nw_data
        integer  :: fermi_col, y_label_width, gap_width, y_label_interval
        integer  :: i1, i2, ne_vis
        real(dp) :: e_min, e_max, e_range, y_center, y_half, y_min, y_max
        real(dp) :: y_scale, x_scale, x_val, dos_val, e_center, half_range
        character(len=1), allocatable :: grid(:,:)
        character(len=64)  :: fmt_label, tmp_str
        character(len=16)  :: xl, xu
        integer  :: last_ix, last_iy
        logical  :: no_fermi

        no_fermi = .false.
        if (present(hide_fermi)) no_fermi = hide_fermi

        xl = 'Energy'; if (present(xlabel)) xl = trim(xlabel)
        xu = 'eV';     if (present(xunit))  xu = trim(xunit)

        ne = size(energy_grid)
        if (ne < 2 .or. size(dos_data, 1) /= ne) then
            write(*, '(a)') '  No DOS data to plot.'
            return
        end if

        ! default x view: full range
        e_center = 0.0_dp
        half_range = energy_grid(ne) - energy_grid(1)
        if (present(e_center_in))    e_center   = e_center_in
        if (present(half_range_in))  half_range = half_range_in

        if (half_range < 0.25_dp) half_range = 0.25_dp
        e_min = e_center - half_range
        e_max = e_center + half_range
        e_range = e_max - e_min
        if (e_range < 1.0e-12_dp) e_range = 1.0_dp

        call find_visible_range(energy_grid, e_min, e_max, i1, i2)
        ne_vis = i2 - i1 + 1

        call get_term_size(nw, nh)
        if (term_w_in > 0) nw = term_w_in
        if (term_h_in > 0) nh = term_h_in
        nw = max(30, min(300, nw))
        nh = max(15, min(100, nh))
        nh = max(15, nh - 9)

        ! y-axis: center + half-range; default from full data: [0, y_max0]
        y_max = 1.0_dp
        do is = 1, nspin
            do ie = 1, ne
                if (dos_data(ie, is) > y_max) y_max = dos_data(ie, is)
            end do
        end do
        y_max = y_max * 1.15_dp
        if (y_max < 1.0e-12_dp) y_max = 1.0_dp
        y_center = y_max / 2.0_dp
        y_half   = y_max / 2.0_dp
        if (present(y_center_in)) y_center = y_center_in
        if (present(y_half_in))   y_half   = y_half_in
        if (y_half < 1.0e-12_dp) y_half = 0.5_dp
        y_min = y_center - y_half
        y_max = y_center + y_half

        y_label_width = 0
        write(tmp_str, '(f12.4)') y_max
        y_label_width = max(y_label_width, len_trim(adjustl(tmp_str)))
        write(tmp_str, '(f12.4)') y_min
        y_label_width = max(y_label_width, len_trim(adjustl(tmp_str)))
        y_label_width = max(y_label_width, 6)
        gap_width = 1
        nw_data = max(20, nw - y_label_width - gap_width)

        allocate(grid(nw_data, nh))
        grid = ' '

        x_scale = real(nw_data - 3, dp) / e_range
        y_scale = real(nh - 1, dp) / (y_max - y_min)

        fermi_col = nint((0.0_dp - e_min) * x_scale) + 2
        fermi_col = max(2, min(nw_data - 1, fermi_col))

        ! plot DOS curves
        do is = 1, nspin
            last_ix = -1
            do ie = i1, i2
                x_val = energy_grid(ie)
                if (x_val < e_min .or. x_val > e_max) cycle
                dos_val = dos_data(ie, is)
                ix = nint((x_val - e_min) * x_scale) + 2
                iy = nh - nint((dos_val - y_min) * y_scale)
                ix = max(2, min(nw_data - 1, ix))
                iy = max(1, min(nh, iy))
                if (is == 1) then
                    grid(ix, iy) = 'U'
                else
                    grid(ix, iy) = 'D'
                end if
                if (last_ix > 0) call draw_line(nw_data, nh, grid, last_ix, last_iy, ix, iy)
                last_ix = ix; last_iy = iy
            end do
        end do

        ! Fermi level vertical line
        if (.not. no_fermi) then
            do iy = 1, nh
                if (fermi_col >= 2 .and. fermi_col <= nw_data - 1) then
                    if (grid(fermi_col, iy) == ' ') grid(fermi_col, iy) = '|'
                end if
            end do
        end if

        ! y=0 horizontal reference line
        iy = nh - nint((0.0_dp - y_min) * y_scale)
        iy = max(1, min(nh, iy))
        do ix = 2, nw_data - 1
            if (grid(ix, iy) == ' ') grid(ix, iy) = '-'
        end do
        if (grid(1, iy) == '|') grid(1, iy) = 'Y'
        if (grid(nw_data, iy) == '|') grid(nw_data, iy) = 'Z'

        ! borders
        do iy = 1, nh
            if (grid(1, iy) == ' ') grid(1, iy) = '|'
            if (grid(nw_data, iy) == ' ') grid(nw_data, iy) = '|'
        end do
        do ix = 2, nw_data - 1
            if (grid(ix, nh) == ' ' .or. grid(ix, nh) == '-') grid(ix, nh) = '-'
        end do
        grid(1, nh) = 'B'
        grid(nw_data, nh) = 'R'
        if (.not. no_fermi) then
            if (grid(fermi_col, nh) == '-') grid(fermi_col, nh) = 'T'
        end if

        ! header
        write(*, '(a)') ''
        if (present(title)) then
            write(*, '(a)') C_BOLD // '  ' // trim(title) // '  ' // C_RESET
        else
            write(*, '(a)') C_BOLD // '  Density of States  ' // C_RESET
        end if
        write(*, '(a,i0,a,i0,a,f5.2,a)') &
            C_RESET // C_CYAN, ne_vis, C_RESET // ' pts, ' // C_CYAN, &
            nspin, C_RESET // ' spin(s),  smearing=' // C_CYAN, smearing, &
            ' ' // trim(xu) // C_RESET
        write(tmp_str, '(f8.4)') e_min
        write(*, '(a)', advance='no') '  ' // C_DIM // trim(xl) // ': [' // C_RESET // &
            trim(adjustl(tmp_str)) // C_DIM // ' to ' // C_RESET
        write(tmp_str, '(f8.4)') e_max
        if (no_fermi) then
            write(*, '(a,i0,1x,i0)') trim(adjustl(tmp_str)) // C_DIM // '] ' // trim(xu) // &
                C_RESET // '  grid:' // C_RESET, nw_data, nh
        else
            write(*, '(a,i0,1x,i0,a)') trim(adjustl(tmp_str)) // C_DIM // '] ' // trim(xu) // &
                C_RESET // '  grid:' // C_RESET, nw_data, nh, &
                '  ' // C_DIM // 'E_F=0 (dashed)' // C_RESET
        end if
        write(tmp_str, '(f8.4)') y_max
        write(*, '(a)') '  ' // C_DIM // 'Y max: ' // C_RESET // trim(adjustl(tmp_str))
        write(*, '(a)') ''

        ! top border
        write(*, '(a)') repeat(' ', y_label_width + gap_width) // C_AXIS &
            // '┌' // repeat('─', nw_data - 2) // '┐' // C_RESET

        y_label_interval = max(1, nh / 6)
        do iy = 1, nh
            fmt_label = ''
            if (mod(iy - 1, y_label_interval) == 0 .or. iy == 1 .or. iy == nh) then
                dos_val = y_max - real(iy - 1, dp) / real(nh - 1, dp) * (y_max - y_min)
                write(fmt_label, '(f12.4)') dos_val
                fmt_label = adjustl(fmt_label)
            end if
            call print_dos_row(grid, nw_data, iy, nh, fmt_label, y_label_width, gap_width)
        end do

        ! x-axis
        if (present(xlabel)) then
            call write_x_axis(y_label_width, gap_width, nw_data, e_min, e_max)
        else
            write(tmp_str, '(f8.2)') e_max
            tmp_str = adjustl(tmp_str)
            write(fmt_label, '(f8.2)') e_min
            fmt_label = adjustl(fmt_label)
            write(*, '(a,a,a)') repeat(' ', y_label_width + gap_width) // C_DIM, &
                trim(fmt_label) // repeat(' ', max(0, nw_data - 3 - len_trim(tmp_str))) &
                // trim(tmp_str), C_RESET
        end if

        ! E_F marker
        if (.not. no_fermi) then
            write(*, '(a,a)') repeat(' ', y_label_width + gap_width) // C_DIM, &
                repeat(' ', fermi_col - 2) // C_RED // 'E_F' // C_RESET
        end if

        deallocate(grid)
    end subroutine plot_dos_ascii

    subroutine write_x_axis(label_w, gap, nw, v_min, v_max)
        integer, intent(in) :: label_w, gap, nw
        real(dp), intent(in) :: v_min, v_max
        integer :: n_ticks, i, pos, j, lbl_len
        real(dp) :: tick_val
        character(len=8) :: lbl
        character(len=512) :: line

        ! Multi-tick mode: one tick per ~16 chars of plot width
        n_ticks = max(3, nw / 16)
        line = repeat(' ', nw)

        do i = 0, n_ticks - 1
            tick_val = v_min + i * (v_max - v_min) / real(n_ticks - 1, dp)
            write(lbl, '(f8.1)') tick_val
            lbl = adjustl(lbl)
            lbl_len = len_trim(lbl)
            pos = 1 + i * (nw - 1) / max(1, n_ticks - 1)
            pos = pos - lbl_len / 2
            pos = max(1, min(nw - lbl_len + 1, pos))
            do j = 1, lbl_len
                if (pos + j - 1 <= nw) line(pos + j - 1:pos + j - 1) = lbl(j:j)
            end do
        end do

        write(*, '(a,a,a)') repeat(' ', label_w + gap) // C_DIM, &
            trim(line), C_RESET
    end subroutine write_x_axis

    ! find indices of energy_grid within [e_min, e_max]
    subroutine find_visible_range(grid, v_min, v_max, i1, i2)
        real(dp), intent(in)  :: grid(:), v_min, v_max
        integer,  intent(out) :: i1, i2
        integer :: n, lo, hi, mid
        n = size(grid)
        ! binary search for first index >= v_min
        lo = 1; hi = n
        do while (lo < hi)
            mid = (lo + hi) / 2
            if (grid(mid) < v_min) then; lo = mid + 1
            else;                        hi = mid
            end if
        end do
        i1 = max(1, lo - 1)
        ! binary search for last index <= v_max
        lo = 1; hi = n
        do while (lo < hi)
            mid = (lo + hi + 1) / 2
            if (grid(mid) > v_max) then; hi = mid - 1
            else;                        lo = mid
            end if
        end do
        i2 = min(n, lo + 1)
        i1 = max(1, min(n, i1))
        i2 = max(1, min(n, i2))
    end subroutine find_visible_range

    subroutine print_dos_row(grid, nw, iy, nh, label, lbl_w, gap_w)
        integer, intent(in) :: nw, iy, nh, lbl_w, gap_w
        character(len=1), intent(in) :: grid(nw, nh)
        character(len=*), intent(in) :: label
        character(len=4096) :: buf
        character(len=1)   :: ch
        character(len=3)   :: dc
        integer :: ix, bp, labellen, n, slen, run_start, run_len, rt
        character(len=16) :: cur_color, clr
        logical :: is_fermi_line

        bp = 1
        labellen = len_trim(label)
        if (labellen > 0) then
            n = max(0, lbl_w - labellen)
            if (n > 0) then; buf(bp:bp+n-1) = repeat(' ', n); bp = bp + n; end if
            slen = len_trim(C_DIM); buf(bp:bp+slen-1) = trim(C_DIM); bp = bp + slen
            buf(bp:bp+labellen-1) = label(1:labellen); bp = bp + labellen
            slen = len_trim(C_RESET); buf(bp:bp+slen-1) = C_RESET; bp = bp + slen
            buf(bp:bp+gap_w-1) = repeat(' ', gap_w); bp = bp + gap_w
        else
            n = lbl_w + gap_w
            buf(bp:bp+n-1) = repeat(' ', n); bp = bp + n
        end if

        if (nw >= 1) then
            run_start = 1
            do while (run_start <= nw)
                ch = grid(run_start, iy)
                is_fermi_line = (ch == '|' .and. run_start > 1 .and. run_start < nw)
                rt = dos_char_type(ch, iy == nh)
                run_len = 1
                do ix = run_start + 1, nw
                    if (dos_char_type(grid(ix, iy), iy == nh) == rt) then
                        run_len = run_len + 1
                    else
                        exit
                    end if
                end do
                cur_color = dos_type_color(rt, is_fermi_line)
                if (len_trim(cur_color) > 0) then
                    clr = trim(cur_color)
                    slen = len_trim(clr)
                    buf(bp:bp+slen-1) = clr; bp = bp + slen
                end if
                do ix = 1, run_len
                    dc = dos_type_char(grid(run_start + ix - 1, iy), rt)
                    slen = len_trim(dc)
                    if (slen == 0) slen = 1
                    buf(bp:bp+slen-1) = dc(1:slen); bp = bp + slen
                end do
                if (len_trim(cur_color) > 0) then
                    slen = len_trim(C_RESET)
                    buf(bp:bp+slen-1) = C_RESET; bp = bp + slen
                end if
                run_start = run_start + run_len
            end do
        end if

        write(*, '(a)') buf(1:bp-1)
    end subroutine print_dos_row

    pure function dos_char_type(ch, is_bottom) result(ct)
        character(len=1), intent(in) :: ch
        logical, intent(in) :: is_bottom
        integer :: ct
        ct = 0
        if (ch == '|') then; ct = 1
        else if (ch == 'B' .or. ch == 'R' .or. ch == 'T' .or. ch == 'Y' .or. ch == 'Z') then; ct = 1
        else if (ch == '-') then; ct = 2
        else if (ch == 'U') then; ct = 3
        else if (ch == 'D') then; ct = 4
        else if (ch == '.') then; ct = 5
        else if (ch == 'S') then; ct = 6    ! s-orbital PDOS
        else if (ch == 'P') then; ct = 7    ! p-orbital PDOS
        else if (ch == 'L') then; ct = 8    ! d-orbital PDOS ('D' already used)
        else if (ch == 'F') then; ct = 9    ! f-orbital PDOS
        end if
    end function dos_char_type

    pure function dos_type_color(ct, is_fermi) result(clr)
        integer, intent(in) :: ct
        logical, intent(in) :: is_fermi
        character(len=16) :: clr
        select case (ct)
        case (1)
            if (is_fermi) then; clr = C_RED
            else;               clr = C_AXIS
            end if
        case (2);  clr = C_AXIS
        case (3);  clr = C_CYAN
        case (4);  clr = C_YELLOW
        case (5);  clr = C_DIM
        case (6);  clr = C_CYAN      ! s
        case (7);  clr = C_YELLOW    ! p
        case (8);  clr = C_GREEN     ! d
        case (9);  clr = C_RED       ! f
        case default
            if (is_fermi) then; clr = C_RED
            else;               clr = ''
            end if
        end select
    end function dos_type_color

    pure function dos_type_char(ch, ct) result(dc)
        character(len=1), intent(in) :: ch
        integer, intent(in) :: ct
        character(len=3) :: dc
        select case (ct)
        case (1)
            select case (ch)
            case ('B'); dc = '└'
            case ('R'); dc = '┘'
            case ('T'); dc = '┴'
            case ('Y'); dc = '├'    ! y=0 left junction
            case ('Z'); dc = '┤'    ! y=0 right junction
            case default; dc = '│'
            end select
        case (2);  dc = '─'
        case (3);  dc = '●'
        case (4);  dc = '○'
        case (5);  dc = '·'
        case (6);  dc = '●'     ! s-orbital
        case (7);  dc = '○'     ! p-orbital
        case (8);  dc = '△'     ! d-orbital
        case (9);  dc = '▽'     ! f-orbital
        case default; dc = ' '
        end select
    end function dos_type_char



    ! ----------------------------------------------------------------
    !  CSV export for Origin / other plotting software
    ! ----------------------------------------------------------------
    subroutine write_dos_csv(energy_grid, dos_data, nspin, csv_file, iostat, iomsg)
        real(dp), intent(in) :: energy_grid(:)
        real(dp), intent(in) :: dos_data(:,:)
        integer,  intent(in) :: nspin
        character(len=*), intent(in) :: csv_file
        integer,  intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer :: unit, ios, ie, ne
        character(len=128) :: line

        iostat = 0
        if (present(iomsg)) iomsg = ''
        ne = size(energy_grid)
        if (ne < 1) then
            iostat = 1; if (present(iomsg)) iomsg = 'No data to export'; return
        end if

        open(newunit=unit, file=trim(csv_file), status='replace', &
             action='write', iostat=ios)
        if (ios /= 0) then
            iostat = ios
            if (present(iomsg)) iomsg = 'Cannot write: ' // trim(csv_file)
            return
        end if

        ! header
        if (nspin == 1) then
            write(unit, '(a)') '# Energy(eV),DOS'
        else
            write(unit, '(a)') '# Energy(eV),DOS_up,DOS_down'
        end if

        do ie = 1, ne
            if (nspin == 1) then
                write(line, '(f12.6,a,es14.6)') energy_grid(ie), ',', dos_data(ie, 1)
            else
                write(line, '(f12.6,a,es14.6,a,es14.6)') &
                    energy_grid(ie), ',', dos_data(ie, 1), ',', dos_data(ie, 2)
            end if
            ! strip internal spaces from Fortran formatted output
            line = adjustl(line)
            write(unit, '(a)') trim(line)
        end do

        close(unit)
    end subroutine write_dos_csv
    ! ----------------------------------------------------------------
    !  PDOS multi-channel ASCII plot (s, p, d, f)
    ! ----------------------------------------------------------------
    subroutine plot_pdos_ascii(energy_grid, pdos_data, e_fermi, smearing, &
            term_w_in, term_h_in, y_center_in, y_half_in, e_center_in, half_range_in)
        real(dp), intent(in) :: energy_grid(:)
        real(dp), intent(in) :: pdos_data(:,:)
        real(dp), intent(in) :: e_fermi, smearing
        integer,  intent(in) :: term_w_in, term_h_in
        real(dp), intent(in), optional :: y_center_in, y_half_in, e_center_in, half_range_in

        integer, parameter :: MAX_CH = 5
        character(len=1), parameter :: ch_mark(5) = ['U', 'S', 'P', 'L', 'F']  ! tot,s,p,d,f (tot uses 'U'=type 3=●)
        integer  :: ne, nw, nh, ix, iy, ie, ich, nw_data
        integer  :: fermi_col, y_label_width, gap_width, y_label_interval
        integer  :: i1, i2
        real(dp) :: e_min, e_max, e_range, y_center, y_half, y_min, y_max
        real(dp) :: y_scale, x_scale, x_val, dos_val
        real(dp) :: e_center, half_range
        character(len=1), allocatable :: grid(:,:)
        character(len=64)  :: fmt_label, tmp_str
        integer  :: last_ix, last_iy

        ne = size(energy_grid)
        if (ne < 2 .or. size(pdos_data, 1) /= ne) then
            write(*, '(a)') '  No PDOS data to plot.'; return
        end if

        e_center = 0.0_dp
        half_range = energy_grid(ne) - energy_grid(1)
        if (present(e_center_in))   e_center   = e_center_in
        if (present(half_range_in)) half_range = half_range_in

        if (half_range < 0.25_dp) half_range = 0.25_dp
        e_min = e_center - half_range
        e_max = e_center + half_range
        e_range = e_max - e_min
        if (e_range < 1.0e-12_dp) e_range = 1.0_dp

        call find_visible_range(energy_grid, e_min, e_max, i1, i2)

        call get_term_size(nw, nh)
        if (term_w_in > 0) nw = term_w_in
        if (term_h_in > 0) nh = term_h_in
        nw = max(30, min(300, nw))
        nh = max(15, min(100, nh))
        nh = max(15, nh - 10)

        y_max = 1.0_dp
        do ich = 2, MAX_CH
            do ie = 1, ne
                if (pdos_data(ie, ich) > y_max) y_max = pdos_data(ie, ich)
            end do
        end do
        y_max = y_max * 1.15_dp
        if (y_max < 1.0e-12_dp) y_max = 1.0_dp
        y_center = y_max / 2.0_dp
        y_half   = y_max / 2.0_dp
        if (present(y_center_in)) y_center = y_center_in
        if (present(y_half_in))   y_half   = y_half_in
        if (y_half < 1.0e-12_dp) y_half = 0.5_dp
        y_min = y_center - y_half
        y_max = y_center + y_half

        y_label_width = 0
        write(tmp_str, '(f12.4)') y_max
        y_label_width = max(y_label_width, len_trim(adjustl(tmp_str)))
        write(tmp_str, '(f12.4)') y_min
        y_label_width = max(y_label_width, len_trim(adjustl(tmp_str)))
        y_label_width = max(y_label_width, 6)
        gap_width = 1
        nw_data = max(20, nw - y_label_width - gap_width)

        allocate(grid(nw_data, nh))
        grid = ' '

        x_scale = real(nw_data - 3, dp) / e_range
        y_scale = real(nh - 1, dp) / (y_max - y_min)

        fermi_col = nint((0.0_dp - e_min) * x_scale) + 2
        fermi_col = max(2, min(nw_data - 1, fermi_col))

        ! Plot s, p, d, f channels (skip total = index 1)
        do ich = 2, MAX_CH
            last_ix = -1
            do ie = i1, i2
                x_val = energy_grid(ie)
                if (x_val < e_min .or. x_val > e_max) cycle
                dos_val = pdos_data(ie, ich)
                ix = nint((x_val - e_min) * x_scale) + 2
                iy = nh - nint((dos_val - y_min) * y_scale)
                ix = max(2, min(nw_data - 1, ix))
                iy = max(1, min(nh, iy))
                grid(ix, iy) = ch_mark(ich)   ! U=total, S=s, P=p, L=d, F=f
                if (last_ix > 0) call draw_line(nw_data, nh, grid, last_ix, last_iy, ix, iy)
                last_ix = ix; last_iy = iy
            end do
        end do

        ! Fermi level
        do iy = 1, nh
            if (fermi_col >= 2 .and. fermi_col <= nw_data - 1) then
                if (grid(fermi_col, iy) == ' ') grid(fermi_col, iy) = '|'
            end if
        end do

        ! y=0 horizontal reference line
        iy = nh - nint((0.0_dp - y_min) * y_scale)
        iy = max(1, min(nh, iy))
        do ix = 2, nw_data - 1
            if (grid(ix, iy) == ' ') grid(ix, iy) = '-'
        end do
        if (grid(1, iy) == '|') grid(1, iy) = 'Y'
        if (grid(nw_data, iy) == '|') grid(nw_data, iy) = 'Z'

        ! borders
        do iy = 1, nh
            if (grid(1, iy) == ' ') grid(1, iy) = '|'
            if (grid(nw_data, iy) == ' ') grid(nw_data, iy) = '|'
        end do
        do ix = 2, nw_data - 1
            if (grid(ix, nh) == ' ' .or. grid(ix, nh) == '-') grid(ix, nh) = '-'
        end do
        grid(1, nh) = 'B'; grid(nw_data, nh) = 'R'
        if (grid(fermi_col, nh) == '-') grid(fermi_col, nh) = 'T'

        ! header
        write(*, '(a)') ''
        write(*, '(a,i0,a,f5.2,a)') C_BOLD // '  Projected DOS  ' // &
            C_RESET // C_CYAN, i2 - i1 + 1, C_RESET // ' pts,  smearing=' // &
            C_CYAN, smearing, ' eV' // C_RESET
        write(tmp_str, '(f8.4)') e_min
        write(*, '(a)', advance='no') '  ' // C_DIM // 'Energy: [' // C_RESET // &
            trim(adjustl(tmp_str)) // C_DIM // ' to ' // C_RESET
        write(tmp_str, '(f8.4)') e_max
        write(*, '(a,i0,1x,i0,a)') trim(adjustl(tmp_str)) // C_DIM // '] eV' // &
            C_RESET // '  grid:' // C_RESET, nw_data, nh
        write(tmp_str, '(f8.4)') y_max
        write(*, '(a)') '  ' // C_DIM // 'DOS max: ' // C_RESET // trim(adjustl(tmp_str))
        write(*, '(a)') ''

        ! top border
        write(*, '(a)') repeat(' ', y_label_width + gap_width) // C_AXIS &
            // '┌' // repeat('─', nw_data - 2) // '┐' // C_RESET

        y_label_interval = max(1, nh / 6)
        do iy = 1, nh
            fmt_label = ''
            if (mod(iy - 1, y_label_interval) == 0 .or. iy == 1 .or. iy == nh) then
                dos_val = y_max - real(iy - 1, dp) / real(nh - 1, dp) * (y_max - y_min)
                write(fmt_label, '(f12.4)') dos_val
                fmt_label = adjustl(fmt_label)
            end if
            call print_dos_row(grid, nw_data, iy, nh, fmt_label, y_label_width, gap_width)
        end do

        ! x-axis
        write(tmp_str, '(f8.2)') e_max
        tmp_str = adjustl(tmp_str)
        write(fmt_label, '(f8.2)') e_min
        fmt_label = adjustl(fmt_label)
        write(*, '(a,a,a)') repeat(' ', y_label_width + gap_width) // C_DIM, &
            trim(fmt_label) // repeat(' ', max(0, nw_data - 3 - len_trim(tmp_str))) &
            // trim(tmp_str), C_RESET

        write(*, '(a,a)') repeat(' ', y_label_width + gap_width) // C_DIM, &
            repeat(' ', fermi_col - 2) // C_RED // 'E_F' // C_RESET

        ! legend
        write(*, '(a)') ''
        write(*, '(a)') '  ' // C_CYAN   // '● s'   // C_RESET // '  ' &
                          // C_YELLOW // '○ p'   // C_RESET // '  ' &
                          // C_GREEN  // '△ d'   // C_RESET // '  ' &
                          // C_RED    // '▽ f'   // C_RESET

        deallocate(grid)
    end subroutine plot_pdos_ascii

    ! ----------------------------------------------------------------
    !  PDOS CSV export (multi-channel)
    ! ----------------------------------------------------------------
    subroutine write_pdos_csv(energy_grid, pdos_data, csv_file, iostat, iomsg)
        real(dp), intent(in) :: energy_grid(:)
        real(dp), intent(in) :: pdos_data(:,:)
        character(len=*), intent(in) :: csv_file
        integer,  intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer :: unit, ios, ie, ne

        iostat = 0
        if (present(iomsg)) iomsg = ''
        ne = size(energy_grid)
        if (ne < 1) then
            iostat = 1; if (present(iomsg)) iomsg = 'No data'; return
        end if

        open(newunit=unit, file=trim(csv_file), status='replace', &
             action='write', iostat=ios)
        if (ios /= 0) then
            iostat = ios
            if (present(iomsg)) iomsg = 'Cannot write: ' // trim(csv_file)
            return
        end if

        write(unit, '(a)') '# Energy(eV),Total,s,p,d,f'
        do ie = 1, ne
            write(unit, '(f12.6,5(a,es14.6))') &
                energy_grid(ie), ',', pdos_data(ie, 1), ',', pdos_data(ie, 2), &
                ',', pdos_data(ie, 3), ',', pdos_data(ie, 4), ',', pdos_data(ie, 5)
        end do

        close(unit)
    end subroutine write_pdos_csv

end module dos_plotter
