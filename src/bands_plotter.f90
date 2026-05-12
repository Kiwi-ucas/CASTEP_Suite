module bands_plotter
    !! Band structure visualization: ASCII terminal plot + SVG output
    !! Zero external dependencies — pure Fortran 2008
    use castep_config, only: dp, bands_data_t, HARTREE_TO_EV, &
        IO_SUCCESS, IO_WRITE_FAIL, int2str
    use term_utils, only: C_RED, C_GREEN, C_YELLOW, C_CYAN, C_BOLD, C_DIM, C_RESET, C_AXIS, &
        get_term_size, draw_line
    implicit none
    private

    integer, parameter, public :: &
        BANDS_MODE_ASCII   = 1, &
        BANDS_MODE_SVG     = 2

    public :: plot_bands_ascii, write_bands_svg

contains

    ! ----------------------------------------------------------------
    !  ASCII terminal plot (auto-sized, ANSI colors)
    ! ----------------------------------------------------------------
    subroutine plot_bands_ascii(bands, term_w_in, term_h_in, e_center, half_range, &
            k_pct_in, k_width_pct_in)
        type(bands_data_t), intent(in) :: bands
        integer, intent(in) :: term_w_in, term_h_in
        real(dp), intent(in) :: half_range
        real(dp), intent(inout) :: e_center
        real(dp), intent(in), optional :: k_pct_in, k_width_pct_in

        character(len=1), allocatable :: grid(:,:)
        real(dp), allocatable :: ev_sym(:,:), kdist_sym(:)
        real(dp) :: e_min, e_max, e_range, k_max, e_val, band_step
        real(dp) :: x_scale, y_scale, fermi_ev
        integer  :: ix, iy, ik, ie, nbands, nk, nk_sym, nh, nw
        integer  :: fermi_row, last_ix, last_iy, ib_plot, n_plot
        integer  :: label_interval, label_width, gap_width, k_max_len
        character(len=32)  :: fmt_label, tmp_str, sym_char
        integer  :: band_symbol

        ! VBM/CBM detection
        real(dp) :: vbm_energy, cbm_energy, direct_cbm
        real(dp) :: gap_direct, gap_indirect
        integer  :: k_vbm, k_cbm, band_vbm, band_cbm
        integer  :: ix_vbm, iy_vbm, ix_cbm, iy_cbm
        logical  :: has_gap
        character(len=16) :: gap_label

        ! k-point path boundary detection
        integer, allocatable :: boundary_kpts(:)
        character(len=8), allocatable :: kpt_labels(:)
        integer :: nbounds, ib
        real(dp) :: dir1(3), dir2(3), d1, d2, cos_angle
        real(dp) :: kcart_prev(3), kcart_curr(3), kcart_next(3)
        real(dp) :: cell_inv(3,3), det
        integer  :: ii
        character(len=256) :: label_line
        integer  :: lbl_x, lbl_start
        real(dp) :: k_wmin, k_wmax, k_pct, k_width_pct
        nk = bands%num_kpoints
        nbands = bands%num_eigenvalues
        if (nk < 1 .or. nbands < 1) then
            write(*, '(a)') '  No band data to plot.'
            return
        end if

        ! auto-detect terminal size
        call get_term_size(nw, nh)
        if (term_w_in > 0) nw = term_w_in
        if (term_h_in > 0) nh = term_h_in
        nw = max(30, min(300, nw))
        nh = max(15, min(100, nh))
        nh = max(15, nh - 9)  ! reserve space for headers + footer

        ! ===== Part 1: symmetric band structure =====
        if (nk >= 2) then
            nk_sym = nk + 1
            allocate(ev_sym(nbands, nk_sym))
            allocate(kdist_sym(nk_sym))
            ev_sym(:, 1:nk) = bands%eigenvalues(:, 1:nk, 1)
            ev_sym(:, nk_sym) = bands%eigenvalues(:, 1, 1)
            kdist_sym(1:nk) = bands%kpath_dist(1:nk)
            kdist_sym(nk_sym) = bands%kpath_dist(nk) + (bands%kpath_dist(2) - bands%kpath_dist(1))
        else
            nk_sym = nk
            allocate(ev_sym(nbands, nk_sym))
            allocate(kdist_sym(nk_sym))
            ev_sym(:, 1:nk) = bands%eigenvalues(:, 1:nk, 1)
            kdist_sym(1:nk) = bands%kpath_dist(1:nk)
        end if

        ! fermi energy in eV
        fermi_ev = bands%fermi_energy * HARTREE_TO_EV

        ! find VBM (highest occupied) and CBM (lowest unoccupied) in eV
        vbm_energy = -huge(1.0_dp)
        cbm_energy =  huge(1.0_dp)
        k_vbm = 1; k_cbm = 1
        band_vbm = 1; band_cbm = 1
        do ik = 1, nk_sym
            do ie = 1, nbands
                e_val = ev_sym(ie, ik) * HARTREE_TO_EV
                if (e_val < fermi_ev .and. e_val > vbm_energy) then
                    vbm_energy = e_val; k_vbm = ik; band_vbm = ie
                end if
                if (e_val > fermi_ev .and. e_val < cbm_energy) then
                    cbm_energy = e_val; k_cbm = ik; band_cbm = ie
                end if
            end do
        end do

        has_gap = (vbm_energy > -huge(1.0_dp) .and. cbm_energy < huge(1.0_dp))
        if (has_gap) then
            gap_indirect = cbm_energy - vbm_energy
            direct_cbm = huge(1.0_dp)
            do ie = 1, nbands
                e_val = ev_sym(ie, k_vbm) * HARTREE_TO_EV
                if (e_val > fermi_ev .and. e_val < direct_cbm) direct_cbm = e_val
            end do
            gap_direct = direct_cbm - vbm_energy
        else
            gap_indirect = 0.0_dp; gap_direct = 0.0_dp
        end if

        ! energy range: e_center +/- half_range
        e_min = e_center - half_range
        e_max = e_center + half_range
        e_range = e_max - e_min
        if (e_range < 1.0e-12_dp) e_range = 1.0_dp

        k_max = kdist_sym(nk_sym)
        if (k_max < 1.0e-12_dp) k_max = 1.0_dp

        ! k-path window (percentage-based, for horizontal scrolling)
        k_pct = 0.5_dp
        k_width_pct = 1.0_dp
        if (present(k_pct_in)) k_pct = k_pct_in
        if (present(k_width_pct_in)) k_width_pct = k_width_pct_in
        k_wmin = max(0.0_dp, (k_pct - k_width_pct / 2.0_dp) * k_max)
        k_wmax = min(k_max, (k_pct + k_width_pct / 2.0_dp) * k_max)

        ! compute label width early so grid width fits terminal
        label_width = 0
        write(tmp_str, '(f10.4)') e_max
        label_width = max(label_width, len_trim(adjustl(tmp_str)))
        write(tmp_str, '(f10.4)') e_min
        label_width = max(label_width, len_trim(adjustl(tmp_str)))
        label_width = max(label_width, 6)
        gap_width = 1
        nw = max(20, nw - label_width - gap_width)

        ! ===== Part 2: detect k-path segment boundaries =====
        allocate(boundary_kpts(nk_sym))
        allocate(kpt_labels(nk_sym))
        nbounds = 0

        ! invert cell vectors for Cartesian k-coordinate conversion
        cell_inv = bands%cell_vectors
        det = cell_inv(1,1) * (cell_inv(2,2)*cell_inv(3,3) - cell_inv(2,3)*cell_inv(3,2)) &
            - cell_inv(1,2) * (cell_inv(2,1)*cell_inv(3,3) - cell_inv(2,3)*cell_inv(3,1)) &
            + cell_inv(1,3) * (cell_inv(2,1)*cell_inv(3,2) - cell_inv(2,2)*cell_inv(3,1))
        if (abs(det) < 1.0e-12_dp) det = 1.0_dp
        ! simple 3x3 inverse using cofactors
        cell_inv(1,1) = (bands%cell_vectors(2,2)*bands%cell_vectors(3,3) - bands%cell_vectors(2,3)*bands%cell_vectors(3,2)) / det
        cell_inv(1,2) = -(bands%cell_vectors(1,2)*bands%cell_vectors(3,3) - bands%cell_vectors(1,3)*bands%cell_vectors(3,2)) / det
        cell_inv(1,3) = (bands%cell_vectors(1,2)*bands%cell_vectors(2,3) - bands%cell_vectors(1,3)*bands%cell_vectors(2,2)) / det
        cell_inv(2,1) = -(bands%cell_vectors(2,1)*bands%cell_vectors(3,3) - bands%cell_vectors(2,3)*bands%cell_vectors(3,1)) / det
        cell_inv(2,2) = (bands%cell_vectors(1,1)*bands%cell_vectors(3,3) - bands%cell_vectors(1,3)*bands%cell_vectors(3,1)) / det
        cell_inv(2,3) = -(bands%cell_vectors(1,1)*bands%cell_vectors(2,3) - bands%cell_vectors(1,3)*bands%cell_vectors(2,1)) / det
        cell_inv(3,1) = (bands%cell_vectors(2,1)*bands%cell_vectors(3,2) - bands%cell_vectors(2,2)*bands%cell_vectors(3,1)) / det
        cell_inv(3,2) = -(bands%cell_vectors(1,1)*bands%cell_vectors(3,2) - bands%cell_vectors(1,2)*bands%cell_vectors(3,1)) / det
        cell_inv(3,3) = (bands%cell_vectors(1,1)*bands%cell_vectors(2,2) - bands%cell_vectors(1,2)*bands%cell_vectors(2,1)) / det

        ! first k-point is always a boundary
        nbounds = 1
        boundary_kpts(1) = 1

        ! detect direction changes along k-path (skip symmetric tail)
        do ik = 2, nk-1
            ! Cartesian k at ik-1, ik, ik+1
            do ii = 1, 3
                kcart_prev(ii) = cell_inv(ii,1) * bands%kpoint_coords(1,ik-1) &
                                + cell_inv(ii,2) * bands%kpoint_coords(2,ik-1) &
                                + cell_inv(ii,3) * bands%kpoint_coords(3,ik-1)
                kcart_curr(ii) = cell_inv(ii,1) * bands%kpoint_coords(1,ik) &
                                + cell_inv(ii,2) * bands%kpoint_coords(2,ik) &
                                + cell_inv(ii,3) * bands%kpoint_coords(3,ik)
                kcart_next(ii) = cell_inv(ii,1) * bands%kpoint_coords(1,ik+1) &
                                + cell_inv(ii,2) * bands%kpoint_coords(2,ik+1) &
                                + cell_inv(ii,3) * bands%kpoint_coords(3,ik+1)
            end do
            dir1 = kcart_curr - kcart_prev
            dir2 = kcart_next - kcart_curr
            d1 = sqrt(dir1(1)**2 + dir1(2)**2 + dir1(3)**2)
            d2 = sqrt(dir2(1)**2 + dir2(2)**2 + dir2(3)**2)
            if (d1 > 1.0e-12_dp .and. d2 > 1.0e-12_dp) then
                cos_angle = (dir1(1)*dir2(1) + dir1(2)*dir2(2) + dir1(3)*dir2(3)) / (d1 * d2)
                if (cos_angle < 0.966_dp) then  ! angle > 15 degrees
                    nbounds = nbounds + 1
                    boundary_kpts(nbounds) = ik
                end if
            end if
        end do

        ! last original k-point is always a boundary
        if (boundary_kpts(nbounds) /= nk) then
            nbounds = nbounds + 1
            boundary_kpts(nbounds) = nk
        end if

        ! generate labels
        do ib = 1, nbounds
            write(kpt_labels(ib), '(a,i0)') 'K', ib
        end do

        ! downsample bands if too many for terminal height
        n_plot = nbands
        if (nbands > nh) n_plot = nh
        band_step = real(nbands - 1, dp) / real(max(1, n_plot - 1), dp)

        ! allocate grid: col 1=left border, cols 2..nw-1=data, col nw=right border
        allocate(grid(nw, nh))
        grid = ' '

        ! map fermi level row
        fermi_row = nh - int((fermi_ev - e_min) / e_range * (nh - 1))
        fermi_row = max(1, min(nh, fermi_row))

        ! x_scale maps kpath_dist window to [0, nw-3], data cols 2..nw-1
        x_scale = real(nw - 3, dp) / (k_wmax - k_wmin)
        y_scale = real(nh - 1, dp) / e_range

        ! plot sampled bands with connecting lines and multi-symbol
        do ib_plot = 1, n_plot
            ie = 1 + nint((ib_plot - 1) * band_step)
            ie = max(1, min(nbands, ie))
            band_symbol = mod(ib_plot - 1, 5)
            write(sym_char, '(i1)') band_symbol
            last_ix = -1
            do ik = 1, nk_sym
                e_val = ev_sym(ie, ik) * HARTREE_TO_EV
                ix = nint((kdist_sym(ik) - k_wmin) * x_scale) + 2
                iy = nh - nint((e_val - e_min) * y_scale)
                ix = max(2, min(nw-1, ix)); iy = max(1, min(nh, iy))
                grid(ix, iy) = sym_char(1:1)
                if (last_ix > 0) call draw_line(nw, nh, grid, last_ix, last_iy, ix, iy)
                last_ix = ix; last_iy = iy
            end do
        end do

        ! VBM/CBM markers on grid (overwrite band points)
        if (has_gap) then
            ix_vbm = nint((kdist_sym(k_vbm) - k_wmin) * x_scale) + 2
            iy_vbm = nh - nint((vbm_energy - e_min) * y_scale)
            ix_vbm = max(2, min(nw-1, ix_vbm)); iy_vbm = max(1, min(nh, iy_vbm))
            grid(ix_vbm, iy_vbm) = 'V'

            ix_cbm = nint((kdist_sym(k_cbm) - k_wmin) * x_scale) + 2
            iy_cbm = nh - nint((cbm_energy - e_min) * y_scale)
            ix_cbm = max(2, min(nw-1, ix_cbm)); iy_cbm = max(1, min(nh, iy_cbm))
            grid(ix_cbm, iy_cbm) = 'C'
        end if

        ! draw fermi level (cols 2..nw-1)
        do ix = 2, nw-1
            if (grid(ix, fermi_row) == ' ') grid(ix, fermi_row) = '-'
        end do

        ! draw left and right borders
        do iy = 1, nh
            if (grid(1, iy) == ' ') grid(1, iy) = '|'
            if (grid(nw, iy) == ' ') grid(nw, iy) = '|'
        end do

        ! draw bottom border (x-axis, cols 2..nw-1)
        do ix = 2, nw-1
            if (grid(ix, nh) == ' ' .or. grid(ix, nh) == '-') grid(ix, nh) = '-'
        end do

        ! corners and junctions (top border is printed separately)
        grid(1, nh) = 'B'          ! bottom-left └
        grid(nw, nh) = 'R'         ! bottom-right ┘
        grid(1, fermi_row) = 'L'   ! Fermi left junction ├
        grid(nw, fermi_row) = 'J'  ! Fermi right junction ┤

        ! x-axis ticks at boundary k-points
        do ib = 1, nbounds
            ix = nint((kdist_sym(boundary_kpts(ib)) - k_wmin) * x_scale) + 2
            ix = max(2, min(nw-1, ix))
            if (grid(ix, nh) == '-') grid(ix, nh) = 'X'
        end do

        ! --- print ---
        write(*, '(a)') ''
        write(*, '(a,i0,a,i0,a,i0,a)') C_BOLD // '  Band Structure  ' // &
            C_RESET // C_CYAN, nk, C_RESET // ' k-points, ' // C_CYAN, &
            n_plot, C_RESET // '/' // C_CYAN, nbands, C_RESET // ' bands'
        if (nk >= 2) write(*, '(a)') '  ' // C_DIM // '(symmetric display)' // C_RESET

        write(tmp_str, '(f8.4)') e_min
        write(*, '(a)', advance='no') '  ' // C_DIM // 'Window: [' // C_RESET // &
            trim(adjustl(tmp_str)) // C_DIM // ' to ' // C_RESET
        write(tmp_str, '(f8.4)') e_max
        write(fmt_label, '(f6.2)') e_range
        write(*, '(a,i0,1x,i0,a)') trim(adjustl(tmp_str)) // C_DIM // '] eV' // &
            C_RESET // '  term:' // C_RESET, nw, nh, &
            '  ' // C_DIM // 'range=' // C_RESET // trim(adjustl(fmt_label)) // ' eV'

        label_interval = max(1, nh / 8)

        ! band gap info above plot
        if (has_gap) then
            if (abs(gap_direct - gap_indirect) < 0.005_dp .or. k_vbm == k_cbm) then
                gap_label = 'direct'
            else
                gap_label = 'indirect'
            end if
            write(fmt_label, '(f8.4)') gap_indirect
            write(*, '(a,f8.4,a,a)') '  ' // C_RED // 'E_F = ' // C_RESET, &
                fermi_ev, ' eV  ' // C_DIM // '(dashed)' // C_RESET &
                // '  |  ' // C_GREEN // 'Band Gap = ' // C_RESET // C_BOLD &
                // trim(adjustl(fmt_label)) // ' eV  (' // trim(gap_label) // ')' // C_RESET
            write(fmt_label, '(f8.4)') vbm_energy
            write(tmp_str, '(f8.4)') cbm_energy
            if (gap_label == 'indirect') then
                write(sym_char, '(f8.4)') gap_direct
                write(*, '(a,i0,a,i0,a)') &
                    '  VBM: ' // trim(adjustl(fmt_label)) // ' eV (k-pt ', k_vbm, &
                    ')  |  CBM: ' // trim(adjustl(tmp_str)) // ' eV (k-pt ', k_cbm, &
                    ')  |  Dir: ' // trim(adjustl(sym_char)) // ' eV'
            else
                write(*, '(a,i0,a,i0,a)') &
                    '  VBM: ' // trim(adjustl(fmt_label)) // ' eV (k-pt ', k_vbm, &
                    ')  |  CBM: ' // trim(adjustl(tmp_str)) // ' eV (k-pt ', k_cbm, ')'
            end if
        else
            write(*, '(a,f8.4,a)') '  ' // C_RED // 'E_F = ' // C_RESET, &
                fermi_ev, ' eV  |  ' // C_YELLOW // 'Metallic (no band gap)' // C_RESET
        end if
        write(*, '(a)') ''

        ! top border line
        write(*, '(a)') repeat(' ', label_width + gap_width) // C_AXIS &
            // '┌' // repeat('─', nw-2) // '┐' // C_RESET

        do iy = 1, nh
            fmt_label = ''
            if (mod(iy-1, label_interval) == 0 .or. iy == fermi_row .or. &
                iy == 1 .or. iy == nh) then
                if (iy == fermi_row) then
                    write(fmt_label, '(f10.4)') fermi_ev
                else
                    e_val = e_max - real(iy-1, dp) / real(nh-1, dp) * e_range
                    write(fmt_label, '(f10.4)') e_val
                end if
                fmt_label = adjustl(fmt_label)
            end if
            call print_grid_row(grid, nw, iy, fermi_row, nh, &
                fmt_label, label_width, gap_width)
        end do

        ! x-axis distance labels (window min at left, window max at right)
        write(tmp_str, '(f12.6)') k_wmax
        tmp_str = adjustl(tmp_str)
        k_max_len = len_trim(tmp_str)
        write(fmt_label, '(f10.4)') k_wmin
        fmt_label = adjustl(fmt_label)
        write(*, '(a,a,a)') repeat(' ', label_width + gap_width) // C_DIM, &
            trim(fmt_label) // repeat(' ', max(0, nw - 3 - k_max_len)) &
            // trim(tmp_str), C_RESET

        ! k-point label line
        label_line = repeat(' ', label_width + gap_width)
        do ib = 1, nbounds
            ix = nint((kdist_sym(boundary_kpts(ib)) - k_wmin) * x_scale) + 2
            ix = max(2, min(nw-1, ix))
            if (ix < 1 .or. ix > nw) cycle  ! skip labels outside window
            lbl_x = label_width + gap_width + ix  ! screen column
            lbl_start = lbl_x - len_trim(kpt_labels(ib)) / 2
            lbl_start = max(1, lbl_start)
            if (lbl_start + len_trim(kpt_labels(ib)) - 1 <= len(label_line)) then
                label_line(lbl_start:lbl_start+len_trim(kpt_labels(ib))-1) = trim(kpt_labels(ib))
            end if
        end do
        write(*, '(a)') C_DIM // trim(label_line) // C_RESET

        deallocate(grid)
        deallocate(ev_sym, kdist_sym)
        deallocate(boundary_kpts, kpt_labels)
    end subroutine plot_bands_ascii

    ! detect terminal size from environment
    ! print one grid row with run-length-encoded ANSI coloring
    subroutine print_grid_row(grid, nw, iy, fermi_row, nh, label, lbl_w, gap_w)
        integer, intent(in) :: nw, iy, fermi_row, nh, lbl_w, gap_w
        character(len=1), intent(in) :: grid(nw, nh)
        character(len=*), intent(in) :: label
        character(len=4096) :: buf
        character(len=1)   :: ch
        character(len=3)   :: dc
        integer            :: run_type
        character(len=16)  :: cur_color, clr
        integer :: ix, run_start, run_len, labellen, bp, n, slen
        logical :: is_fermi

        is_fermi = (iy == fermi_row)
        bp = 1       ! buffer position (avoids trim which strips spaces)

        ! right-aligned label prefix + gap
        labellen = len_trim(label)
        if (labellen > 0) then
            n = max(0, lbl_w - labellen)
            if (n > 0) then; buf(bp:bp+n-1) = repeat(' ', n); bp = bp + n; end if
            if (is_fermi) then; clr = C_RED; else; clr = C_DIM; end if
            slen = len_trim(clr)
            buf(bp:bp+slen-1) = trim(clr); bp = bp + slen
            buf(bp:bp+labellen-1) = label(1:labellen); bp = bp + labellen
            slen = len_trim(C_RESET)
            buf(bp:bp+slen-1) = C_RESET; bp = bp + slen
            buf(bp:bp+gap_w-1) = repeat(' ', gap_w); bp = bp + gap_w
        else
            n = lbl_w + gap_w
            buf(bp:bp+n-1) = repeat(' ', n); bp = bp + n
        end if

        ! build grid content with position-based writes (no trim)
        if (nw >= 1) then
            run_start = 1
            do while (run_start <= nw)
                ch = grid(run_start, iy)
                run_type = char_type(ch, is_fermi, iy == nh)
                run_len = 1
                do ix = run_start + 1, nw
                    if (char_type(grid(ix, iy), is_fermi, iy == nh) == run_type) then
                        run_len = run_len + 1
                    else
                        exit
                    end if
                end do
                cur_color = type_color(run_type)
                if (len_trim(cur_color) > 0) then
                    clr = trim(cur_color)
                    slen = len_trim(clr)
                    buf(bp:bp+slen-1) = clr; bp = bp + slen
                end if
                do ix = 1, run_len
                    dc = type_char(grid(run_start + ix - 1, iy), run_type)
                    slen = len_trim(dc)
                    if (slen == 0) slen = 1   ! space character
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
    end subroutine print_grid_row

    ! classify a grid cell into a run type
    pure function char_type(ch, is_fermi, is_bottom) result(ct)
        character(len=1), intent(in) :: ch
        logical, intent(in) :: is_fermi, is_bottom
        integer :: ct
        ct = 0
        if (ch == '|' .or. ch == 'B' .or. ch == 'R') then
            ct = 1     ! border / corners
        else if (ch == '-') then
            if (is_bottom) then; ct = 5         ! x-axis (check before fermi)
            else if (is_fermi) then; ct = 2     ! fermi
            end if
        else if (ch == 'L' .or. ch == 'J') then
            ct = 2      ! fermi junctions
        else if (ch == '*') then
            ct = 3      ! band point (legacy)
        else if (ch >= '0' .and. ch <= '9') then
            ct = 3      ! band point (multi-symbol)
        else if (ch == '.') then
            ct = 4      ! connector
        else if (ch == 'V') then
            ct = 6      ! VBM marker
        else if (ch == 'C') then
            ct = 7      ! CBM marker
        else if (ch == 'X') then
            ct = 5      ! x-axis tick
        end if
    end function char_type

    pure function type_color(ct) result(clr)
        integer, intent(in) :: ct
        character(len=16) :: clr
        select case (ct)
        case (1,5);  clr = C_AXIS
        case (2);    clr = C_RED
        case (3,8,9,10,11,12); clr = C_CYAN
        case (4);    clr = C_DIM
        case (6);    clr = C_GREEN
        case (7);    clr = C_YELLOW
        case default; clr = ''
        end select
    end function type_color

    pure function type_char(ch, ct) result(dc)
        character(len=1), intent(in) :: ch
        integer, intent(in) :: ct
        character(len=3) :: dc
        select case (ct)
        case (1)  ! border / y-axis
            select case (ch)
            case ('B'); dc = '└'  ! bottom-left corner
            case ('R'); dc = '┘'  ! bottom-right corner
            case default; dc = '│'  ! vertical border
            end select
        case (2)  ! fermi level
            select case (ch)
            case ('L'); dc = '├'  ! left junction
            case ('J'); dc = '┤'  ! right junction
            case default; dc = '─'
            end select
        case (3)  ! band point (multi-symbol)
            select case (ch)
            case ('0'); dc = '●'
            case ('1'); dc = '○'
            case ('2'); dc = '□'
            case ('3'); dc = '△'
            case ('4'); dc = '▽'
            case default; dc = '●'
            end select
        case (4);  dc = '·'  ! connector
        case (5)  ! x-axis
            if (ch == 'X') then; dc = '┴'   ! tick mark
            else;               dc = '─'    ! horizontal
            end if
        case (6);  dc = '◆'  ! VBM
        case (7);  dc = '◈'  ! CBM
        case default; dc = ' '
        end select
    end function type_char



    ! ----------------------------------------------------------------
    !  SVG output
    ! ----------------------------------------------------------------
    subroutine write_bands_svg(bands, svg_file, w_px, h_px, iostat, iomsg)
        type(bands_data_t), intent(in)  :: bands
        character(len=*), intent(in)    :: svg_file
        integer, intent(in)             :: w_px, h_px
        integer, intent(out)            :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer, parameter :: MARGIN_LEFT  = 80, MARGIN_RIGHT  = 20
        integer, parameter :: MARGIN_TOP   = 40, MARGIN_BOTTOM = 60
        real(dp) :: e_min, e_max, e_range, k_max, e_val, fermi_ev
        real(dp) :: kx_scale, ky_scale
        integer  :: unit, ios, ik, ie, nk, nbands
        integer  :: plot_w, plot_h, x, y
        integer  :: fermi_y
        character(len=64)   :: fermi_label
        character(len=2048) :: points, line_buf
        character(len=32)   :: rgb

        iostat = IO_SUCCESS
        if (present(iomsg)) iomsg = ''
        nk = bands%num_kpoints
        nbands = bands%num_eigenvalues
        if (nk < 1 .or. nbands < 1) then
            iostat = 1
            if (present(iomsg)) iomsg = 'No band data to plot'
            return
        end if

        open(newunit=unit, file=trim(svg_file), status='replace', &
             action='write', iostat=ios)
        if (ios /= 0) then
            iostat = IO_WRITE_FAIL
            if (present(iomsg)) iomsg = 'Cannot write: ' // trim(svg_file)
            return
        end if

        plot_w = w_px - MARGIN_LEFT - MARGIN_RIGHT
        plot_h = h_px - MARGIN_TOP  - MARGIN_BOTTOM

        fermi_ev = bands%fermi_energy * HARTREE_TO_EV

        ! energy range: Fermi +/- 10 eV, 5% padding
        e_min = fermi_ev - 10.0_dp
        e_max = fermi_ev + 10.0_dp
        e_range = e_max - e_min
        if (e_range < 1.0e-12_dp) e_range = 1.0_dp
        e_min = e_min - 0.05_dp * e_range
        e_max = e_max + 0.05_dp * e_range
        e_range = e_max - e_min

        k_max = bands%kpath_dist(nk)
        if (k_max < 1.0e-12_dp) k_max = 1.0_dp
        kx_scale = real(plot_w, dp) / k_max
        ky_scale = real(plot_h, dp) / e_range

        fermi_y = MARGIN_TOP + nint((e_max - fermi_ev) * ky_scale)

        ! --- SVG header ---
        write(unit, '(a)') '<?xml version="1.0" encoding="UTF-8"?>'
        write(unit, '(a,i0,a,i0,a,i0,a,i0,a)') &
            '<svg xmlns="http://www.w3.org/2000/svg" width="', w_px, &
            '" height="', h_px, &
            '" viewBox="0 0 ', w_px, ' ', h_px, '">'
        write(unit, '(a)') '  <rect width="100%" height="100%" fill="white"/>'

        ! --- plot area background ---
        write(unit, '(a,i0,a,i0,a,i0,a,i0,a)') &
            '  <rect x="', MARGIN_LEFT, '" y="', MARGIN_TOP, &
            '" width="', plot_w, '" height="', plot_h, &
            '" fill="none" stroke="black" stroke-width="1"/>'

        ! --- Fermi level line ---
        write(unit, '(a,i0,a,i0,a,i0,a,i0,a)') &
            '  <line x1="', MARGIN_LEFT, '" y1="', fermi_y, &
            '" x2="', MARGIN_LEFT+plot_w, '" y2="', fermi_y, &
            '" stroke="red" stroke-width="1.5" stroke-dasharray="6,4"/>'
        write(fermi_label, '(f10.4)') fermi_ev
        fermi_label = adjustl(fermi_label)
        write(unit, '(a,i0,a,i0,3a)') &
            '  <text x="', MARGIN_LEFT+plot_w+5, '" y="', fermi_y+5, &
            '" fill="red" font-size="12">E_F=', trim(fermi_label), '</text>'

        ! --- bands as polylines ---
        do ie = 1, nbands
            points = ''
            do ik = 1, nk
                x = MARGIN_LEFT + nint(bands%kpath_dist(ik) * kx_scale)
                e_val = bands%eigenvalues(ie, ik, 1) * HARTREE_TO_EV
                y = MARGIN_TOP + nint((e_max - e_val) * ky_scale)
                x = max(MARGIN_LEFT, min(MARGIN_LEFT+plot_w, x))
                y = max(MARGIN_TOP, min(MARGIN_TOP+plot_h, y))
                write(line_buf, '(i0,a,i0)') x, ',', y
                line_buf = trim(adjustl(line_buf))
                if (ik > 1) then
                    points = trim(points) // ' ' // trim(line_buf)
                else
                    points = trim(line_buf)
                end if
            end do
            call band_color(ie, nbands, rgb)
            write(unit, '(a,a,a,a,a)') &
                '  <polyline points="', trim(points), &
                '" fill="none" stroke="', trim(rgb), &
                '" stroke-width="0.8"/>'
        end do

        ! --- axis labels ---
        ! y-axis
        write(unit, '(a,i0,a,i0,a)') &
            '  <text x="15" y="', h_px/2, &
            '" transform="rotate(-90,15,', h_px/2, &
            ')" text-anchor="middle" font-size="14">Energy (eV)</text>'
        ! x-axis
        write(unit, '(a,i0,a,i0,a)') &
            '  <text x="', w_px/2, '" y="', h_px-15, &
            '" text-anchor="middle" font-size="14">k-path distance</text>'

        ! y-axis tick labels
        call write_svg_text(unit, MARGIN_LEFT-5, MARGIN_TOP, &
            real2str_short(e_max), 'end', 'black')
        call write_svg_text(unit, MARGIN_LEFT-5, MARGIN_TOP+plot_h, &
            real2str_short(e_min), 'end', 'black')

        ! title
        write(unit, '(a,i0,a,i0,a,a,i0,a,i0,a)') &
            '  <text x="', w_px/2, '" y="', MARGIN_TOP-10, &
            '" text-anchor="middle" font-size="16" font-weight="bold">', &
            'Band Structure (', nk, ' k-points, ', nbands, ' bands)</text>'

        write(unit, '(a)') '</svg>'
        close(unit)

        write(*, '(a,i0,a,i0,a)') '  SVG written: ', nk, &
            ' k-points, ', nbands, ' bands'
    end subroutine write_bands_svg

    subroutine band_color(iband, ntotal, rgb)
        integer, intent(in) :: iband, ntotal
        character(len=*), intent(out) :: rgb
        integer :: r, g, b
        real(dp) :: hue
        hue = real(iband-1, dp) / real(max(1, ntotal-1), dp)
        call hsl_to_rgb(hue, 0.7_dp, 0.5_dp, r, g, b)
        write(rgb, '(a,i0,a,i0,a,i0,a)') 'rgb(', r, ',', g, ',', b, ')'
    end subroutine band_color

    subroutine hsl_to_rgb(h, s, l, r, g, b)
        real(dp), intent(in)  :: h, s, l
        integer,  intent(out) :: r, g, b
        real(dp) :: c, x, m, r1, g1, b1, hp
        c = (1.0_dp - abs(2.0_dp*l - 1.0_dp)) * s
        hp = h * 6.0_dp
        x = c * (1.0_dp - abs(mod(hp, 2.0_dp) - 1.0_dp))
        if (hp < 1.0_dp) then
            r1 = c; g1 = x; b1 = 0.0_dp
        else if (hp < 2.0_dp) then
            r1 = x; g1 = c; b1 = 0.0_dp
        else if (hp < 3.0_dp) then
            r1 = 0.0_dp; g1 = c; b1 = x
        else if (hp < 4.0_dp) then
            r1 = 0.0_dp; g1 = x; b1 = c
        else if (hp < 5.0_dp) then
            r1 = x; g1 = 0.0_dp; b1 = c
        else
            r1 = c; g1 = 0.0_dp; b1 = x
        end if
        m = l - c / 2.0_dp
        r = nint((r1 + m) * 255.0_dp)
        g = nint((g1 + m) * 255.0_dp)
        b = nint((b1 + m) * 255.0_dp)
        r = max(0, min(255, r))
        g = max(0, min(255, g))
        b = max(0, min(255, b))
    end subroutine hsl_to_rgb

    subroutine write_svg_text(unit, x, y, txt, anchor, color)
        integer, intent(in) :: unit, x, y
        character(len=*), intent(in) :: txt, anchor, color
        write(unit, '(a,i0,a,i0,a,a,a,a,3a)') &
            '  <text x="', x, '" y="', y, &
            '" text-anchor="', trim(anchor), &
            '" font-size="11" fill="', trim(color), '">', &
            trim(txt), '</text>'
    end subroutine write_svg_text

    function real2str_short(val) result(s)
        real(dp), intent(in) :: val
        character(len=32) :: s
        if (abs(val) < 1.0e-2_dp .and. abs(val) > 1.0e-12_dp) then
            write(s, '(es10.3)') val
        else if (abs(val) >= 1000.0_dp) then
            write(s, '(f10.1)') val
        else
            write(s, '(f10.4)') val
        end if
        s = adjustl(s)
    end function real2str_short

end module bands_plotter
