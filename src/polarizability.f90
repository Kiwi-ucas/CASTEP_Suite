module polarizability
    !! Static polarizability via AIMD polarization fluctuation method
    !! Combines CASTEP DFPT optical dielectric tensor with CP2K dipole trajectory
    use castep_config, only: dp, pi, MAX_LINE_LEN
    implicit none
    private

    public :: pol_data_t, parse_castep_file, &
              parse_cp2k_dipoles, unwrap_dipoles, &
              compute_static_dielectric_windowed, &
              compute_polarizability, free_pol_data
    ! Public for drift_analysis module (not compiled into main program)
    public :: EPSILON_0, KBOLTZMANN, DEBYE_TO_CM, ANG3_TO_M3, DEBYE_PER_ANG

    ! Physical constants
    real(dp), parameter :: EPSILON_0 = 8.854187817e-12_dp     ! F/m
    real(dp), parameter :: KBOLTZMANN = 1.380649e-23_dp       ! J/K
    real(dp), parameter :: DEBYE_TO_CM = 3.33564e-30_dp       ! Debye → C·m
    real(dp), parameter :: ANG3_TO_M3 = 1.0e-30_dp            ! Å³ → m³
    real(dp), parameter :: DEBYE_PER_ANG = 4.803204_dp        ! e·Å → Debye

    ! Error codes
    integer, parameter :: IO_EPS_NOT_FOUND  = 112
    integer, parameter :: IO_EPS_PARSE_ERROR = 113
    integer, parameter :: IO_DIPOLE_ERROR   = 114

    type :: pol_data_t
        ! CASTEP optical dielectric tensor ε_∞ (3×3)
        real(dp) :: eps_inf(3,3) = 0.0_dp
        ! CP2K metadata
        integer  :: n_frames = 0
        real(dp) :: temperature = 0.0_dp
        real(dp) :: volume_ang3 = 0.0_dp
        real(dp) :: cell_abc(3) = 0.0_dp  ! a, b, c in Angstrom
        ! Dipole trajectory (n_frames, 3) in Debye (unwrapped in-place)
        real(dp), allocatable :: dipoles(:,:)
        ! Results (3×3 tensors)
        real(dp) :: eps_static(3,3) = 0.0_dp
        real(dp) :: eps_ion(3,3) = 0.0_dp
        real(dp) :: alpha_static(3,3) = 0.0_dp
        real(dp) :: alpha_ion(3,3) = 0.0_dp
        real(dp) :: alpha_inf(3,3) = 0.0_dp
        ! Unwrap statistics
        integer :: n_unwraps = 0
    end type

contains

    subroutine parse_castep_file(filename, eps_inf, cell_abc, cell_volume, iostat, iomsg)
        !! Extract optical dielectric tensor ε_∞ and lattice parameters from .castep
        !! Single-pass parser — opens the file once and collects all data.
        !!
        !! Returns:
        !!   eps_inf(3,3)   — optical dielectric tensor
        !!   cell_abc(3)    — lattice parameters a, b, c (Å)
        !!   cell_volume    — cell volume (Å³) from general parallelepiped formula
        character(len=*), intent(in)  :: filename
        real(dp), intent(out)         :: eps_inf(3,3)
        real(dp), intent(out)         :: cell_abc(3)
        real(dp), intent(out)         :: cell_volume
        integer, intent(out)          :: iostat
        character(len=*), intent(out), optional :: iomsg
        character(len=MAX_LINE_LEN) :: line, lattice_lines(3)
        real(dp) :: angles(3)  ! alpha, beta, gamma in degrees
        real(dp) :: cos_a, cos_b, cos_g, cos2_sum, deg2rad
        integer :: unit, ios, idx, i
        logical :: found_eps, found_cell

        iostat = 0
        if (present(iomsg)) iomsg = ''
        eps_inf = 0.0_dp
        cell_abc = 0.0_dp
        cell_volume = 0.0_dp
        angles = 90.0_dp  ! default orthogonal

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            iostat = IO_EPS_NOT_FOUND
            if (present(iomsg)) iomsg = 'Cannot open file: ' // trim(filename)
            return
        end if

        found_eps = .false.
        found_cell = .false.

        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit

            ! --- Optical Permittivity (ε_∞) ---
            if (index(line, 'Optical Permittivity') > 0) then
                read(unit, '(a)', iostat=ios) line  ! skip separator
                if (ios /= 0) exit
                do i = 1, 3
                    read(unit, '(a)', iostat=ios) line
                    if (ios /= 0) exit
                    read(line, *, iostat=ios) eps_inf(i, 1:3)
                    if (ios /= 0) exit
                end do
                if (ios == 0) found_eps = .true.
                cycle
            end if

            ! --- Lattice parameters ---
            if (index(line, 'Lattice parameters') > 0) then
                ! Read and save the 3 lattice lines for dual parsing
                do i = 1, 3
                    read(unit, '(a)', iostat=ios) lattice_lines(i)
                    if (ios /= 0) exit
                end do
                if (ios /= 0) cycle

                ! Parse cell lengths and angles from saved lines
                ! Format: " a =   10.128028   alpha =  90.000000"
                !         " b =   10.128028   beta  =  90.000000"
                !         " c =   10.128028   gamma =  90.000000"
                do i = 1, 3
                    line = lattice_lines(i)
                    ! Cell length: use " a =" / " b =" / " c =" (leading space)
                    if (i == 1) idx = index(line, ' a =')
                    if (i == 2) idx = index(line, ' b =')
                    if (i == 3) idx = index(line, ' c =')
                    if (idx > 0) read(line(idx+4:), *, iostat=ios) cell_abc(i)
                    if (ios /= 0) cycle

                    ! Angle: look for "alpha =" / "beta  =" / "gamma ="
                    if (i == 1) idx = index(line, 'alpha =')
                    if (i == 2) idx = index(line, 'beta  =')
                    if (i == 3) idx = index(line, 'gamma =')
                    if (idx > 0) read(line(idx+7:), *, iostat=ios) angles(i)
                end do
                if (any(cell_abc > 0.0_dp)) found_cell = .true.
                cycle
            end if
        end do
        close(unit)

        if (.not. found_eps) then
            iostat = IO_EPS_PARSE_ERROR
            if (present(iomsg)) iomsg = 'No optical dielectric tensor found in: ' // trim(filename)
            return
        end if

        if (.not. found_cell) then
            iostat = IO_EPS_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Cell parameters not found in: ' // trim(filename)
            return
        end if

        ! Compute general parallelepiped volume
        deg2rad = acos(-1.0_dp) / 180.0_dp
        cos_a = cos(angles(1) * deg2rad)
        cos_b = cos(angles(2) * deg2rad)
        cos_g = cos(angles(3) * deg2rad)
        cos2_sum = cos_a**2 + cos_b**2 + cos_g**2

        cell_volume = cell_abc(1) * cell_abc(2) * cell_abc(3) * &
                      sqrt(max(0.0_dp, 1.0_dp - cos2_sum + 2.0_dp * cos_a * cos_b * cos_g))
    end subroutine parse_castep_file


    subroutine parse_cp2k_dipoles(dir_path, data, iostat, iomsg)
        !! Read CP2K dipole files from directory
        !! Each file has format:
        !!   Dipole moment [Debye]
        !!     X= xxxxx Y= yyyyy Z= zzzzz  Total= ttttt
        character(len=*), intent(in)  :: dir_path
        type(pol_data_t), intent(inout) :: data
        integer, intent(out)          :: iostat
        character(len=*), intent(out), optional :: iomsg
        character(len=MAX_LINE_LEN) :: line, fname, cmd
        character(len=32) :: file_list
        integer :: list_unit, dip_unit, ios, iframe, n_files
        integer :: cmd_stat, cmd_exit
        real(dp) :: dx, dy, dz

        iostat = 0
        if (present(iomsg)) iomsg = ''

        ! Generate file list via `find` (not shell glob).
        ! `ls dir/*cell_dipole*` expands the glob via argv, which exceeds ARG_MAX
        ! when the trajectory contains thousands of dipole files.
        ! `find` walks the directory without argv expansion.
        file_list = '.dipole_list_tmp'
        ! Use `sort -V` (version sort) to handle filenames with mixed digit
        ! widths correctly (e.g. dipole_9999.dat before dipole_10000.dat).
        ! Plain `sort` is lexicographic and breaks at digit-width boundaries.
        cmd = 'find "' // trim(dir_path) // '" -maxdepth 1 -name "*cell_dipole*" ' // &
              '2>/dev/null | sort -V > ' // trim(file_list)
        call execute_command_line(trim(cmd), wait=.true., &
                                  exitstat=cmd_exit, cmdstat=cmd_stat)
        if (cmd_stat /= 0 .or. cmd_exit /= 0) then
            iostat = IO_DIPOLE_ERROR
            if (present(iomsg)) iomsg = 'Failed to list dipole files in: ' // trim(dir_path)
            return
        end if

        open(newunit=list_unit, file=trim(file_list), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            iostat = IO_DIPOLE_ERROR
            if (present(iomsg)) iomsg = 'Cannot list dipole files in: ' // trim(dir_path)
            return
        end if

        ! Count files
        n_files = 0
        do
            read(list_unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) > 0) n_files = n_files + 1
        end do

        if (n_files == 0) then
            close(list_unit, status='delete')
            iostat = IO_DIPOLE_ERROR
            if (present(iomsg)) iomsg = 'No dipole files found in: ' // trim(dir_path)
            return
        end if

        ! Allocate arrays
        data%n_frames = n_files
        allocate(data%dipoles(n_files, 3))
        data%dipoles = 0.0_dp

        ! Read all files
        rewind(list_unit)
        do iframe = 1, n_files
            read(list_unit, '(a)', iostat=ios) fname
            if (ios /= 0) then
                iostat = IO_DIPOLE_ERROR
                if (present(iomsg)) iomsg = 'Error reading file list at frame'
                close(list_unit, status='delete')
                return
            end if

            open(newunit=dip_unit, file=trim(fname), status='old', action='read', iostat=ios)
            if (ios /= 0) then
                iostat = IO_DIPOLE_ERROR
                if (present(iomsg)) iomsg = 'Cannot open: ' // trim(fname)
                close(list_unit, status='delete')
                return
            end if

            ! Skip lines until "Dipole moment [Debye]"
            do
                read(dip_unit, '(a)', iostat=ios) line
                if (ios /= 0) exit
                if (index(line, 'Dipole moment') > 0) exit
            end do
            if (ios /= 0 .or. index(line, 'Dipole moment') == 0) then
                close(dip_unit)
                close(list_unit, status='delete')
                iostat = IO_DIPOLE_ERROR
                if (present(iomsg)) iomsg = 'Dipole moment not found in: ' // trim(fname)
                return
            end if

            ! Read dipole values: "  X= xxxxx Y= yyyyy Z= zzzzz  Total= ttttt"
            read(dip_unit, '(a)', iostat=ios) line
            if (ios /= 0) then
                close(dip_unit)
                close(list_unit, status='delete')
                iostat = IO_DIPOLE_ERROR
                if (present(iomsg)) iomsg = 'Error reading dipole values in: ' // trim(fname)
                return
            end if
            close(dip_unit)

            ! Parse X= Y= Z= values
            call parse_dipole_line(line, dx, dy, dz, ios)
            if (ios /= 0) then
                close(dip_unit)
                close(list_unit, status='delete')
                iostat = IO_DIPOLE_ERROR
                if (present(iomsg)) iomsg = 'Cannot parse dipole values in: ' // trim(fname)
                return
            end if
            data%dipoles(iframe, 1) = dx
            data%dipoles(iframe, 2) = dy
            data%dipoles(iframe, 3) = dz
        end do
        close(list_unit, status='delete')
    end subroutine parse_cp2k_dipoles


    subroutine parse_dipole_line(line, dx, dy, dz, iostat)
        !! Parse "  X= xxxxx Y= yyyyy Z= zzzzz  Total= ttttt" line
        character(len=*), intent(in) :: line
        real(dp), intent(out) :: dx, dy, dz
        integer, intent(out)  :: iostat
        character(len=MAX_LINE_LEN) :: tmp
        integer :: ix, iy, iz

        iostat = 0
        dx = 0.0_dp; dy = 0.0_dp; dz = 0.0_dp

        ! Find positions of X= Y= Z=
        ix = index(line, 'X=')
        iy = index(line, 'Y=')
        iz = index(line, 'Z=')

        if (ix == 0 .or. iy == 0 .or. iz == 0) then
            iostat = -1
            return
        end if

        ! Extract values
        read(line(ix+2:iy-1), *, iostat=iostat) dx
        if (iostat /= 0) return
        read(line(iy+2:iz-1), *, iostat=iostat) dy
        if (iostat /= 0) return
        ! Z extends to end of line or "Total="
        tmp = line(iz+2:)
        iz = index(tmp, 'Total')
        if (iz > 0) then
            read(tmp(1:iz-1), *, iostat=iostat) dz
        else
            read(tmp, *, iostat=iostat) dz
        end if
    end subroutine parse_dipole_line


    subroutine unwrap_dipoles(data)
        !! Unwrap Berry phase polarization quantum jumps
        !! Pq,α = DEBYE_PER_ANG × a_α [Debye]
        !! Compare consecutive RAW values to detect boundary crossings;
        !! use nint for multi-quantum correction.  Must compare raw-to-raw
        !! (not raw-vs-unwrapped), otherwise the comparison basis shifts
        !! after each real wrap, triggering cascading false wraps.
        type(pol_data_t), intent(inout) :: data
        real(dp) :: pq(3), delta, offset(3)
        real(dp), allocatable :: raw(:, :)
        integer :: i, j, jumps

        do j = 1, 3
            pq(j) = DEBYE_PER_ANG * data%cell_abc(j)
        end do

        ! Keep a copy of the original raw values for comparison
        allocate(raw(data%n_frames, 3))
        raw = data%dipoles

        offset = 0.0_dp
        data%n_unwraps = 0

        do i = 2, data%n_frames
            do j = 1, 3
                ! Compare consecutive raw values (before any offset)
                delta = raw(i, j) - raw(i-1, j)
                jumps = nint(delta / pq(j))
                if (jumps /= 0) then
                    offset(j) = offset(j) - real(jumps, dp) * pq(j)
                    data%n_unwraps = data%n_unwraps + abs(jumps)
                end if
                data%dipoles(i, j) = raw(i, j) + offset(j)
            end do
        end do

        deallocate(raw)
    end subroutine unwrap_dipoles


    subroutine compute_static_dielectric_windowed(data, time_step_fs, iostat, iomsg, verbose)
        !! Window-based static dielectric with extrapolation to W→0
        !! Removes Li+ diffusion via per-window detrend + median
        !! Extrapolates to zero window size for vibrational limit
        type(pol_data_t), intent(inout) :: data
        real(dp), intent(in) :: time_step_fs
        integer, intent(out) :: iostat
        character(len=*), intent(out), optional :: iomsg
        logical, intent(in), optional :: verbose

        real(dp) :: ws_ps(4) = [0.1_dp, 0.2_dp, 0.5_dp, 1.0_dp]
        integer, parameter :: N_WS = 4
        real(dp) :: vol_m3, kt, denom, dt_ps
        real(dp) :: med_eps_ion(N_WS) = 0.0_dp
        integer :: i, j, iw, n_win, nf, win_start, win_end
        real(dp) :: slope, intercept, dn
        real(dp), allocatable :: dw(:, :), cov(:, :), eps_ion(:, :)
        real(dp) :: iso_vals(1000), med_iso, sum_x, sum_y, sum_xx, sum_xy

        iostat = 0
        if (present(iomsg)) iomsg = ''

        if (data%n_frames < 2) then
            iostat = IO_DIPOLE_ERROR
            if (present(iomsg)) iomsg = 'Need at least 2 frames'
            return
        end if

        dt_ps = time_step_fs / 1000.0_dp  ! fs → ps
        vol_m3 = data%volume_ang3 * ANG3_TO_M3
        kt = KBOLTZMANN * data%temperature
        denom = EPSILON_0 * kt * vol_m3

        if (present(verbose)) then
            if (verbose) then
                write(*, '(a)') ''
                write(*, '(a)') '  Window-size convergence:'
                write(*, '(a)') '  Window (ps)   N_wins   ε_ion (median)'
                write(*, '(a)') '  ----------   ------   --------------'
            end if
        end if

        do iw = 1, N_WS
            nf = nint(ws_ps(iw) / dt_ps)
            if (nf < 4 .or. nf > data%n_frames / 2) cycle
            n_win = data%n_frames / nf
            if (n_win > 1000) n_win = 1000

            allocate(dw(nf, 3), cov(3, 3), eps_ion(3, 3))
            do i = 1, n_win
                win_start = (i - 1) * nf + 1
                win_end = win_start + nf - 1

                ! Copy window data and detrend
                do j = 1, 3
                    dw(:, j) = data%dipoles(win_start:win_end, j)
                end do
                call detrend_window(dw, nf)

                ! Covariance of this window
                call compute_covariance_single(dw, nf, cov)
                cov = cov * DEBYE_TO_CM**2
                eps_ion = cov / denom
                iso_vals(i) = (eps_ion(1,1) + eps_ion(2,2) + eps_ion(3,3)) / 3.0_dp
            end do

            ! Median over windows
            med_iso = median(iso_vals(1:n_win))
            med_eps_ion(iw) = med_iso
            if (present(verbose)) then
                if (verbose) write(*, '(a,f6.1,a,i8,a,f8.3)') '  ', ws_ps(iw), '         ', n_win, '     ', med_iso
            end if
            deallocate(dw, cov, eps_ion)
        end do

        ! Linear extrapolation to W→0 using all valid windows
        sum_x = 0.0_dp; sum_y = 0.0_dp; sum_xx = 0.0_dp; sum_xy = 0.0_dp
        dn = 0.0_dp
        do iw = 1, N_WS
            if (med_eps_ion(iw) > 0.0_dp) then
                dn = dn + 1.0_dp
                sum_x = sum_x + ws_ps(iw)
                sum_y = sum_y + med_eps_ion(iw)
                sum_xx = sum_xx + ws_ps(iw)**2
                sum_xy = sum_xy + ws_ps(iw) * med_eps_ion(iw)
            end if
        end do

        if (dn >= 2.0_dp) then
            slope = (dn * sum_xy - sum_x * sum_y) / (dn * sum_xx - sum_x**2)
            intercept = (sum_y - slope * sum_x) / dn
            if (present(verbose)) then
                if (verbose) then
                    write(*, '(a)') ''
                    write(*, '(a,f8.3)') '  ε_ion (W→0 extrapolation): ', max(0.0_dp, intercept)
                    write(*, '(a,f8.3)') '  Diffusion slope dε/dW:     ', slope
                end if
            end if

            ! Write back to pol for downstream display
            data%eps_ion = 0.0_dp
            do j = 1, 3
                data%eps_ion(j, j) = max(0.0_dp, intercept)
            end do
        end if
    end subroutine compute_static_dielectric_windowed


    subroutine detrend_window(dw, nf)
        !! Linear detrend of a single window's dipole data
        real(dp), intent(inout) :: dw(:,:)
        integer, intent(in) :: nf
        real(dp) :: sum_t, sum_t2, sum_d, sum_td, slope, intercept, dn
        integer :: i, j

        dn = real(nf, dp)
        sum_t = dn * (dn + 1.0_dp) * 0.5_dp
        sum_t2 = dn * (dn + 1.0_dp) * (2.0_dp * dn + 1.0_dp) / 6.0_dp

        do j = 1, 3
            sum_d = sum(dw(1:nf, j))
            sum_td = 0.0_dp
            do i = 1, nf
                sum_td = sum_td + real(i, dp) * dw(i, j)
            end do
            slope = (dn * sum_td - sum_t * sum_d) / (dn * sum_t2 - sum_t * sum_t)
            intercept = (sum_d - slope * sum_t) / dn
            do i = 1, nf
                dw(i, j) = dw(i, j) - (slope * real(i, dp) + intercept)
            end do
        end do
    end subroutine detrend_window


    subroutine compute_covariance_single(dw, nf, cov)
        !! Compute 3×3 covariance of detrended window data
        real(dp), intent(in) :: dw(:,:)
        integer, intent(in) :: nf
        real(dp), intent(out) :: cov(3,3)
        real(dp) :: m_mean(3)
        integer :: i, j, k

        m_mean = 0.0_dp
        do i = 1, nf
            m_mean = m_mean + dw(i, :)
        end do
        m_mean = m_mean / real(nf, dp)

        cov = 0.0_dp
        do i = 1, nf
            do j = 1, 3
                do k = 1, 3
                    cov(j, k) = cov(j, k) + (dw(i, j) - m_mean(j)) * (dw(i, k) - m_mean(k))
                end do
            end do
        end do
        cov = cov / real(nf - 1, dp)
    end subroutine compute_covariance_single


    pure function median(arr) result(m)
        !! Median of an array (in-place sorting via simple selection for small N)
        real(dp), intent(in) :: arr(:)
        real(dp) :: m
        real(dp), allocatable :: tmp(:)
        integer :: n, i, j, min_idx
        real(dp) :: temp

        n = size(arr)
        allocate(tmp(n))
        tmp = arr
        ! Selection sort to median position
        do i = 1, n / 2 + 1
            min_idx = i
            do j = i + 1, n
                if (tmp(j) < tmp(min_idx)) min_idx = j
            end do
            temp = tmp(i); tmp(i) = tmp(min_idx); tmp(min_idx) = temp
        end do
        if (mod(n, 2) == 1) then
            m = tmp(n / 2 + 1)
        else
            m = (tmp(n / 2) + tmp(n / 2 + 1)) * 0.5_dp
        end if
    end function median


    subroutine compute_polarizability(data)
        !! Convert dielectric tensor to polarizability (Å³)
        !! α_αβ = (ε_αβ − δ_αβ) × Ω / (4π)
        type(pol_data_t), intent(inout) :: data
        real(dp) :: factor, kronecker
        integer :: i, j

        factor = data%volume_ang3 / (4.0_dp * pi)

        do i = 1, 3
            do j = 1, 3
                kronecker = merge(1.0_dp, 0.0_dp, i == j)
                data%alpha_static(i, j) = (data%eps_static(i, j) - kronecker) * factor
                data%alpha_ion(i, j) = data%eps_ion(i, j) * factor
                data%alpha_inf(i, j) = (data%eps_inf(i, j) - kronecker) * factor
            end do
        end do
    end subroutine compute_polarizability


    subroutine free_pol_data(data)
        !! Deallocate arrays and reset all fields
        type(pol_data_t), intent(inout) :: data
        if (allocated(data%dipoles)) deallocate(data%dipoles)
        data%eps_inf = 0.0_dp
        data%n_frames = 0
        data%temperature = 0.0_dp
        data%volume_ang3 = 0.0_dp
        data%cell_abc = 0.0_dp
        data%eps_static = 0.0_dp
        data%eps_ion = 0.0_dp
        data%alpha_static = 0.0_dp
        data%alpha_ion = 0.0_dp
        data%alpha_inf = 0.0_dp
        data%n_unwraps = 0
    end subroutine free_pol_data

end module polarizability
