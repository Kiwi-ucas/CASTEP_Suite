module castep_vib
    !! Parse CASTEP .castep vibrational frequencies (Gamma q-pt) for IR and Raman.
    !!
    !! .castep files have explicit column headers in the "Vibrational Frequencies"
    !! section that unambiguously label IR intensity, Raman activity, and
    !! symmetry-based active flags — a cleaner source than .phonon for IR/Raman.
    use castep_config, only: dp, pi
    implicit none
    private

    public :: vib_data_t, parse_castep_vib, compute_ir_spectrum, &
        compute_raman_spectrum, free_vib_data

    integer, parameter :: MAX_LINE = 1024

    type :: vib_data_t
        integer :: n_modes = 0
        real(dp), allocatable :: freq(:)            ! frequencies (cm-1)
        character(len=4), allocatable :: irrep(:)    ! irreducible representation
        real(dp), allocatable :: ir_intensity(:)     ! IR intensity ((D/A)^2/amu)
        logical, allocatable :: ir_active(:)         ! symmetry-based IR active flag
        real(dp), allocatable :: raman_activity(:)   ! Raman activity (A^4/amu), only if computed
        logical, allocatable :: raman_active(:)      ! symmetry-based Raman active flag
        logical :: has_raman_numeric = .false.       ! true if Raman activity is numeric column
        real(dp) :: smearing = 1.0_dp
        real(dp) :: freq_grid(4001)                  ! output frequency grid (cm-1)
        real(dp), allocatable :: ir_spectrum(:)      ! broadened IR spectrum
        real(dp), allocatable :: raman_spectrum(:)   ! broadened Raman spectrum
    end type vib_data_t

contains

    subroutine parse_castep_vib(filename, data, iostat, iomsg)
        !! Parse Gamma-point vibrational frequencies from .castep file.
        !! Detects column layout from explicit headers.
        character(len=*), intent(in) :: filename
        type(vib_data_t), intent(out) :: data
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer :: unit, ios, i, n_tokens
        character(len=MAX_LINE) :: line, inner
        character(len=MAX_LINE) :: hdr1, hdr2
        logical :: found_section, found_gamma, in_data
        real(dp) :: freq_val, ir_val, ram_val
        character(len=4) :: irrep_val
        character(len=1) :: ir_active_ch, ram_active_ch
        real(dp), allocatable :: temp_freq(:), temp_ir(:), temp_ram(:)
        character(len=4), allocatable :: temp_irrep(:)
        logical, allocatable :: temp_ir_act(:), temp_ram_act(:)
        integer :: n_alloc, n_read

        iostat = 0
        found_section = .false.
        found_gamma = .false.
        data%has_raman_numeric = .false.

        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            iostat = 100
            if (present(iomsg)) iomsg = 'Cannot open file: ' // trim(filename)
            return
        end if

        ! ── Locate "Vibrational Frequencies" section ──
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'Vibrational Frequencies') > 0) then
                found_section = .true.
                exit
            end if
        end do

        if (.not. found_section) then
            iostat = 101
            if (present(iomsg)) iomsg = 'Vibrational Frequencies section not found in .castep'
            close(unit); return
        end if

        ! ── Find Gamma q-point (q-pt=1 or q-pt = 1) ──
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'q-pt=') > 0 .or. index(line, 'q-pt =') > 0) then
                ! Check if this is q-pt=1 (Gamma)
                if (index(line, 'q-pt=    1') > 0 .or. index(line, 'q-pt =    1') > 0) then
                    found_gamma = .true.
                    exit
                else
                    ! Hit a non-Gamma q-pt before finding Gamma — shouldn't happen
                    iostat = 102
                    if (present(iomsg)) iomsg = 'Gamma q-point not found before other q-pts'
                    close(unit); return
                end if
            end if
        end do

        if (.not. found_gamma) then
            iostat = 103
            if (present(iomsg)) iomsg = 'Gamma q-point (q-pt=1) not found'
            close(unit); return
        end if

        ! ── Skip to column headers (past acoustic sum rule line) ──
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) then
                iostat = 104
                if (present(iomsg)) iomsg = 'Unexpected EOF before column headers'
                close(unit); return
            end if
            ! Acoustic sum rule line contains "Acoustic sum rule"
            if (index(line, 'Acoustic sum rule') > 0) then
                ! Next non-empty lines are column headers
                read(unit, '(a)', iostat=ios) hdr1
                if (ios /= 0) then
                    iostat = 105
                    if (present(iomsg)) iomsg = 'Missing column header line 1'
                    close(unit); return
                end if
                read(unit, '(a)', iostat=ios) hdr2
                if (ios /= 0) then
                    iostat = 106
                    if (present(iomsg)) iomsg = 'Missing column header line 2'
                    close(unit); return
                end if
                exit
            end if
        end do

        ! ── Parse column headers to determine layout ──
        data%has_raman_numeric = (index(hdr1, 'raman activity') > 0 .and. &
                                   index(hdr1, 'raman activity  active') > 0)

        ! ── Skip blank line after headers ──
        read(unit, '(a)', iostat=ios) line

        ! ── Read data lines (until next section or non-data line) ──
        n_alloc = 500
        allocate(temp_freq(n_alloc), temp_irrep(n_alloc), temp_ir(n_alloc), &
                 temp_ram(n_alloc), temp_ir_act(n_alloc), temp_ram_act(n_alloc))
        temp_ir = 0.0_dp
        temp_ram = 0.0_dp
        temp_ir_act = .false.
        temp_ram_act = .false.
        n_read = 0

        in_data = .true.
        do while (in_data)
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit

            ! Check for end of Gamma data block
            if (len_trim(line) == 0) cycle  ! skip blank lines
            if (index(line, 'q-pt=') > 0 .or. index(line, 'q-pt =') > 0) exit
            if (index(line, '-----') > 0) cycle  ! skip separator lines
            if (index(line, 'Rep  Mul') > 0) exit  ! character table header
            if (index(line, 'Character table') > 0) exit

            ! Must be a data line — extract content between '+' delimiters
            inner = strip_castep_line(line)
            if (len_trim(inner) == 0) cycle

            n_tokens = count_tokens(inner)
            if (n_tokens < 3) cycle  ! need at least N, freq, irrep

            n_read = n_read + 1
            if (n_read > n_alloc) then
                call expand_arrays(temp_freq, temp_irrep, temp_ir, temp_ram, &
                                   temp_ir_act, temp_ram_act, n_alloc)
            end if

            ! Parse based on detected column layout
            if (data%has_raman_numeric .and. n_tokens >= 7) then
                read(inner, *, iostat=ios) i, freq_val, irrep_val, ir_val, &
                    ir_active_ch, ram_val, ram_active_ch
                if (ios == 0) then
                    temp_ram(n_read) = ram_val
                    temp_ram_act(n_read) = (ram_active_ch == 'Y' .or. ram_active_ch == 'y')
                end if
            else if (n_tokens >= 5) then
                ! US/NC format: N, freq, irrep, ir_intensity, ir_active, ram_active
                ! Read all tokens as scalars, handle irrep as char
                read(inner, *, iostat=ios) i, freq_val, irrep_val, ir_val, &
                    ir_active_ch, ram_active_ch
                if (ios /= 0) then
                    ! Try without ram_active_ch (5 tokens: N, freq, irrep, ir_int, ir_active)
                    read(inner, *, iostat=ios) i, freq_val, irrep_val, ir_val, ir_active_ch
                    ram_active_ch = 'N'
                end if
            else
                ! Minimal: N, freq, irrep
                read(inner, *, iostat=ios) i, freq_val, irrep_val
                ir_val = 0.0_dp
                ir_active_ch = 'N'
                ram_active_ch = 'N'
            end if

            if (ios /= 0) then
                ! Parse failed — skip this line
                n_read = n_read - 1
                cycle
            end if

            temp_freq(n_read) = freq_val
            temp_irrep(n_read) = irrep_val
            temp_ir(n_read) = ir_val
            temp_ir_act(n_read) = (ir_active_ch == 'Y' .or. ir_active_ch == 'y')
            if (.not. data%has_raman_numeric) then
                temp_ram_act(n_read) = (ram_active_ch == 'Y' .or. ram_active_ch == 'y')
            end if
        end do

        close(unit)

        if (n_read == 0) then
            iostat = 107
            if (present(iomsg)) iomsg = 'No vibrational frequency data found at Gamma'
            deallocate(temp_freq, temp_irrep, temp_ir, temp_ram, temp_ir_act, temp_ram_act)
            return
        end if

        ! ── Trim to actual size ──
        data%n_modes = n_read
        allocate(data%freq(n_read), data%irrep(n_read), data%ir_intensity(n_read), &
                 data%ir_active(n_read), data%raman_activity(n_read), data%raman_active(n_read))
        data%freq(1:n_read) = temp_freq(1:n_read)
        data%irrep(1:n_read) = temp_irrep(1:n_read)
        data%ir_intensity(1:n_read) = temp_ir(1:n_read)
        data%ir_active(1:n_read) = temp_ir_act(1:n_read)
        data%raman_activity(1:n_read) = temp_ram(1:n_read)
        data%raman_active(1:n_read) = temp_ram_act(1:n_read)

        deallocate(temp_freq, temp_irrep, temp_ir, temp_ram, temp_ir_act, temp_ram_act)

    end subroutine parse_castep_vib


    function strip_castep_line(line) result(inner)
        !! Extract content between first and last '+' on a .castep banner line
        character(len=*), intent(in) :: line
        character(len=MAX_LINE) :: inner
        integer :: i1, i2, n
        inner = ''
        n = len_trim(line)
        if (n == 0) return
        i1 = index(line(1:n), '+')
        i2 = index(line(1:n), '+', back=.true.)
        if (i1 > 0 .and. i2 > i1) then
            inner = adjustl(line(i1+1:i2-1))
        else
            inner = adjustl(line(1:n))
        end if
    end function strip_castep_line


    pure function count_tokens(line) result(n)
        !! Count whitespace-separated tokens in a line
        character(len=*), intent(in) :: line
        integer :: n, i
        logical :: in_token
        n = 0
        in_token = .false.
        do i = 1, len_trim(line)
            if (line(i:i) /= ' ') then
                if (.not. in_token) then
                    n = n + 1
                    in_token = .true.
                end if
            else
                in_token = .false.
            end if
        end do
    end function count_tokens


    subroutine expand_arrays(freq, irrep, ir, ram, ir_act, ram_act, n_alloc)
        real(dp), allocatable, intent(inout) :: freq(:), ir(:), ram(:)
        character(len=4), allocatable, intent(inout) :: irrep(:)
        logical, allocatable, intent(inout) :: ir_act(:), ram_act(:)
        integer, intent(inout) :: n_alloc
        real(dp), allocatable :: tmp_freq(:), tmp_ir(:), tmp_ram(:)
        character(len=4), allocatable :: tmp_irrep(:)
        logical, allocatable :: tmp_ir_act(:), tmp_ram_act(:)
        integer :: old_size, new_size

        old_size = n_alloc
        new_size = old_size * 2

        allocate(tmp_freq(old_size), tmp_irrep(old_size), tmp_ir(old_size), &
                 tmp_ram(old_size), tmp_ir_act(old_size), tmp_ram_act(old_size))
        tmp_freq = freq; tmp_irrep = irrep; tmp_ir = ir
        tmp_ram = ram; tmp_ir_act = ir_act; tmp_ram_act = ram_act

        deallocate(freq, irrep, ir, ram, ir_act, ram_act)
        allocate(freq(new_size), irrep(new_size), ir(new_size), &
                 ram(new_size), ir_act(new_size), ram_act(new_size))
        freq(1:old_size) = tmp_freq
        irrep(1:old_size) = tmp_irrep
        ir(1:old_size) = tmp_ir
        ram(1:old_size) = tmp_ram
        ir_act(1:old_size) = tmp_ir_act
        ram_act(1:old_size) = tmp_ram_act
        freq(old_size+1:new_size) = 0.0_dp
        ir(old_size+1:new_size) = 0.0_dp
        ram(old_size+1:new_size) = 0.0_dp
        ir_act(old_size+1:new_size) = .false.
        ram_act(old_size+1:new_size) = .false.

        deallocate(tmp_freq, tmp_irrep, tmp_ir, tmp_ram, tmp_ir_act, tmp_ram_act)
        n_alloc = new_size
    end subroutine expand_arrays


    subroutine compute_ir_spectrum(data, freq_range_min, freq_range_max, n_points, &
                                    smearing, iostat, iomsg)
        !! Gaussian-broadened IR absorption spectrum from Gamma-point intensities
        type(vib_data_t), intent(inout) :: data
        real(dp), intent(in) :: freq_range_min, freq_range_max
        integer, intent(in) :: n_points
        real(dp), intent(in) :: smearing
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        real(dp) :: freq, sigma, norm, dE, gauss, ir_val
        integer :: i, j, n_alloc
        real(dp), parameter :: GAUSS_CUTOFF = 5.0_dp

        iostat = 0

        if (.not. allocated(data%ir_intensity)) then
            iostat = 1
            if (present(iomsg)) iomsg = 'No IR intensity data loaded'
            return
        end if

        n_alloc = min(n_points, 4001)

        if (allocated(data%ir_spectrum)) deallocate(data%ir_spectrum)
        allocate(data%ir_spectrum(n_alloc))
        data%ir_spectrum = 0.0_dp
        data%smearing = smearing

        do i = 1, n_alloc
            data%freq_grid(i) = freq_range_min + &
                (freq_range_max - freq_range_min) * real(i - 1, dp) / real(n_alloc - 1, dp)
        end do

        sigma = smearing
        norm = 1.0_dp / (sigma * sqrt(2.0_dp * pi))

        do i = 1, n_alloc
            freq = data%freq_grid(i)
            do j = 1, data%n_modes
                ir_val = data%ir_intensity(j)
                if (ir_val <= 0.0_dp) cycle
                dE = freq - data%freq(j)
                if (abs(dE) > GAUSS_CUTOFF * sigma) cycle
                gauss = norm * exp(-0.5_dp * (dE / sigma)**2)
                data%ir_spectrum(i) = data%ir_spectrum(i) + gauss * ir_val
            end do
        end do
    end subroutine compute_ir_spectrum


    subroutine compute_raman_spectrum(data, freq_range_min, freq_range_max, n_points, &
                                       smearing, iostat, iomsg)
        !! Gaussian-broadened Raman scattering spectrum from Gamma-point activities
        type(vib_data_t), intent(inout) :: data
        real(dp), intent(in) :: freq_range_min, freq_range_max
        integer, intent(in) :: n_points
        real(dp), intent(in) :: smearing
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        real(dp) :: freq, sigma, norm, dE, gauss, ram_val
        integer :: i, j, n_alloc
        real(dp), parameter :: GAUSS_CUTOFF = 5.0_dp

        iostat = 0

        if (.not. allocated(data%raman_activity)) then
            iostat = 1
            if (present(iomsg)) iomsg = 'No Raman activity data loaded'
            return
        end if

        n_alloc = min(n_points, 4001)

        if (allocated(data%raman_spectrum)) deallocate(data%raman_spectrum)
        allocate(data%raman_spectrum(n_alloc))
        data%raman_spectrum = 0.0_dp
        data%smearing = smearing

        do i = 1, n_alloc
            data%freq_grid(i) = freq_range_min + &
                (freq_range_max - freq_range_min) * real(i - 1, dp) / real(n_alloc - 1, dp)
        end do

        sigma = smearing
        norm = 1.0_dp / (sigma * sqrt(2.0_dp * pi))

        do i = 1, n_alloc
            freq = data%freq_grid(i)
            do j = 1, data%n_modes
                ram_val = data%raman_activity(j)
                if (ram_val <= 0.0_dp) cycle
                dE = freq - data%freq(j)
                if (abs(dE) > GAUSS_CUTOFF * sigma) cycle
                gauss = norm * exp(-0.5_dp * (dE / sigma)**2)
                data%raman_spectrum(i) = data%raman_spectrum(i) + gauss * ram_val
            end do
        end do
    end subroutine compute_raman_spectrum


    subroutine free_vib_data(data)
        type(vib_data_t), intent(inout) :: data
        if (allocated(data%freq))           deallocate(data%freq)
        if (allocated(data%irrep))          deallocate(data%irrep)
        if (allocated(data%ir_intensity))   deallocate(data%ir_intensity)
        if (allocated(data%ir_active))      deallocate(data%ir_active)
        if (allocated(data%raman_activity)) deallocate(data%raman_activity)
        if (allocated(data%raman_active))   deallocate(data%raman_active)
        if (allocated(data%ir_spectrum))    deallocate(data%ir_spectrum)
        if (allocated(data%raman_spectrum)) deallocate(data%raman_spectrum)
        data%n_modes = 0
        data%has_raman_numeric = .false.
    end subroutine free_vib_data

end module castep_vib
