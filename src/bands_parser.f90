module bands_parser
    !! Parse CASTEP .bands output files into structured bands_data_t
    use castep_config, only: dp, bands_data_t, MAX_LINE_LEN, &
        IO_SUCCESS, IO_BANDS_NOT_FOUND, IO_BANDS_PARSE_ERROR
    implicit none
    private

    public :: parse_bands_file, free_bands_data

contains

    subroutine parse_bands_file(filename, bands, iostat, iomsg)
        character(len=*), intent(in)  :: filename
        type(bands_data_t), intent(out) :: bands
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg
        character(len=MAX_LINE_LEN) :: line, msg
        integer  :: unit, ios, ik, is, ie, ne_read
        logical  :: exist

        iostat = IO_SUCCESS
        if (present(iomsg)) iomsg = ''

        inquire(file=trim(filename), exist=exist)
        if (.not. exist) then
            iostat = IO_BANDS_NOT_FOUND
            if (present(iomsg)) iomsg = 'File not found: ' // trim(filename)
            return
        end if

        open(newunit=unit, file=trim(filename), status='old', &
             action='read', iostat=ios)
        if (ios /= 0) then
            iostat = IO_BANDS_NOT_FOUND
            if (present(iomsg)) iomsg = 'Cannot open: ' // trim(filename)
            return
        end if

        ! --- header ---
        call parse_header_int(unit, 'Number of k-points', bands%num_kpoints, msg)
        if (bands%num_kpoints <= 0) then
            iostat = IO_BANDS_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Bad num_kpoints: ' // trim(msg)
            close(unit); return
        end if

        call parse_header_int(unit, 'Number of spin components', bands%num_spin, msg)
        if (bands%num_spin < 1) bands%num_spin = 1

        call parse_header_real(unit, 'Number of electrons', bands%num_electrons, msg)

        call parse_header_int(unit, 'Number of eigenvalues', bands%num_eigenvalues, msg)
        if (bands%num_eigenvalues <= 0) then
            iostat = IO_BANDS_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Bad num_eigenvalues: ' // trim(msg)
            close(unit); return
        end if

        call parse_header_real(unit, 'Fermi energy', bands%fermi_energy, msg)

        ! skip "Unit cell vectors" label line, then read 3 lines of vectors
        read(unit, '(a)', iostat=ios) line
        if (ios /= 0) then
            close(unit); return
        end if
        bands%cell_vectors = 0.0_dp
        do ik = 1, 3
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            read(line, *, iostat=ios) bands%cell_vectors(1:3, ik)
            if (ios /= 0) bands%cell_vectors(1:3, ik) = 0.0_dp
        end do

        ! --- allocate arrays ---
        allocate(bands%kpoint_indices(bands%num_kpoints), &
                 bands%kpoint_coords(4, bands%num_kpoints), &
                 bands%kpath_dist(bands%num_kpoints), &
                 bands%eigenvalues(bands%num_eigenvalues, &
                                  bands%num_kpoints, bands%num_spin), &
                 stat=ios)
        if (ios /= 0) then
            iostat = IO_BANDS_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Allocation failed'
            close(unit); return
        end if
        bands%kpoint_indices = 0
        bands%kpoint_coords  = 0.0_dp
        bands%kpath_dist     = 0.0_dp
        bands%eigenvalues    = 0.0_dp

        ! --- read k-point blocks ---
        do ik = 1, bands%num_kpoints
            ! read "K-point" line
            call read_next_nonempty(unit, line, ios)
            if (ios /= 0) then
                iostat = IO_BANDS_PARSE_ERROR
                if (present(iomsg)) iomsg = 'Unexpected EOF at k-point'
                close(unit); return
            end if
            read(line, *, iostat=ios) bands%kpoint_indices(ik), &
                bands%kpoint_coords(1,ik), bands%kpoint_coords(2,ik), &
                bands%kpoint_coords(3,ik), bands%kpoint_coords(4,ik)
            if (ios /= 0) then
                ! try extracting numeric tokens after "K-point"
                call parse_kpoint_line(line, bands%kpoint_indices(ik), &
                    bands%kpoint_coords(1,ik), bands%kpoint_coords(2,ik), &
                    bands%kpoint_coords(3,ik), bands%kpoint_coords(4,ik), ios)
                if (ios /= 0) then
                    iostat = IO_BANDS_PARSE_ERROR
                    if (present(iomsg)) iomsg = 'Bad K-point line'
                    close(unit); return
                end if
            end if

            do is = 1, bands%num_spin
                ! read "Spin component" line
                call read_next_nonempty(unit, line, ios)
                if (ios /= 0) then
                    iostat = IO_BANDS_PARSE_ERROR
                    if (present(iomsg)) iomsg = 'Missing spin component line'
                    close(unit); return
                end if
                ! skip spin label, just verify it's there
                if (index(line, 'Spin component') == 0) then
                    iostat = IO_BANDS_PARSE_ERROR
                    if (present(iomsg)) iomsg = 'Expected Spin component line'
                    close(unit); return
                end if

                ! read eigenvalues
                ne_read = 0
                ie = 0
                do while (ie < bands%num_eigenvalues)
                    read(unit, '(a)', iostat=ios) line
                    if (ios /= 0) exit
                    if (len_trim(line) == 0) cycle
                    ie = ie + 1
                    read(line, *, iostat=ios) bands%eigenvalues(ie, ik, is)
                    if (ios /= 0) exit
                    ne_read = ne_read + 1
                end do
                if (ne_read < bands%num_eigenvalues) then
                    write(*, '(a,i0,a,i0)') '  Warning: expected ', &
                        bands%num_eigenvalues, ' eigenvalues, got ', ne_read
                    bands%num_eigenvalues = ne_read
                end if
            end do
        end do

        close(unit)

        ! --- compute cumulative k-path distance ---
        call build_kpath_dist(bands%num_kpoints, bands%kpoint_coords, &
                              bands%kpath_dist)
    end subroutine parse_bands_file


    subroutine parse_header_int(unit, key, val, msg)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: key
        integer, intent(out) :: val
        character(len=*), intent(out) :: msg
        character(len=MAX_LINE_LEN) :: line
        integer :: ios, pos, nlen
        val = 0; msg = ''
        read(unit, '(a)', iostat=ios) line
        if (ios /= 0) then
            msg = 'EOF reading ' // key; return
        end if
        nlen = len_trim(key)
        pos = index(line, trim(key))
        if (pos == 0) then
            msg = 'Key not found: ' // key; return
        end if
        read(line(pos+nlen:), *, iostat=ios) val
        if (ios /= 0) msg = 'Bad integer for ' // key
    end subroutine parse_header_int


    subroutine parse_header_real(unit, key, val, msg)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: key
        real(dp), intent(out) :: val
        character(len=*), intent(out) :: msg
        character(len=MAX_LINE_LEN) :: line, suffix
        integer :: ios, pos, nlen
        val = 0.0_dp; msg = ''
        read(unit, '(a)', iostat=ios) line
        if (ios /= 0) then
            msg = 'EOF reading ' // key; return
        end if
        nlen = len_trim(key)
        pos = index(line, trim(key))
        if (pos == 0) then
            msg = 'Key not found: ' // key; return
        end if
        suffix = adjustl(line(pos+nlen:))
        ! strip parenthesized text like "(in atomic units)"
        call strip_parens(suffix)
        read(suffix, *, iostat=ios) val
        if (ios /= 0) msg = 'Bad real for ' // key
    end subroutine parse_header_real

    subroutine strip_parens(s)
        character(len=*), intent(inout) :: s
        integer :: i, j, depth
        character(len=MAX_LINE_LEN) :: tmp
        tmp = ''
        j = 1; depth = 0
        do i = 1, len_trim(s)
            if (s(i:i) == '(') then
                depth = depth + 1
            else if (s(i:i) == ')') then
                depth = depth - 1
            else if (depth == 0) then
                tmp(j:j) = s(i:i)
                j = j + 1
            end if
        end do
        s = adjustl(tmp)
    end subroutine strip_parens


    subroutine parse_kpoint_line(line, idx, kx, ky, kz, wt, ok)
        character(len=*), intent(in)  :: line
        integer,  intent(out) :: idx, ok
        real(dp), intent(out) :: kx, ky, kz, wt
        character(len=MAX_LINE_LEN) :: tmp
        integer :: ios, pos
        pos = index(line, 'K-point')
        if (pos == 0) then
            ok = 1; return
        end if
        tmp = adjustl(line(pos+7:))
        read(tmp, *, iostat=ios) idx, kx, ky, kz, wt
        if (ios /= 0) then
            ok = 1; return
        end if
        ok = 0
    end subroutine parse_kpoint_line


    subroutine build_kpath_dist(nkp, coords, dist)
        integer,  intent(in)  :: nkp
        real(dp), intent(in)  :: coords(4, nkp)
        real(dp), intent(out) :: dist(nkp)
        integer  :: i
        real(dp) :: dx, dy, dz
        if (nkp < 1) return
        dist(1) = 0.0_dp
        do i = 2, nkp
            dx = coords(1,i) - coords(1,i-1)
            dy = coords(2,i) - coords(2,i-1)
            dz = coords(3,i) - coords(3,i-1)
            dist(i) = dist(i-1) + sqrt(dx*dx + dy*dy + dz*dz)
        end do
    end subroutine build_kpath_dist


    subroutine read_next_nonempty(unit, line, ios)
        integer, intent(in) :: unit
        character(len=*), intent(out) :: line
        integer, intent(out) :: ios
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) return
            if (len_trim(line) > 0) return
        end do
    end subroutine read_next_nonempty


    subroutine free_bands_data(bands)
        type(bands_data_t), intent(inout) :: bands
        if (allocated(bands%kpoint_indices)) deallocate(bands%kpoint_indices)
        if (allocated(bands%kpoint_coords))  deallocate(bands%kpoint_coords)
        if (allocated(bands%kpath_dist))     deallocate(bands%kpath_dist)
        if (allocated(bands%eigenvalues))    deallocate(bands%eigenvalues)
        bands%num_kpoints     = 0
        bands%num_spin        = 1
        bands%num_eigenvalues = 0
        bands%num_electrons   = 0.0_dp
        bands%fermi_energy    = 0.0_dp
        bands%cell_vectors    = 0.0_dp
    end subroutine free_bands_data

end module bands_parser
