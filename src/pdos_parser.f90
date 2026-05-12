module pdos_parser
    !! Parse CASTEP .pdos_weights / .pdos_bin binary files into pdos_data_t
    !! Binary format: big-endian, record-delimited
    !!   each record: [4B u32 size] [data] [4B u32 end_marker (=size)]
    !! .pdos_bin has 2 extra prefix records (f64 version + version string)
    use iso_fortran_env, only: int64
    use castep_config, only: dp, pdos_data_t, IO_SUCCESS, &
        IO_PDOS_NOT_FOUND, IO_PDOS_PARSE_ERROR
    implicit none
    private

    public :: parse_pdos_file, free_pdos_data

    integer, parameter :: MAX_PAYLOAD = 8192

contains

    ! ----------------------------------------------------------------
    !  Big-endian byte readers
    ! ----------------------------------------------------------------
    subroutine read_be_u32(unit, val, ios)
        integer, intent(in)  :: unit
        integer, intent(out) :: val
        integer, intent(out) :: ios
        character(1) :: b(4)
        integer :: i
        read(unit, iostat=ios) b
        if (ios /= 0) return
        val = 0
        do i = 1, 4
            val = ishft(val, 8) + iachar(b(i))
        end do
    end subroutine read_be_u32

    subroutine read_be_f64(unit, val, ios)
        integer, intent(in)  :: unit
        real(dp), intent(out) :: val
        integer, intent(out) :: ios
        integer(int64) :: raw
        character(1) :: b(8)
        integer :: i
        read(unit, iostat=ios) b
        if (ios /= 0) return
        raw = 0_int64
        do i = 1, 8
            raw = ishft(raw, 8) + int(iachar(b(i)), kind=int64)
        end do
        val = transfer(raw, val)
    end subroutine read_be_f64

    ! ----------------------------------------------------------------
    !  Record helpers (use local fixed-size buffer)
    ! ----------------------------------------------------------------
    subroutine read_payload(unit, nbytes, buf, ios)
        integer, intent(in) :: unit, nbytes
        character(len=*), intent(out) :: buf
        integer, intent(out) :: ios
        integer :: end_marker
        read(unit, iostat=ios) buf(1:nbytes)
        if (ios /= 0) return
        call read_be_u32(unit, end_marker, ios)
        if (ios /= 0) return
        if (end_marker /= nbytes) ios = 1
    end subroutine read_payload

    subroutine skip_record(unit, ios)
        integer, intent(in) :: unit
        integer, intent(out) :: ios
        integer :: n
        character(4096) :: tmp
        call read_be_u32(unit, n, ios)
        if (ios /= 0) return
        do while (n > 0)
            if (n <= 4096) then
                read(unit, iostat=ios) tmp(1:n)
                if (ios /= 0) return
                n = 0
            else
                read(unit, iostat=ios) tmp
                if (ios /= 0) return
                n = n - 4096
            end if
        end do
        ! skip end marker
        call read_be_u32(unit, n, ios)
    end subroutine skip_record

    ! ----------------------------------------------------------------
    !  Utility for u32/f64 from buffer
    ! ----------------------------------------------------------------
    function be_u32_at(buf, pos) result(v)
        character(len=*), intent(in) :: buf
        integer, intent(in) :: pos
        integer :: v, j
        v = 0
        do j = 0, 3
            v = ishft(v, 8) + iachar(buf(pos+j:pos+j))
        end do
    end function be_u32_at

    function be_f64_at(buf, pos) result(v)
        character(len=*), intent(in) :: buf
        integer, intent(in) :: pos
        real(dp) :: v
        integer(int64) :: raw
        integer :: j
        raw = 0_int64
        do j = 0, 7
            raw = ishft(raw, 8) + int(iachar(buf(pos+j:pos+j)), kind=int64)
        end do
        v = transfer(raw, v)
    end function be_f64_at

    ! ----------------------------------------------------------------
    !  Main parser — all payloads use local fixed-size buffer
    ! ----------------------------------------------------------------
    subroutine parse_pdos_file(filename, pdos, iostat, iomsg)
        character(len=*), intent(in)  :: filename
        type(pdos_data_t), intent(out) :: pdos
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        integer  :: unit, ios, ik, isp, ib, j
        integer  :: total_kp, nspin, norbs, maxb, rec_size, nbands_occ
        logical  :: exist, is_pdos_bin
        real(dp) :: rval
        integer  :: end_check
        character(MAX_PAYLOAD) :: lbuf   ! local fixed-size buffer for all payloads

        iostat = IO_SUCCESS
        if (present(iomsg)) iomsg = ''

        inquire(file=trim(filename), exist=exist)
        if (.not. exist) then
            iostat = IO_PDOS_NOT_FOUND
            if (present(iomsg)) iomsg = 'File not found: ' // trim(filename)
            return
        end if

        open(newunit=unit, file=trim(filename), access='stream', &
             form='unformatted', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            iostat = IO_PDOS_NOT_FOUND
            if (present(iomsg)) iomsg = 'Cannot open: ' // trim(filename)
            return
        end if

        ! --- detect .pdos_bin vs .pdos_weights ---
        is_pdos_bin = .false.
        call read_be_u32(unit, rec_size, ios)
        if (ios /= 0) then
            iostat = IO_PDOS_PARSE_ERROR
            if (present(iomsg)) iomsg = 'Empty or truncated file'
            close(unit); return
        end if

        if (rec_size == 8) then
            call read_be_f64(unit, rval, ios)
            if (ios == 0) then
                call read_be_u32(unit, end_check, ios)
                if (ios == 0 .and. end_check == 8) then
                    is_pdos_bin = .true.
                    ! skip version string record
                    call skip_record(unit, ios)
                    if (ios /= 0) then
                        iostat = IO_PDOS_PARSE_ERROR; close(unit); return
                    end if
                end if
            end if
        end if

        if (.not. is_pdos_bin) rewind(unit)

        ! --- parse header ---
        call read_u32_record(unit, total_kp, ios)
        if (ios /= 0) goto 999
        if (total_kp <= 0) goto 999

        call read_u32_record(unit, nspin, ios)
        if (ios /= 0) goto 999
        if (nspin < 1 .or. nspin > 2) goto 999

        call read_u32_record(unit, norbs, ios)
        if (ios /= 0) goto 999
        if (norbs <= 0) goto 999

        call read_u32_record(unit, maxb, ios)
        if (ios /= 0) goto 999
        if (maxb <= 0) goto 999

        pdos%total_kpoints = total_kp
        pdos%num_spins     = nspin
        pdos%num_orbitals  = norbs
        pdos%max_bands     = maxb

        allocate(pdos%orbital_species(norbs), pdos%orbital_ion(norbs), &
                 pdos%orbital_am(norbs), pdos%kpoint_coords(3, total_kp), &
                 pdos%kpoint_indices(total_kp), &
                 pdos%orbital_weights(norbs, maxb, total_kp, nspin), stat=ios)
        if (ios /= 0) goto 999
        pdos%orbital_weights = 0.0_dp

        ! orbital_species
        call read_int_vec(unit, 4*norbs, pdos%orbital_species, norbs, ios)
        if (ios /= 0) goto 999

        ! orbital_ion
        call read_int_vec(unit, 4*norbs, pdos%orbital_ion, norbs, ios)
        if (ios /= 0) goto 999

        ! orbital_am
        call read_int_vec(unit, 4*norbs, pdos%orbital_am, norbs, ios)
        if (ios /= 0) goto 999

        ! --- parse k-point blocks ---
        do ik = 1, total_kp
            ! k-point header: 28 bytes = index + 3*f64
            call read_be_u32(unit, rec_size, ios)
            if (ios /= 0 .or. rec_size /= 28) goto 999
            call read_payload(unit, 28, lbuf, ios)
            if (ios /= 0) goto 999
            pdos%kpoint_indices(ik)   = be_u32_at(lbuf, 1)
            pdos%kpoint_coords(1, ik) = be_f64_at(lbuf, 5)
            pdos%kpoint_coords(2, ik) = be_f64_at(lbuf, 13)
            pdos%kpoint_coords(3, ik) = be_f64_at(lbuf, 21)

            do isp = 1, nspin
                ! spin index
                call read_u32_record(unit, rec_size, ios)
                if (ios /= 0) goto 999

                ! nbands_occ
                call read_u32_record(unit, nbands_occ, ios)
                if (ios /= 0) goto 999
                if (nbands_occ > maxb) nbands_occ = maxb

                ! orbital weights per band
                do ib = 1, nbands_occ
                    call read_be_u32(unit, rec_size, ios)
                    if (ios /= 0 .or. rec_size /= 8 * norbs) goto 999
                    if (rec_size > MAX_PAYLOAD) goto 999
                    call read_payload(unit, rec_size, lbuf, ios)
                    if (ios /= 0) goto 999
                    do j = 1, norbs
                        pdos%orbital_weights(j, ib, ik, isp) = &
                            be_f64_at(lbuf, 8*(j-1)+1)
                    end do
                end do
            end do
        end do

        close(unit)
        return

999     continue
        iostat = IO_PDOS_PARSE_ERROR
        if (present(iomsg)) iomsg = 'Failed to parse PDOS file: ' // trim(filename)
        close(unit)
        call free_pdos_data(pdos)
    end subroutine parse_pdos_file

    ! ----------------------------------------------------------------
    !  Helper: read a single u32 record
    ! ----------------------------------------------------------------
    subroutine read_u32_record(unit, val, ios)
        integer, intent(in)  :: unit
        integer, intent(out) :: val
        integer, intent(out) :: ios
        integer               :: rec_size
        character(4) :: buf4
        call read_be_u32(unit, rec_size, ios)
        if (ios /= 0) return
        if (rec_size /= 4) then; ios = 1; return; end if
        call read_payload(unit, 4, buf4, ios)
        if (ios /= 0) return
        val = be_u32_at(buf4, 1)
    end subroutine read_u32_record

    ! ----------------------------------------------------------------
    !  Helper: read vector of u32
    ! ----------------------------------------------------------------
    subroutine read_int_vec(unit, expect_bytes, arr, n, ios)
        integer, intent(in)  :: unit, expect_bytes, n
        integer, intent(out) :: arr(n)
        integer, intent(out) :: ios
        integer :: rec_size, j
        character(MAX_PAYLOAD) :: tbuf
        call read_be_u32(unit, rec_size, ios)
        if (ios /= 0) return
        if (rec_size /= expect_bytes .or. rec_size > MAX_PAYLOAD) then
            ios = 1; return
        end if
        call read_payload(unit, rec_size, tbuf, ios)
        if (ios /= 0) return
        do j = 1, n
            arr(j) = be_u32_at(tbuf, 4*(j-1)+1)
        end do
    end subroutine read_int_vec

    ! ----------------------------------------------------------------
    !  Cleanup
    ! ----------------------------------------------------------------
    subroutine free_pdos_data(pdos)
        type(pdos_data_t), intent(inout) :: pdos
        if (allocated(pdos%orbital_species))  deallocate(pdos%orbital_species)
        if (allocated(pdos%orbital_ion))      deallocate(pdos%orbital_ion)
        if (allocated(pdos%orbital_am))       deallocate(pdos%orbital_am)
        if (allocated(pdos%kpoint_coords))    deallocate(pdos%kpoint_coords)
        if (allocated(pdos%kpoint_indices))   deallocate(pdos%kpoint_indices)
        if (allocated(pdos%orbital_weights))  deallocate(pdos%orbital_weights)
        pdos%total_kpoints = 0
        pdos%num_spins     = 1
        pdos%num_orbitals  = 0
        pdos%max_bands     = 0
    end subroutine free_pdos_data

end module pdos_parser
