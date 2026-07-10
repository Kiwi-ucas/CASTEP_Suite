module phonon_dos
    !! Parse CASTEP .phonon files and compute phonon DOS via Gaussian smearing.
    !! Handles only frequencies, q-point weights, and eigenvectors.
    !! IR and Raman data are parsed from .castep via the castep_vib module.
    use castep_config, only: dp, pi
    implicit none
    private

    public :: phonon_dos_data_t, parse_phonon_file, compute_phonon_dos, &
        free_phonon_dos_data

    integer, parameter :: MAX_LINE = 1024

    type :: phonon_dos_data_t
        integer :: n_ions = 0, n_branches = 0, n_qpoints = 0
        real(dp) :: freq_min = 0.0_dp, freq_max = 0.0_dp
        real(dp) :: smearing = 1.0_dp       ! cm⁻¹
        real(dp) :: freq_grid(4001)          ! output frequency grid (cm⁻¹)
        real(dp), allocatable :: phdos(:)    ! PHDOS on grid
        real(dp), allocatable :: freqs(:)    ! raw frequencies (n_branches * n_qpoints)
        real(dp), allocatable :: weights(:)  ! q-point weights (n_qpoints)
    end type phonon_dos_data_t

contains

    subroutine parse_phonon_file(filename, data, iostat, iomsg)
        !! Parse .phonon file: header metadata, q-point frequencies + weights.
        !! Eigenvector blocks are skipped. No IR/Raman data is read —
        !! use the castep_vib module to parse those from .castep instead.
        character(len=*), intent(in) :: filename
        type(phonon_dos_data_t), intent(out) :: data
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg
        integer :: unit, ios, i, j, nb, nq
        character(len=MAX_LINE) :: line
        character(len=20) :: key
        real(dp) :: qx, qy, qz, weight, freq
        integer :: n_read

        iostat = 0

        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            iostat = 100
            if (present(iomsg)) iomsg = 'Cannot open file: ' // trim(filename)
            return
        end if

        ! ── Parse header ──
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'END header') > 0) exit
            if (index(line, 'Number of ions') > 0) then
                read(line, *, iostat=ios) key, key, key, data%n_ions
            else if (index(line, 'Number of branches') > 0) then
                read(line, *, iostat=ios) key, key, key, data%n_branches
            else if (index(line, 'Number of wavevectors') > 0) then
                read(line, *, iostat=ios) key, key, key, data%n_qpoints
            end if
        end do

        if (data%n_ions <= 0 .or. data%n_branches <= 0 .or. data%n_qpoints <= 0) then
            iostat = 101
            if (present(iomsg)) iomsg = 'Invalid header in .phonon file'
            close(unit); return
        end if

        nq = data%n_qpoints
        nb = data%n_branches
        allocate(data%freqs(nb * nq), data%weights(nq))

        ! ── Parse q-point blocks ──
        data%freq_min = huge(1.0_dp)
        data%freq_max = -huge(1.0_dp)

        do j = 1, nq
            ! Read q-point header: "q-pt= N qx qy qz weight"
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) then
                iostat = 102; if (present(iomsg)) iomsg = 'Unexpected end of .phonon file'
                close(unit); return
            end if
            read(line, *, iostat=ios) key, key, qx, qy, qz, weight
            if (ios /= 0) then
                iostat = 103; if (present(iomsg)) iomsg = 'Error reading q-point header'
                close(unit); return
            end if
            data%weights(j) = weight

            ! Read frequency lines: "N  freq" (2 fields, same for all q-points)
            do i = 1, nb
                read(unit, '(a)', iostat=ios) line
                if (ios /= 0) then
                    iostat = 104; if (present(iomsg)) iomsg = 'Error reading phonon frequencies'
                    close(unit); return
                end if
                read(line, *, iostat=ios) n_read, freq
                if (ios /= 0) then
                    iostat = 105; if (present(iomsg)) iomsg = 'Error parsing phonon frequency line'
                    close(unit); return
                end if
                data%freqs((j - 1) * nb + i) = freq
                if (freq < data%freq_min) data%freq_min = freq
                if (freq > data%freq_max) data%freq_max = freq
            end do

            ! Skip eigenvector block — read until next "q-pt=" line or EOF
            do
                read(unit, '(a)', iostat=ios) line
                if (ios /= 0) exit  ! EOF — last q-point
                if (index(line, 'q-pt=') > 0) then
                    backspace(unit)
                    exit
                end if
            end do
        end do

        close(unit)
    end subroutine parse_phonon_file


    subroutine compute_phonon_dos(data, freq_range_min, freq_range_max, n_points, &
                                   smearing, iostat, iomsg)
        !! Compute phonon DOS via Gaussian smearing of mode frequencies.
        type(phonon_dos_data_t), intent(inout) :: data
        real(dp), intent(in) :: freq_range_min, freq_range_max
        integer, intent(in) :: n_points
        real(dp), intent(in) :: smearing
        integer, intent(out) :: iostat
        character(len=*), optional, intent(out) :: iomsg

        real(dp) :: freq, sigma, norm, dE, gauss
        integer :: i, j, n_modes, n_alloc
        real(dp), parameter :: GAUSS_CUTOFF = 5.0_dp

        iostat = 0

        if (.not. allocated(data%freqs) .or. .not. allocated(data%weights)) then
            iostat = 102
            if (present(iomsg)) iomsg = 'No frequency data loaded'
            return
        end if

        n_modes = data%n_branches * data%n_qpoints
        n_alloc = min(n_points, 4001)

        if (allocated(data%phdos)) deallocate(data%phdos)
        allocate(data%phdos(n_alloc))
        data%phdos = 0.0_dp
        data%smearing = smearing

        ! Build frequency grid
        do i = 1, n_alloc
            data%freq_grid(i) = freq_range_min + &
                (freq_range_max - freq_range_min) * real(i - 1, dp) / real(n_alloc - 1, dp)
        end do

        sigma = smearing
        norm = 1.0_dp / (sigma * sqrt(2.0_dp * pi))

        ! Gaussian smearing with 5σ cutoff
        do i = 1, n_alloc
            freq = data%freq_grid(i)
            do j = 1, n_modes
                dE = freq - data%freqs(j)
                if (abs(dE) > GAUSS_CUTOFF * sigma) cycle
                gauss = norm * exp(-0.5_dp * (dE / sigma)**2)
                data%phdos(i) = data%phdos(i) + gauss * data%weights((j - 1) / data%n_branches + 1)
            end do
        end do
    end subroutine compute_phonon_dos


    subroutine free_phonon_dos_data(data)
        type(phonon_dos_data_t), intent(inout) :: data
        if (allocated(data%freqs))   deallocate(data%freqs)
        if (allocated(data%weights)) deallocate(data%weights)
        if (allocated(data%phdos))   deallocate(data%phdos)
        data%n_ions = 0; data%n_branches = 0; data%n_qpoints = 0
    end subroutine free_phonon_dos_data

end module phonon_dos
